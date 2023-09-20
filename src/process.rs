use crate::genotype::{MutationalBurden, NeutralMutationPoisson, Sfs};
use crate::stemcell::mutate;
use crate::subclone::{
    assign, proliferating_cell, save_variant_fraction, CloneId, Distributions, SubClones, Variants,
};
use crate::MAX_SUBCLONES;
use anyhow::Context;
use rand::Rng;
use rand_distr::{Distribution, WeightedIndex};
use sosa::{AdvanceStep, CurrentState, NextReaction};
use std::collections::VecDeque;
use std::fs;
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct ProcessOptions {
    pub distributions: Distributions,
    pub path: PathBuf,
    pub cells2subsample: Option<Vec<usize>>,
    pub neutral_poisson: NeutralMutationPoisson,
}

#[derive(Hash, PartialEq, Eq)]
pub enum Stats2Save {
    Burden,
    BurdenEntropy,
    Genotypes,
    Sfs,
    SfsEntropy,
    Stats,
    VariantFraction,
}

#[derive(Debug, Clone)]
pub struct CellDivisionProbabilities {
    pub p_asymmetric: f64,
    /// Rate of neutral mutations per cell-division
    pub lambda_poisson: f32,
    /// Probability of getting a fit mutant upon cell division
    pub p: f64,
}

/// Exponential growing process
#[derive(Clone, Debug)]
pub struct Exponential {
    /// A collection of clones having a proliferative advantage.
    pub subclones: SubClones,
    /// The counter for the number of proliferative events.
    pub counter_divisions: usize,
    id: usize,
    pub verbosity: u8,
    pub distributions: Distributions,
    pub neutral_mutations: NeutralMutationPoisson,
}

impl Exponential {
    pub fn new(
        process_options: ProcessOptions,
        initial_subclones: SubClones,
        id: usize,
        verbosity: u8,
    ) -> Exponential {
        let hsc = Exponential {
            subclones: initial_subclones,
            distributions: process_options.distributions,
            counter_divisions: 0,
            id,
            verbosity,
            neutral_mutations: process_options.neutral_poisson,
        };
        if verbosity > 1 {
            println!("process created: {:#?}", hsc);
        }
        hsc
    }

    pub fn switch_to_moran(
        self,
        process_options: ProcessOptions,
        snapshot: VecDeque<f32>,
    ) -> Moran {
        Moran {
            subclones: self.subclones,
            counter_divisions: self.counter_divisions,
            id: self.id,
            time: 0.,
            cells2subsample: process_options.cells2subsample,
            snapshot,
            path2dir: process_options.path,
            verbosity: self.verbosity,
            distributions: process_options.distributions,
            neutral_mutations: process_options.neutral_poisson,
        }
    }
}

impl AdvanceStep<MAX_SUBCLONES> for Exponential {
    /// The id of a subclone [`CloneId`].
    type Reaction = CloneId;

    fn advance_step(&mut self, reaction: NextReaction<Self::Reaction>, rng: &mut impl Rng) {
        //! Update the process by simulating the next proliferative event
        //! according to the next `reaction` determined by the [Gillespie
        //! algorithm](`sosa`).
        //!
        //! The proliferation step is implemented as following:
        //!
        //! 1. select the cell that will proliferate next from the clone
        //! with id `reaction` determined by the Gillespie algorithm
        //!
        //! 2. for the proliferating cell:
        //!     * mutate genome by storing a neutral number of TODO
        //!
        //!     * assign to new subclone with a probability determined by the
        //!     rate of mutations conferring a proliferative advantage
        // The Gillespie sampler samples the clone that will proliferate next,
        // that is the clone with id `reaction.event`.
        // Pick random proliferating cells from this clone. **Note that this
        // removes the cells from the clone with id `reaction.event`.**
        let mut stem_cell =
            proliferating_cell(&mut self.subclones, reaction.event, self.verbosity, rng);
        self.counter_divisions += 1;
        self.subclones
            .get_mut_clone_unchecked(reaction.event)
            .assign_cell(stem_cell.clone());

        if self.verbosity > 2 {
            println!("assigning mutations to cell {:#?}", stem_cell)
        }

        let division = self.neutral_mutations.new_muts_upon_division(rng);
        if let Some(mutations) = division {
            mutate(&mut stem_cell, mutations);
        }

        // assign cell to clone. This can have two outcomes based on a
        // Bernouilli trial see `assign` and `self.assign`:
        // 1. the stem cell is re-assigned to the old clone with id
        // `reaction.event` (remember that `self.proliferating_cells` has
        // removed the cells from the clone)
        // 2. the stem cell is assigned to a new clone with id different
        // from `reaction.event`.
        assign(
            &mut self.subclones,
            reaction.event,
            stem_cell,
            &self.distributions,
            rng,
            self.verbosity,
        );
    }

    fn update_state(&self, state: &mut CurrentState<MAX_SUBCLONES>) {
        state.population = Variants::variant_counts(&self.subclones);
    }
}

/// The Moran process saves the state of the agents and simulates new
/// proliferative events at each timestep according to the Gillespie algorithm,
/// see [`Moran::advance_step`].
#[derive(Debug, Clone)]
pub struct Moran {
    /// A collection of clones having a proliferative advantage.
    pub subclones: SubClones,
    /// The counter for the number of proliferative events.
    pub counter_divisions: usize,
    id: usize,
    time: f32,
    pub snapshot: VecDeque<f32>,
    pub path2dir: PathBuf,
    pub verbosity: u8,
    pub distributions: Distributions,
    pub neutral_mutations: NeutralMutationPoisson,
    pub cells2subsample: Option<Vec<usize>>,
}

impl Default for Moran {
    fn default() -> Self {
        let process_options = ProcessOptions {
            distributions: Distributions::default(),
            path: PathBuf::from("./output"),
            cells2subsample: None,
            neutral_poisson: NeutralMutationPoisson::new(10., 1.).unwrap(),
        };

        Moran::new(
            process_options,
            SubClones::default(),
            vec![0.01, 0.1],
            1,
            0.,
            1,
        )
    }
}

impl Moran {
    /// A Moran process with wild-type subclone with neutral fitness, to which
    /// all cells will be assigned.
    pub fn new(
        process_options: ProcessOptions,
        initial_subclones: SubClones,
        mut snapshot: Vec<f32>,
        id: usize,
        time: f32,
        verbosity: u8,
    ) -> Moran {
        snapshot.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let snapshot = VecDeque::from(snapshot);
        let hsc = Moran {
            subclones: initial_subclones,
            distributions: process_options.distributions,
            counter_divisions: 0,
            id,
            path2dir: process_options.path,
            time,
            snapshot,
            cells2subsample: process_options.cells2subsample,
            verbosity,
            neutral_mutations: process_options.neutral_poisson,
        };
        if verbosity > 1 {
            println!("process created: {:#?}", hsc);
        }
        hsc
    }

    fn keep_const_population_upon_symmetric_division(
        &mut self,
        subclone_id: usize,
        rng: &mut impl Rng,
    ) {
        //! If an symmetric division is performed, need to remove a random cell
        //! from the population.
        //!
        //! We proceed as following:
        //!     1. check if there is only one clone in the population, then
        //!     remove a cell from this clone
        //!     2. else, compute the variant counts (removing one cell to avoid
        //!     sampling the cell that will proliferate)
        //!     3. and sample from any clone based on the weights defined by
        //!     the variant counts
        // remove a cell from a random subclone based on the frequencies of
        // the clones at the current state
        let id2remove = if let Some(id) = self.subclones.the_only_one_subclone_present() {
            id
        } else {
            let mut variants = Variants::variant_counts(&self.subclones);
            // do not sample the cell that will proliferate
            variants[subclone_id] -= 1;
            WeightedIndex::new(variants).unwrap().sample(rng)
        };
        // remove a cell from the random subclone
        let _ = self
            .subclones
            .get_mut_clone_unchecked(id2remove)
            .random_cell(rng)
            .with_context(|| "found empty subclone")
            .unwrap();
        if self.verbosity > 2 {
            println!("removing one cell from clone {}", id2remove);
        }
    }

    pub fn make_path(
        &self,
        tosave: Stats2Save,
        cells: usize,
        timepoint: usize,
    ) -> anyhow::Result<PathBuf> {
        let path2dir = self.path2dir.join(format!("{}cells", cells));
        let path2file = match tosave {
            Stats2Save::VariantFraction => path2dir.join("variant_fraction"),
            Stats2Save::Genotypes => path2dir.join("genotypes"),
            Stats2Save::Burden => path2dir.join("burden"),
            Stats2Save::BurdenEntropy => path2dir.join("burden_entropy"),
            Stats2Save::Sfs => path2dir.join("sfs"),
            Stats2Save::SfsEntropy => path2dir.join("sfs_entropy"),
            Stats2Save::Stats => path2dir.join("stats"),
        };
        let path2file = path2file.join(timepoint.to_string());
        fs::create_dir_all(&path2file).with_context(|| "Cannot create dir")?;
        if self.verbosity > 1 {
            println!("creating dirs {:#?}", path2file);
        }
        Ok(path2file.join(self.id.to_string()))
    }

    pub fn save(
        &self,
        timepoint: usize,
        nb_cells: usize,
        rng: &mut impl Rng,
    ) -> anyhow::Result<()> {
        save_variant_fraction(
            &self.subclones,
            &self.make_path(Stats2Save::VariantFraction, nb_cells, timepoint)?,
            self.verbosity,
        )?;

        let cells = if nb_cells == self.subclones.compute_tot_cells() as usize {
            self.subclones.get_cells()
        } else {
            self.subclones.get_cells_subsampled(nb_cells, rng)
        };

        if self.verbosity > 1 {
            println!("saving {} cells", cells.len());
        }

        Sfs::from_cells(&cells, self.verbosity)
            .unwrap_or_else(|_| panic!("cannot create SFS for timepoint {}", timepoint))
            .save(&self.make_path(Stats2Save::Sfs, nb_cells, timepoint)?)?;

        MutationalBurden::from_cells(&cells, self.verbosity)
            .unwrap_or_else(|_| panic!("cannot create burden for the timepoint {}", timepoint))
            .save(&self.make_path(Stats2Save::Burden, nb_cells, timepoint)?)?;

        if self.verbosity > 0 {
            println!(
                "saving measurements after {} mutational events",
                self.counter_divisions
            );
        }

        Ok(())
    }
}

impl AdvanceStep<MAX_SUBCLONES> for Moran {
    /// The id of a subclone [`CloneId`].
    type Reaction = CloneId;

    fn advance_step(&mut self, reaction: NextReaction<Self::Reaction>, rng: &mut impl Rng) {
        //! Update the process by simulating the next proliferative event
        //! according to the next `reaction` determined by the [Gillespie
        //! algorithm](`sosa`).
        //!
        //! The proliferation step is implemented as following:
        //!
        //! 1. select the cell that will proliferate next from the clone
        //! with id `reaction` determined by the Gillespie algorithm
        //!
        //! 2. if the next proliferation is performed asymmetricaly, continue
        //! to step 3. Else, clone the proliferating cell and assign it to the
        //! subclone with id `reaction` (from step 1); next randomly select
        //! another cell from another clone to keep to population constant
        //!
        //! 3. for the proliferating cell only (not the clone from 2):
        //!     * mutate genome by TODO
        //!
        //!     * assign to new subclone with a probability determined by the
        //!     rate of mutations conferring a proliferative advantage
        // The Gillespie sampler samples the clone that will proliferate next,
        // that is the clone with id `reaction.event`.
        // Pick random proliferating cells from this clone. **Note that this
        // removes the cells from the clone with id `reaction.event`.**
        self.time += reaction.time;
        let mut stem_cell =
            proliferating_cell(&mut self.subclones, reaction.event, self.verbosity, rng);
        self.counter_divisions += 1;
        let interdivison_time = self.time - stem_cell.last_division_t;
        if self.verbosity > 1 {
            println!("cell {:#?} is dividing", stem_cell);
            println!("at time {}", self.time);
        }

        let background =
            self.neutral_mutations
                .new_muts_background(interdivison_time, rng, self.verbosity);
        if self.verbosity > 2 {
            println!(
                "assigning {} background mutations to cell {:#?}",
                background.as_ref().unwrap_or(&vec![]).len(),
                stem_cell
            )
        }
        if let Some(mutations) = background {
            mutate(&mut stem_cell, mutations);
        }
        stem_cell.last_division_t = self.time;

        if self.distributions.can_only_be_symmetric()
            || !self.distributions.bern_asymmetric.sample(rng)
        {
            // remove a cell from the population
            if self.verbosity > 2 {
                println!("keeping the cell population constant");
            }
            self.keep_const_population_upon_symmetric_division(reaction.event, rng);
            self.subclones
                .get_mut_clone_unchecked(reaction.event)
                .assign_cell(stem_cell.clone());
        }

        let division = self.neutral_mutations.new_muts_upon_division(rng);
        if self.verbosity > 2 {
            println!(
                "assigning {} mutations upon cell division to cell {:#?}",
                division.as_ref().unwrap_or(&vec![]).len(),
                stem_cell
            )
        }
        if let Some(mutations) = division {
            mutate(&mut stem_cell, mutations);
        }

        // assign cell to clone. This can have two outcomes based on a
        // Bernouilli trial see `assign` and `self.assign`:
        // 1. the stem cell is re-assigned to the old clone with id
        // `reaction.event` (remember that `self.proliferating_cells` has
        // removed the cells from the clone)
        // 2. the stem cell is assigned to a new clone with id different
        // from `reaction.event`.
        assign(
            &mut self.subclones,
            reaction.event,
            stem_cell,
            &self.distributions,
            rng,
            self.verbosity,
        );

        // take snapshot
        if self.verbosity > 2 {
            println!(
                "{:#?} timepoints to save, time now is {}",
                self.snapshot, self.time
            );
        }
        if let Some(&time) = self.snapshot.front() {
            if self.time >= time {
                if self.verbosity > 0 {
                    println!(
                        "saving variant fraction for {:#?} at time {} for timepoint {}",
                        time,
                        self.time,
                        self.snapshot.len()
                    );
                }
                let mut cells2save = vec![self.subclones.compute_tot_cells() as usize];
                if let Some(subsampling) = self.cells2subsample.as_ref() {
                    for cell in subsampling {
                        cells2save.push(*cell);
                    }
                }
                for cells in cells2save {
                    self.save(self.snapshot.len(), cells, rng)
                        .expect("cannot save snapshot");
                }
                self.snapshot.pop_front();
            }
        }
    }

    fn update_state(&self, state: &mut CurrentState<MAX_SUBCLONES>) {
        state.population = Variants::variant_counts(&self.subclones);
    }
}
