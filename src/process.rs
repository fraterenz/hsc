use crate::genotype::{MutationalBurden, Sfs, Variant};
use crate::stemcell::{assign_background_mutations, mutate, StemCell};
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
pub struct Snapshot {
    /// The number of cells to subsample
    pub cells2sample: usize,
    /// The time at which we subsample
    pub time: f32,
}

#[derive(Debug, Clone)]
pub struct SavingOptions {
    pub filename: PathBuf,
    pub save_sfs_only: bool,
}

#[derive(Debug, Clone)]
pub struct ProcessOptions {
    pub path: PathBuf,
    pub snapshots: VecDeque<Snapshot>,
}

#[derive(Hash, PartialEq, Eq)]
pub enum Stats2Save {
    Burden,
    Sfs,
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
    pub verbosity: u8,
    pub distributions: Distributions,
}

impl Exponential {
    pub fn new(
        initial_subclones: SubClones,
        distributions: Distributions,
        verbosity: u8,
    ) -> Exponential {
        let hsc = Exponential {
            subclones: initial_subclones,
            distributions,
            counter_divisions: 0,
            verbosity,
        };
        if verbosity > 1 {
            println!("process created: {:#?}", hsc);
        }
        hsc
    }

    pub fn switch_to_moran(
        self,
        process_options: ProcessOptions,
        distributions: Distributions,
        filename: PathBuf,
        save_sfs_only: bool,
    ) -> Moran {
        Moran {
            subclones: self.subclones,
            counter_divisions: self.counter_divisions,
            time: 0.,
            snapshots: process_options.snapshots,
            path2dir: process_options.path,
            verbosity: self.verbosity,
            filename,
            distributions,
            save_sfs_only,
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
        //! 2. for the proliferating cell and its daughter cell:
        //!     * mutate genome by storing a neutral number of mutations
        //!
        //!     * assign to new subclone with a probability determined by the
        //!     rate of mutations conferring a proliferative advantage
        // The Gillespie sampler samples the clone that will proliferate next,
        // that is the clone with id `reaction.event`.
        // Pick random proliferating cells from this clone. **Note that this
        // removes the cells from the clone with id `reaction.event`.**
        let stem_cell =
            proliferating_cell(&mut self.subclones, reaction.event, self.verbosity, rng);
        self.counter_divisions += 1;
        for mut cell in [stem_cell.clone(), stem_cell] {
            if self.verbosity > 2 {
                println!("assigning mutations to cell {:#?}", cell)
            }

            let division = self
                .distributions
                .neutral_poisson
                .new_muts_upon_division(rng);
            if let Some(mutations) = division {
                mutate(&mut cell, mutations);
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
                cell,
                &self.distributions,
                rng,
                self.verbosity,
            );
        }
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
    pub time: f32,
    pub path2dir: PathBuf,
    pub verbosity: u8,
    pub distributions: Distributions,
    pub snapshots: VecDeque<Snapshot>,
    pub filename: PathBuf,
    pub save_sfs_only: bool,
}

impl Default for Moran {
    fn default() -> Self {
        let process_options = ProcessOptions {
            path: PathBuf::from("./output"),
            snapshots: VecDeque::default(),
        };

        Moran::new(
            process_options,
            SubClones::default(),
            0.,
            SavingOptions {
                filename: PathBuf::default(),
                save_sfs_only: false,
            },
            Distributions::default(),
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
        time: f32,
        saving_options: SavingOptions,
        distributions: Distributions,
        verbosity: u8,
    ) -> Moran {
        let hsc = Moran {
            subclones: initial_subclones,
            distributions,
            counter_divisions: 0,
            path2dir: process_options.path,
            time,
            snapshots: process_options.snapshots,
            filename: saving_options.filename,
            verbosity,
            save_sfs_only: saving_options.save_sfs_only,
        };
        if verbosity > 1 {
            println!("process created: {:#?}", hsc);
        }
        hsc
    }

    fn assign_mutations(&self, stem_cell: &mut StemCell, mutations: Option<Vec<Variant>>) {
        if let Some(mutations) = mutations {
            if self.verbosity > 2 {
                println!("assigning {:#?} to cell {:#?}", mutations, stem_cell);
            }
            mutate(stem_cell, mutations);
        } else if self.verbosity > 2 {
            println!("no mutations to assign to cell {:#?}", stem_cell);
        }
    }

    fn keep_const_population_upon_symmetric_division(&mut self, rng: &mut impl Rng) {
        //! If an symmetric division is performed, need to remove a random cell
        //! from the population.
        //! The cell having just proliferated can be removed as well, any cell
        //! can be removed.
        //!
        //! We proceed as following:
        //!     1. check if there is only one clone in the population, then
        //!     remove a cell from this clone
        //!     2. else, compute the variant counts
        //!     3. and sample from any clone based on the weights defined by
        //!     the variant counts
        if self.verbosity > 2 {
            println!("keeping the cell population constant");
        }
        // remove a cell from a random subclone based on the frequencies of
        // the clones at the current state
        let id2remove = if let Some(id) = self.subclones.the_only_one_subclone_present() {
            id
        } else {
            let variants = Variants::variant_counts(&self.subclones);
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
        time: f32,
    ) -> anyhow::Result<PathBuf> {
        let path2dir = self.path2dir.join(format!("{}cells", cells));
        let path2file = match tosave {
            Stats2Save::VariantFraction => path2dir.join("variant_fraction"),
            Stats2Save::Burden => path2dir.join("burden"),
            Stats2Save::Sfs => path2dir.join("sfs"),
        };
        let mut timepoint = format!("{:.1}", time).replace('.', "dot");
        timepoint.push_str("years");
        let path2file = path2file.join(timepoint);
        fs::create_dir_all(&path2file).with_context(|| "Cannot create dir")?;
        if self.verbosity > 1 {
            println!("creating dirs {:#?}", path2file);
        }
        Ok(path2file.join(self.filename.clone()))
    }

    pub fn save(
        &self,
        time: f32,
        nb_cells: usize,
        save_sfs_only: bool,
        rng: &mut impl Rng,
    ) -> anyhow::Result<()> {
        let cells_with_idx = if nb_cells == self.subclones.compute_tot_cells() as usize {
            self.subclones.get_cells_with_clones_idx()
        } else {
            self.subclones
                .get_cells_subsampled_with_clones_idx(nb_cells, rng)
        };
        assert_eq!(cells_with_idx.len(), nb_cells);
        let cells: Vec<&StemCell> = cells_with_idx.iter().map(|ele| ele.0).collect();

        if self.verbosity > 1 {
            println!("saving {} cells", cells.len());
        }

        Sfs::from_cells(&cells, self.verbosity)
            .unwrap_or_else(|_| panic!("cannot create SFS for timepoint at time {}", time))
            .save(&self.make_path(Stats2Save::Sfs, nb_cells, time)?)?;

        if !save_sfs_only {
            MutationalBurden::from_cells(&cells, self.verbosity)
                .unwrap_or_else(|_| {
                    panic!("cannot create burden for the timepoint at time {}", time)
                })
                .save(&self.make_path(Stats2Save::Burden, nb_cells, time)?)?;
            save_variant_fraction(
                &SubClones::from(
                    cells_with_idx
                        .into_iter()
                        .map(|(cell, id)| (cell.to_owned(), id))
                        .collect::<Vec<(StemCell, usize)>>(),
                ),
                &self.make_path(Stats2Save::VariantFraction, nb_cells, time)?,
                self.verbosity,
            )?;
        }

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
        //! Here we also save the state of the simulation at different
        //! timepoints, see `self.snapshot`.
        //!
        //! For every simulated timestep, the algorithm samples an exponential
        //! waiting time and the non-empty clone that will proliferate.
        //! From the proliferating clone selected, we sample a random cell to
        //! which we assign neutral and fit mutations.
        //!
        //! The proliferation step is implemented as following:
        //!
        //! 1. select the cell `c` that will proliferate next from the clone
        //! with id `reaction` determined by the Gillespie algorithm
        //!
        //! 2. draw and assign mb neutral background mutations to `c` from
        //! Poisson(DeltaT * mub), where mub is the backgound mutation rate
        //!
        //! 3. draw and assign md neutral division mutations to `c` from
        //! Poisson(mud), where mud is the division mutation rate
        //!
        //! 4. draw from a Bernoulli trial a fit mutation and in case assign
        //! `c` to a new clone. This clone must be empty and different from the
        //! old clone which `c` belonged to.
        //!
        //! 5. (optional) if that was a symmetric division, clone the
        //! proliferating cell `c` into `c1` and repeat step 3, 4 with `c1`.

        // 1. select a cell `c` at random
        // the SSA samples the clone that will proliferate next, that is the
        // clone with id `reaction.event`. Pick random proliferating cells from
        // this clone. **Note that this removes the cells from the clone with
        // id `reaction.event`.**
        self.time += reaction.time;
        let mut stem_cell =
            proliferating_cell(&mut self.subclones, reaction.event, self.verbosity, rng);
        self.counter_divisions += 1;

        // 2. draw background mutations and assign them to `c`
        assign_background_mutations(
            &mut stem_cell,
            self.time,
            &self.distributions.neutral_poisson,
            rng,
            self.verbosity,
        );
        if self.verbosity > 1 {
            println!("cell {:#?} is dividing", stem_cell);
            println!("at time {}", self.time);
        }

        // 3., 4. and 5.
        if self.distributions.can_only_be_symmetric()
            || !self.distributions.bern_asymmetric.sample(rng)
        {
            for mut cell in [stem_cell.clone(), stem_cell] {
                let division = self
                    .distributions
                    .neutral_poisson
                    .new_muts_upon_division(rng);
                if self.verbosity > 2 {
                    println!(
                        "assigning {} mutations upon cell division to cell {:#?}",
                        division.as_ref().unwrap_or(&vec![]).len(),
                        cell
                    );
                }
                self.assign_mutations(&mut cell, division);
                assign(
                    &mut self.subclones,
                    reaction.event,
                    cell,
                    &self.distributions,
                    rng,
                    self.verbosity,
                );
            }

            // remove a cell from the population
            self.keep_const_population_upon_symmetric_division(rng);
        } else {
            let division = self
                .distributions
                .neutral_poisson
                .new_muts_upon_division(rng);
            if self.verbosity > 2 {
                println!(
                    "assigning {} mutations upon cell division to cell {:#?}",
                    division.as_ref().unwrap_or(&vec![]).len(),
                    stem_cell
                )
            }
            self.assign_mutations(&mut stem_cell, division);
            assign(
                &mut self.subclones,
                reaction.event,
                stem_cell,
                &self.distributions,
                rng,
                self.verbosity,
            );
        }

        if self.verbosity > 2 {
            println!("{} cells", self.subclones.compute_tot_cells());
        }

        if !self.snapshots.is_empty() {
            // take snapshot
            if self.snapshots.iter().any(|s| self.time >= s.time) {
                let snapshot = self.snapshots.pop_front().unwrap();
                // this is important: we update all background mutations at this
                // time such that all cells are all on the same page.
                // Since background mutations are implemented at each division, and
                // cells do not proliferate at the same rate, we need to correct and
                // update the background mutations at the timepoint corresponding
                // to sampling step, i.e. before saving.
                for stem_cell in self.subclones.get_mut_cells() {
                    assign_background_mutations(
                        stem_cell,
                        self.time,
                        &self.distributions.neutral_poisson,
                        rng,
                        self.verbosity,
                    );
                }
                if self.verbosity > 0 {
                    println!(
                        "saving variant fraction for {:#?} at time {} for timepoint {}th with {} cells",
                        snapshot.time,
                        self.time,
                        self.snapshots.len(),
                        snapshot.cells2sample
                    );
                }
                self.save(self.time, snapshot.cells2sample, self.save_sfs_only, rng)
                    .expect("cannot save snapshot");
            }
        }
    }

    fn update_state(&self, state: &mut CurrentState<MAX_SUBCLONES>) {
        state.population = Variants::variant_counts(&self.subclones);
    }
}
