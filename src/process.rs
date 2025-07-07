use crate::genotype::{MutationalBurden, Sfs};
use crate::proliferation::{NeutralMutations, Proliferation};
use crate::stemcell::{assign_background_mutations, StemCell};
use crate::subclone::{save_variant_fraction, CloneId, Distributions, SubClones, Variants};
use crate::{MAX_SUBCLONES, TIME_AT_BIRTH};
use anyhow::Context;
use rand::Rng;
use rand_distr::{Distribution, WeightedIndex};
use sosa::{AdvanceStep, CurrentState, NextReaction};
use std::collections::VecDeque;
use std::fs;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy)]
pub enum SavingCells {
    WholePopulation,
    Subsampling {
        sample_size: usize,
        population_as_well: bool,
    },
}

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
    pub save_population: bool,
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
    pub proliferation: Proliferation,
    pub time: f32,
}

impl Exponential {
    pub fn new(
        initial_subclones: SubClones,
        distributions: Distributions,
        proliferation: Proliferation,
        verbosity: u8,
    ) -> Exponential {
        let hsc = Exponential {
            subclones: initial_subclones,
            distributions,
            counter_divisions: 0,
            proliferation,
            verbosity,
            time: 0f32,
        };
        if verbosity > 1 {
            println!("process created: {hsc:#?}");
        }
        hsc
    }

    pub fn switch_to_moran(
        mut self,
        process_options: ProcessOptions,
        distributions: Distributions,
        filename: PathBuf,
        save_sfs_only: bool,
        save_population: bool,
        rng: &mut impl Rng,
    ) -> Moran {
        //! End the exponential growing phase and switch to a fixed-size
        //! population phase.
        //!
        //! There is a delay between birth and the end of the exponential phase,
        //! since stem cells stop dividing exponentially before birth.
        //! We add background mutations in this interval of time between the
        //! end of the exponentially growing phase and the Moran process.
        if self.verbosity > 0 {
            println!("switching to Moran at time {}", self.time);
        }
        if let NeutralMutations::UponDivisionAndBackground = self.proliferation.neutral_mutation {
            // this is important: we update all background mutations at this
            // time such that all cells are all on the same page.
            // Since background mutations are implemented at each division, and
            // cells do not proliferate at the same rate, we need to correct and
            // update the background mutations at the timepoint corresponding
            // to sampling step, i.e. before saving.
            //
            // Moreover there is delay between the end of the exp. growing phase
            // and birth, hence we use `TIME_AT_BIRTH`.
            if self.verbosity > 0 {
                println!("updating the neutral background mutations for all cells");
            }
            for stem_cell in self.subclones.get_mut_cells() {
                if stem_cell.get_last_division_time() < &TIME_AT_BIRTH {
                    assign_background_mutations(
                        stem_cell,
                        TIME_AT_BIRTH,
                        &self.distributions.neutral_poisson,
                        rng,
                        self.verbosity,
                    );
                }
                // this is required as we are restarting the time
                stem_cell.set_last_division_time(0f32).unwrap();
            }
        }
        let mut moran = Moran {
            subclones: self.subclones,
            counter_divisions: self.counter_divisions,
            // restart the time
            time: 0.,
            snapshots: process_options.snapshots,
            path2dir: process_options.path,
            verbosity: self.verbosity,
            filename,
            distributions,
            save_sfs_only,
            save_population,
            proliferation: self.proliferation,
        };
        moran
            .save(
                moran.time,
                &SavingCells::WholePopulation,
                moran.save_sfs_only,
                rng,
            )
            .expect("cannot save simulation at the end of the exponential phase");
        moran
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
        //!    with id `reaction` determined by the Gillespie algorithm
        //!
        //! 2. for the proliferating cell and its daughter cell:
        //!     * mutate genome by storing a neutral number of mutations
        //!
        //!     * assign to new subclone with a probability determined by the
        //!       rate of mutations conferring a proliferative advantage
        // The Gillespie sampler samples the clone that will proliferate next,
        // that is the clone with id `reaction.event`.
        // Pick random proliferating cells from this clone. **Note that this
        // removes the cells from the clone with id `reaction.event`.**
        self.counter_divisions += 1;
        self.time += reaction.time;
        self.proliferation.proliferate(
            &mut self.subclones,
            self.time,
            reaction.event,
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
    pub time: f32,
    pub path2dir: PathBuf,
    pub verbosity: u8,
    pub distributions: Distributions,
    pub snapshots: VecDeque<Snapshot>,
    pub filename: PathBuf,
    pub save_sfs_only: bool,
    pub save_population: bool,
    pub proliferation: Proliferation,
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
                save_population: true,
            },
            Distributions::default(),
            Proliferation::default(),
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
        proliferation: Proliferation,
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
            save_population: saving_options.save_population,
            proliferation,
        };
        if verbosity > 1 {
            println!("process created: {hsc:#?}");
        }
        hsc
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
        let variants = Variants::variant_counts(&self.subclones);
        let id2remove = WeightedIndex::new(variants).unwrap().sample(rng);
        let _ = self
            .subclones
            .get_mut_clone_unchecked(id2remove)
            .random_cell(rng)
            .with_context(|| "found empty subclone")
            .unwrap();
        if self.verbosity > 2 {
            println!("removing one cell from clone {id2remove}");
        }
    }

    pub fn make_path(
        &self,
        tosave: Stats2Save,
        cells: usize,
        time: f32,
    ) -> anyhow::Result<PathBuf> {
        let path2dir = self.path2dir.join(format!("{cells}cells"));
        let path2file = match tosave {
            Stats2Save::VariantFraction => path2dir.join("variant_fraction"),
            Stats2Save::Burden => path2dir.join("burden"),
            Stats2Save::Sfs => path2dir.join("sfs"),
        };
        let mut timepoint = format!("{time:.1}").replace('.', "dot");
        timepoint.push_str("years");
        let path2file = path2file.join(timepoint);
        fs::create_dir_all(&path2file).with_context(|| "Cannot create dir")?;
        if self.verbosity > 1 {
            println!("creating dirs {path2file:#?}");
        }
        Ok(path2file.join(self.filename.clone()))
    }

    fn save_it(
        &self,
        time: f32,
        cells_with_idx: Vec<(&StemCell, usize)>,
        save_sfs_only: bool,
    ) -> anyhow::Result<()> {
        if self.verbosity > 0 {
            println!("saving data at time {time}");
        }
        let cells: Vec<&StemCell> = cells_with_idx.iter().map(|ele| ele.0).collect();
        let nb_cells = cells.len();

        if self.verbosity > 0 {
            println!("saving {nb_cells} cells");
        }

        Sfs::from_cells(&cells, self.verbosity)
            .unwrap_or_else(|_| panic!("cannot create SFS for timepoint at time {time}"))
            .save(
                &self.make_path(Stats2Save::Sfs, nb_cells, time)?,
                self.verbosity,
            )?;

        if !save_sfs_only {
            MutationalBurden::from_cells(&cells, self.verbosity)
                .unwrap_or_else(|_| panic!("cannot create burden for the timepoint at time {time}"))
                .save(
                    &self.make_path(Stats2Save::Burden, nb_cells, time)?,
                    self.verbosity,
                )?;
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
                "saved measurements after {} mutational events",
                self.counter_divisions
            );
        }
        Ok(())
    }

    pub fn save(
        &mut self,
        time: f32,
        saving_cells: &SavingCells,
        save_sfs_only: bool,
        rng: &mut impl Rng,
    ) -> anyhow::Result<()> {
        if self.verbosity > 0 {
            println!("saving process at time {time}");
        }
        if let NeutralMutations::UponDivisionAndBackground = self.proliferation.neutral_mutation {
            // this is important: we update all background mutations at this
            // time such that all cells are all on the same page.
            // Since background mutations are implemented at each division, and
            // cells do not proliferate at the same rate, we need to correct and
            // update the background mutations at the timepoint corresponding
            // to sampling step, i.e. before saving.
            if self.verbosity > 0 {
                println!("updating the neutral background mutations for all cells");
            }
            for stem_cell in self.subclones.get_mut_cells() {
                assign_background_mutations(
                    stem_cell,
                    self.time,
                    &self.distributions.neutral_poisson,
                    rng,
                    self.verbosity,
                );
            }
        }
        let population = self.subclones.get_cells_with_clones_idx();
        if self.verbosity > 0 {
            println!("saving {saving_cells:#?}");
        }
        match saving_cells {
            SavingCells::WholePopulation => self
                .save_it(time, population, save_sfs_only)
                .with_context(|| "cannot save the full population")
                .unwrap(),
            SavingCells::Subsampling {
                sample_size,
                population_as_well: populaiton_as_well,
            } => {
                if *populaiton_as_well {
                    self.save_it(time, population, save_sfs_only)
                        .with_context(|| "cannot save the full population")
                        .unwrap();
                }
                let cells_with_idx = self
                    .subclones
                    .get_cells_subsampled_with_clones_idx(*sample_size, rng);
                assert_eq!(cells_with_idx.len(), *sample_size);
                self.save_it(time, cells_with_idx, save_sfs_only)
                    .with_context(|| "cannot save the subsample")
                    .unwrap();
            }
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
        // take snapshot
        while !self.snapshots.is_empty() && self.snapshots.iter().any(|s| self.time >= s.time) {
            let snapshot = self.snapshots.pop_front().unwrap();
            if self.verbosity > 0 {
                println!(
                    "saving state for timepoint at time {:#?} at simulation's time {} with {} cells",
                    snapshot.time, self.time, snapshot.cells2sample
                );
            }
            let saving_cells =
                if snapshot.cells2sample == self.subclones.compute_tot_cells() as usize {
                    SavingCells::WholePopulation
                } else {
                    SavingCells::Subsampling {
                        sample_size: snapshot.cells2sample,
                        population_as_well: self.save_population,
                    }
                };
            self.save(self.time, &saving_cells, self.save_sfs_only, rng)
                .expect("cannot save snapshot");
        }

        // 1. select a cell `c` at random
        // the SSA samples the clone that will proliferate next, that is the
        // clone with id `reaction.event`. Pick random proliferating cells from
        // this clone. **Note that this removes the cells from the clone with
        // id `reaction.event`.**
        self.time += reaction.time;
        self.counter_divisions += 1;
        self.proliferation.proliferate(
            &mut self.subclones,
            self.time,
            reaction.event,
            &self.distributions,
            rng,
            self.verbosity,
        );

        // remove a cell from the population
        self.keep_const_population_upon_symmetric_division(rng);

        if self.verbosity > 2 {
            println!("{} cells", self.subclones.compute_tot_cells());
        }
    }

    fn update_state(&self, state: &mut CurrentState<MAX_SUBCLONES>) {
        state.population = Variants::variant_counts(&self.subclones);
    }
}

#[cfg(test)]
mod tests {
    use std::num::{NonZeroU64, NonZeroU8};

    use quickcheck_macros::quickcheck;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    use crate::{genotype, Probs};

    use super::*;

    struct MoranFitVariant {
        tot_cells: u64,
        process: Moran,
    }

    struct MoranNoFitVariant {
        tot_cells: u64,
        process: Moran,
    }

    fn create_moran_fit_variants(cells: NonZeroU64) -> MoranFitVariant {
        MoranFitVariant {
            tot_cells: cells.get(),
            process: create_moran(true, cells),
        }
    }

    fn create_moran_no_fit_variants(cells: NonZeroU64) -> MoranNoFitVariant {
        MoranNoFitVariant {
            tot_cells: cells.get(),
            process: create_moran(false, cells),
        }
    }

    fn create_moran(fit_variants: bool, cells: NonZeroU64) -> Moran {
        let mu = if fit_variants { 0.999 } else { 0. };
        Moran::new(
            ProcessOptions {
                path: PathBuf::default(),
                snapshots: VecDeque::default(),
            },
            SubClones::new(vec![StemCell::new(); cells.get() as usize], 3, 0),
            0f32,
            SavingOptions {
                filename: PathBuf::default(),
                save_sfs_only: true,
                save_population: true,
            },
            Distributions::new(Probs::new(1., 1., mu, 1., cells.get(), 0), 0),
            Proliferation::default(),
            0,
        )
    }

    fn assert_pop_const(moran: &Moran, cells: u64) {
        assert_eq!(moran.subclones.compute_tot_cells(), cells);
    }

    fn assert_not_all_cells_in_wild_type(moran: &Moran, cells: u64) {
        assert_ne!(moran.subclones.get_clone_unchecked(0).cell_count(), cells);
    }

    fn assert_all_cells_in_wild_type(moran: &Moran, cells: u64) {
        assert_eq!(moran.subclones.get_clone_unchecked(0).cell_count(), cells);
        assert_eq!(moran.subclones.the_only_one_subclone_present().unwrap(), 0);
    }

    #[quickcheck]
    fn advance_moran_no_fit_variant_test(cells: NonZeroU8, seed: u64) {
        let mut moran = create_moran_no_fit_variants(NonZeroU64::new(cells.get() as u64).unwrap());
        let reaction = NextReaction {
            time: 12.2,
            event: 0,
        };
        let burden_before = genotype::MutationalBurden::from_moran(&moran.process, 0).unwrap();
        assert_pop_const(&moran.process, moran.tot_cells);
        assert_all_cells_in_wild_type(&moran.process, moran.tot_cells);
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);

        moran.process.advance_step(reaction, rng);

        assert_pop_const(&moran.process, moran.tot_cells);
        assert_all_cells_in_wild_type(&moran.process, moran.tot_cells);
        let burden_after = genotype::MutationalBurden::from_moran(&moran.process, 0).unwrap();
        assert!(burden_before.0.keys().sum::<u16>() < burden_after.0.keys().sum::<u16>());
    }

    #[quickcheck]
    fn advance_moran_test(cells: NonZeroU8, seed: u64, new_clone: bool) {
        let mut moran = create_moran_fit_variants(NonZeroU64::new(cells.get() as u64).unwrap());
        let reaction = if new_clone {
            NextReaction {
                time: cells.get() as f32 * 0.9999,
                event: 0,
            }
        } else {
            NextReaction {
                time: 0.0002,
                event: 0,
            }
        };
        let burden_before = genotype::MutationalBurden::from_moran(&moran.process, 0).unwrap();
        assert_pop_const(&moran.process, moran.tot_cells);
        assert_all_cells_in_wild_type(&moran.process, moran.tot_cells);
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);

        moran.process.advance_step(reaction, rng);

        assert_pop_const(&moran.process, moran.tot_cells);
        if new_clone {
            assert_not_all_cells_in_wild_type(&moran.process, moran.tot_cells);
        }
        let burden_after = genotype::MutationalBurden::from_moran(&moran.process, 0).unwrap();
        assert!(burden_before.0.keys().sum::<u16>() <= burden_after.0.keys().sum::<u16>());
    }

    struct ExpFitVariant {
        tot_cells: u64,
        process: Exponential,
    }

    struct ExpNoFitVariant {
        tot_cells: u64,
        process: Exponential,
    }

    fn create_exp_fit_variants(cells: NonZeroU64) -> ExpFitVariant {
        ExpFitVariant {
            tot_cells: cells.get(),
            process: create_exp(true, cells),
        }
    }

    fn create_exp_no_fit_variants(cells: NonZeroU64) -> ExpNoFitVariant {
        ExpNoFitVariant {
            tot_cells: cells.get(),
            process: create_exp(false, cells),
        }
    }

    fn create_exp(fit_variants: bool, cells: NonZeroU64) -> Exponential {
        let mu = if fit_variants { 0.999 } else { 0. };
        Exponential::new(
            SubClones::new(vec![StemCell::new(); cells.get() as usize], 3, 0),
            Distributions::new(Probs::new(10., 1., mu, 1., cells.get(), 0), 0),
            Proliferation::default(),
            0,
        )
    }

    fn assert_not_all_cells_in_wild_type_exp(exp: &Exponential, cells: u64) {
        assert_ne!(exp.subclones.get_clone_unchecked(0).cell_count(), cells);
    }

    fn assert_all_cells_in_wild_type_exp(exp: &Exponential, cells: u64) {
        assert_eq!(exp.subclones.get_clone_unchecked(0).cell_count(), cells);
        assert_eq!(exp.subclones.the_only_one_subclone_present().unwrap(), 0);
    }

    #[quickcheck]
    fn advance_exp_no_fit_variant_test(cells: NonZeroU8, seed: u64) {
        let mut exp = create_exp_no_fit_variants(NonZeroU64::new(cells.get() as u64).unwrap());
        let reaction = NextReaction {
            time: 12.2,
            event: 0,
        };
        let burden_before = genotype::MutationalBurden::from_exp(&exp.process, 0).unwrap();
        let old_cells = exp.process.subclones.compute_tot_cells();
        assert_all_cells_in_wild_type_exp(&exp.process, exp.tot_cells);
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);

        exp.process.advance_step(reaction, rng);

        assert_eq!(exp.process.subclones.compute_tot_cells(), old_cells + 1);
        let burden_after = genotype::MutationalBurden::from_exp(&exp.process, 0).unwrap();
        assert!(burden_before.0.keys().sum::<u16>() <= burden_after.0.keys().sum::<u16>());
    }

    #[quickcheck]
    fn advance_exp_test(cells: NonZeroU8, seed: u64, new_clone: bool) {
        let mut exp = create_exp_fit_variants(NonZeroU64::new(cells.get() as u64).unwrap());
        let reaction = if new_clone {
            NextReaction {
                time: cells.get() as f32 * 0.9999,
                event: 0,
            }
        } else {
            NextReaction {
                time: 0.0002,
                event: 0,
            }
        };
        let burden_before = genotype::MutationalBurden::from_exp(&exp.process, 0).unwrap();
        let old_cells = exp.process.subclones.compute_tot_cells();
        assert_all_cells_in_wild_type_exp(&exp.process, exp.tot_cells);
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);

        exp.process.advance_step(reaction, rng);

        assert_eq!(exp.process.subclones.compute_tot_cells(), old_cells + 1);
        if new_clone {
            assert_not_all_cells_in_wild_type_exp(&exp.process, exp.tot_cells);
        }
        let burden_after = genotype::MutationalBurden::from_exp(&exp.process, 0).unwrap();
        assert!(burden_before.0.keys().sum::<u16>() <= burden_after.0.keys().sum::<u16>());
    }
}
