use crate::proliferation::Proliferation;
use crate::snapshots::{SavingCells, SavingOptions, Snapshot, StatsConfig, save_it};
use crate::subclone::{
    CloneId, Distributions, Fitness, FitnessConfig, RateSampler, SubClones, Variants,
};
use crate::{MAX_SUBCLONES, TIME_AT_BIRTH};
use anyhow::Context;
use log::{debug, info, trace};
use rand::{Rng, SeedableRng};
use rand_distr::Distribution;
use rand_distr::weighted::WeightedIndex;
use sosa::{
    AdvanceStep, CurrentState, IterTime, NextReaction, ReactionRates, SimState, StopReason, exprand,
};
use std::collections::VecDeque;
use std::num::NonZeroUsize;
use std::path::PathBuf;

/// Sentinel rates passed to [`sosa::simulate`] for the [`Moran`] phase.
///
/// [`Moran`] overrides [`AdvanceStep::next_reaction`] to read its own
/// `rates` field, so the `&ReactionRates` argument given to `sosa::simulate`
/// is ignored. We still need *some* value to pass — that's this static.
pub static DUMMY_RATES: ReactionRates<MAX_SUBCLONES> = ReactionRates([0.0; MAX_SUBCLONES]);

/// Output/configuration bundle used by [`Moran::new`] and [`switch_to_moran`].
///
/// Groups the three knobs that govern what the Moran phase writes to disk
/// (path, snapshots, filename, save_population, per-snapshot stats).
#[derive(Debug, Clone)]
pub struct MoranConfig {
    pub process_options: ProcessOptions,
    pub saving_options: SavingOptions,
    pub stats: StatsConfig,
}

#[derive(Debug, Clone, Default)]
pub struct ProcessOptions {
    pub path: PathBuf,
    pub snapshots: VecDeque<Snapshot>,
}

#[derive(Debug, Clone)]
pub struct CellDivisionProbabilities {
    /// Rate of neutral mutations per cell-division
    pub lambda_poisson: f32,
    /// Probability of getting a fit mutant upon cell division
    pub p: f64,
}

#[derive(Clone, Debug)]
pub enum ExponentialPhase {
    /// The embryonic phase up to birth
    Development,
    /// A regrowth exponential phase upon subsampling or treatment.
    Regrowth,
}

pub fn switch_to_moran(
    mut exponential: Exponential,
    config: MoranConfig,
    distributions: Distributions,
    rate_sampler: RateSampler,
    rng: &mut impl Rng,
) -> Moran {
    //! End the exponential growing phase and switch to a fixed-size
    //! population phase.
    //!
    //! There is a delay between birth and the end of the exponential phase,
    //! since stem cells stop dividing exponentially before birth.
    //! We add background mutations in this interval of time between the
    //! end of the exponentially growing phase and the Moran process.
    debug!("switching to Moran at time {}", exponential.time);
    let (current_time, time_to_restart) = match exponential.phase {
        ExponentialPhase::Development => (TIME_AT_BIRTH, 0f32),
        ExponentialPhase::Regrowth => (exponential.time, exponential.time),
    };
    // Realise the lazily-deferred background mutations up to `current_time`.
    // For Development that is `TIME_AT_BIRTH` (there is a delay between the
    // end of the exp. growing phase and birth); for Regrowth it is the
    // current process time. The clock-reset below is a separate
    // phase-transition concern that stays in this function.
    exponential.proliferation.realise_background_mutations(
        &mut exponential.subclones,
        current_time,
        &exponential.distributions,
        rng,
    );
    for stem_cell in exponential.subclones.get_mut_cells() {
        stem_cell.set_last_division_time(time_to_restart).unwrap();
    }
    let rates = rate_sampler.rates_for(&exponential.subclones, rng);
    Moran {
        subclones: exponential.subclones,
        counter_divisions: exponential.counter_divisions,
        // restart the time within the process
        time: time_to_restart,
        snapshots: config.process_options.snapshots,
        path2dir: config.process_options.path,
        filename: config.saving_options.filename,
        distributions,
        stats: config.stats,
        save_population: config.saving_options.save_population,
        proliferation: exponential.proliferation,
        rates,
        rate_sampler,
    }
}

/// Exponential growing process
#[derive(Clone, Debug)]
pub struct Exponential {
    /// A collection of clones having a proliferative advantage.
    pub subclones: SubClones,
    /// The counter for the number of proliferative events.
    pub counter_divisions: usize,
    pub distributions: Distributions,
    pub proliferation: Proliferation,
    pub time: f32,
    pub phase: ExponentialPhase,
}

impl Exponential {
    pub fn new(
        initial_subclones: SubClones,
        distributions: Distributions,
        proliferation: Proliferation,
        time: f32,
        phase: ExponentialPhase,
    ) -> Exponential {
        let hsc = Exponential {
            subclones: initial_subclones,
            distributions,
            counter_divisions: 0,
            proliferation,
            phase,
            time,
        };
        debug!("process created: {hsc:#?}");
        hsc
    }
}

impl From<Moran> for Exponential {
    fn from(value: Moran) -> Self {
        //! Construct an exponential-growing process in [`ExponentialPhase::Regrowth`] mode.
        //!
        //! This is useful when we undersample the Moran process to a smaller
        //! population size (for example upon treatment, the number of cells is
        //! reduced) and we want to grow back the population to its size before
        //! the undersampling.
        Exponential::new(
            value.subclones,
            value.distributions,
            value.proliferation,
            value.time,
            ExponentialPhase::Regrowth,
        )
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
        let _ = self.proliferation.proliferate(
            &mut self.subclones,
            self.time,
            reaction.event,
            &self.distributions,
            rng,
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
    pub distributions: Distributions,
    pub snapshots: VecDeque<Snapshot>,
    pub filename: PathBuf,
    pub save_population: bool,
    pub proliferation: Proliferation,
    pub stats: StatsConfig,
    /// Per-clone Gillespie rates owned by this Moran. Read by the overridden
    /// [`AdvanceStep::next_reaction`]; the `&ReactionRates` argument that
    /// [`sosa::simulate`] supplies is ignored — see [`DUMMY_RATES`].
    pub rates: ReactionRates<MAX_SUBCLONES>,
    /// Owns the rate-sampling strategy. Used to draw [`Self::rates`] at
    /// construction and (in a later commit) to resample on births.
    pub rate_sampler: RateSampler,
}

impl Default for Moran {
    fn default() -> Self {
        let config = MoranConfig {
            process_options: ProcessOptions {
                path: PathBuf::from("./output"),
                snapshots: VecDeque::default(),
            },
            saving_options: SavingOptions {
                filename: PathBuf::default(),
                save_population: true,
            },
            stats: StatsConfig::default(),
        };
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(0);
        let rate_sampler = RateSampler::from_config(0.0, &FitnessConfig::Static(Fitness::Neutral));

        Moran::new(
            config,
            SubClones::default(),
            0.,
            Distributions::default(),
            Proliferation::default(),
            rate_sampler,
            &mut rng,
        )
    }
}

impl Moran {
    /// A Moran process with wild-type subclone with neutral fitness, to which
    /// all cells will be assigned.
    pub fn new(
        config: MoranConfig,
        initial_subclones: SubClones,
        time: f32,
        distributions: Distributions,
        proliferation: Proliferation,
        rate_sampler: RateSampler,
        rng: &mut impl Rng,
    ) -> Moran {
        let rates = rate_sampler.rates_for(&initial_subclones, rng);
        let hsc = Moran {
            subclones: initial_subclones,
            distributions,
            counter_divisions: 0,
            path2dir: config.process_options.path,
            time,
            snapshots: config.process_options.snapshots,
            filename: config.saving_options.filename,
            stats: config.stats,
            save_population: config.saving_options.save_population,
            proliferation,
            rates,
            rate_sampler,
        };
        debug!("process created: {hsc:#?}");
        hsc
    }

    fn keep_const_population_upon_symmetric_division(&mut self, rng: &mut impl Rng) {
        //! If an symmetric division is performed, need to remove a random cell
        //! from the population.
        //! The cell having just proliferated can be removed as well, any cell
        //! can be removed.
        //!
        //! Compute compute the variant counts and sample from any clone based
        //! on the weights defined by the variant counts.
        //! This is a complicated way to acces one random cell in the
        //! population.
        trace!("keeping the cell population constant");
        // remove a cell from a random subclone based on the frequencies of
        // the clones at the current state
        let variants = Variants::variant_counts(&self.subclones);
        let id2remove = WeightedIndex::new(variants).unwrap().sample(rng) as CloneId;
        let _ = self
            .subclones
            .get_mut_clone_unchecked(id2remove)
            .random_cell(rng)
            .with_context(|| "found empty subclone")
            .unwrap();
        trace!("removing one cell from clone {id2remove}");
    }

    pub fn save(
        &mut self,
        time: f32,
        saving_cells: &SavingCells,
        rng: &mut impl Rng,
    ) -> anyhow::Result<()> {
        info!("saving process at time {time}");
        // realise the lazily-deferred background mutations so that every cell
        // is up to date before we read the mutation state for output.
        self.proliferation.realise_background_mutations(
            &mut self.subclones,
            self.time,
            &self.distributions,
            rng,
        );
        let population = self.subclones.get_cells_with_clones_idx();
        trace!("saving {saving_cells:#?}");
        match saving_cells {
            SavingCells::WholePopulation => save_it(
                &self.path2dir,
                &self.filename,
                time,
                population,
                &self.subclones,
                &self.rates,
                &self.stats,
            )
            .with_context(|| "cannot save the full population")
            .unwrap(),
            SavingCells::Subsampling {
                sample_size,
                population_as_well: populaiton_as_well,
            } => {
                if *populaiton_as_well {
                    save_it(
                        &self.path2dir,
                        &self.filename,
                        time,
                        population,
                        &self.subclones,
                        &self.rates,
                        &self.stats,
                    )
                    .with_context(|| "cannot save the full population")
                    .unwrap();
                }
                let cells_with_idx = self
                    .subclones
                    .get_cells_subsampled_with_clones_idx(*sample_size, rng);
                assert_eq!(cells_with_idx.len(), *sample_size);
                save_it(
                    &self.path2dir,
                    &self.filename,
                    time,
                    cells_with_idx,
                    &self.subclones,
                    &self.rates,
                    &self.stats,
                )
                .with_context(|| "cannot save the subsample")
                .unwrap();
            }
        }
        Ok(())
    }

    pub fn into_subsampled(self, nb_cells: NonZeroUsize, rng: &mut impl Rng) -> Self {
        //! Take a subsample (without replacement) of the population of this
        //! process using [`choose_multiple`](https://docs.rs/rand/latest/rand/seq/trait.IndexedRandom.html#method.choose_multiple).
        let mut moran: Moran = self;
        moran.subclones = moran.subclones.into_subsampled(nb_cells, rng);
        debug!("Subsampled {} cells", moran.subclones.get_cells().len());

        moran
    }
}

impl AdvanceStep<MAX_SUBCLONES> for Moran {
    /// The id of a subclone [`CloneId`].
    type Reaction = CloneId;

    fn next_reaction(
        &self,
        state: &CurrentState<MAX_SUBCLONES>,
        _rates: &ReactionRates<MAX_SUBCLONES>,
        possible_reactions: &[Self::Reaction; MAX_SUBCLONES],
        max_iter_and_time: IterTime,
        iter_and_time: IterTime,
        rng: &mut impl Rng,
    ) -> (SimState, Option<NextReaction<Self::Reaction>>) {
        //! Override of the default trait impl: reads Gillespie rates from
        //! `self.rates` instead of the `_rates` argument, which lets the
        //! `--multihits` feature mutate rates between steps. Body mirrors
        //! `sosa-5.1.0/src/lib.rs:282-321`; the inlined Gillespie waiting-time
        //! computation matches `ReactionRates::compute_times_events` (which is
        //! crate-private in sosa).
        if state.population.iter().sum::<u64>() == 0u64 {
            return (SimState::Stop(StopReason::NoIndividualsLeft), None);
        };
        if iter_and_time.iter >= max_iter_and_time.iter {
            return (SimState::Stop(StopReason::MaxItersReached), None);
        };
        if iter_and_time.time >= max_iter_and_time.time {
            return (SimState::Stop(StopReason::MaxTimeReached), None);
        };

        let mut times = self.rates.0;
        for (i, t) in times.iter_mut().enumerate() {
            *t = exprand(*t * state.population[i] as f32, rng);
        }
        assert!(
            times.iter().any(|t| t.is_normal()),
            "all Gillespie waiting times are non-normal",
        );

        let mut selected_event = 0_usize;
        let mut smaller_waiting_time = times[0];
        for (idx, &waiting_time) in times.iter().enumerate() {
            if waiting_time <= smaller_waiting_time {
                smaller_waiting_time = waiting_time;
                selected_event = idx;
            }
        }
        (
            SimState::Continue,
            Some(NextReaction {
                time: smaller_waiting_time,
                event: possible_reactions[selected_event],
            }),
        )
    }

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
            debug!(
                "saving state for timepoint at simulation's time {} for timepoint at {:#?} with {} cells",
                self.time, snapshot.time, snapshot.cells2sample
            );
            let saving_cells =
                if snapshot.cells2sample == self.subclones.compute_tot_cells() as usize {
                    SavingCells::WholePopulation
                } else {
                    SavingCells::Subsampling {
                        sample_size: snapshot.cells2sample,
                        population_as_well: self.save_population,
                    }
                };
            self.save(self.time, &saving_cells, rng)
                .expect("cannot save snapshot");
        }

        // 1. select a cell `c` at random
        // the SSA samples the clone that will proliferate next, that is the
        // clone with id `reaction.event`. Pick random proliferating cells from
        // this clone. **Note that this removes the cells from the clone with
        // id `reaction.event`.**
        self.time += reaction.time;
        self.counter_divisions += 1;
        let birth = self.proliferation.proliferate(
            &mut self.subclones,
            self.time,
            reaction.event,
            &self.distributions,
            rng,
        );
        // `--multihits` hook: when the new clone's hit-count tier differs from
        // what `rate_sampler.rates_for` pre-sampled at construction
        // (`RateSamplerMode::Multihits` + 2nd-or-later hit), redraw its rate
        // from the matching tier and write it into `self.rates`. Static mode
        // and 1st hits skip the redraw — the pre-sampled value is already
        // correct.
        if let Some(event) = birth {
            if self.rate_sampler.mode.should_resample(event.hit_count) {
                let new_rate = self.rate_sampler.pick(event.hit_count, rng);
                debug!(
                    "resampling rate for new clone {} (parent={}, hit_count={:?}): {new_rate}",
                    event.new_clone_id, event.parent_clone_id, event.hit_count,
                );
                self.rates.0[event.new_clone_id as usize] = new_rate;
            }
        }

        // remove a cell from the population
        self.keep_const_population_upon_symmetric_division(rng);

        trace!("{} cells", self.subclones.compute_tot_cells());
    }

    fn update_state(&self, state: &mut CurrentState<MAX_SUBCLONES>) {
        state.population = Variants::variant_counts(&self.subclones);
    }
}

#[cfg(test)]
mod tests {
    use std::num::{NonZeroU8, NonZeroU64};

    use quickcheck_macros::quickcheck;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    use crate::{Probs, genotype, stemcell::StemCell, subclone::HitCount};

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
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let rate_sampler = RateSampler::from_config(1.0, &FitnessConfig::Static(Fitness::Neutral));
        Moran::new(
            MoranConfig {
                process_options: ProcessOptions {
                    path: PathBuf::default(),
                    snapshots: VecDeque::default(),
                },
                saving_options: SavingOptions {
                    filename: PathBuf::default(),
                    save_population: true,
                },
                stats: StatsConfig::default(),
            },
            SubClones::new(vec![StemCell::new(); cells.get() as usize], 3),
            0f32,
            Distributions::new(Probs::new(1., 1., mu, 1., cells.get())),
            Proliferation::default(),
            rate_sampler,
            &mut rng,
        )
    }

    fn create_moran_multihits_with_first_hit_seed(half: u64) -> Moran {
        // Build a Moran in Multihits mode with `half` cells in clone 0
        // (wild-type) and `half` cells in clone 1, where clone 1 is marked as
        // a First-hit clone. Both clones are non-empty so `next_clone` cannot
        // select them as the new clone target, and a birth from clone 1
        // produces a child with `HitCount::Second`, which triggers the
        // resample hook.
        let cells = 2 * half;
        let mut subclones = SubClones::new(vec![StemCell::new(); half as usize], cells as usize);
        {
            let clone1 = subclones.get_mut_clone_unchecked(1);
            for _ in 0..half {
                clone1.assign_cell(StemCell::new());
            }
            clone1.set_hit_count(HitCount::First);
            clone1.set_parent_id(0);
        }
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let rate_sampler = RateSampler::from_config(1.0, &FitnessConfig::Multihits);
        Moran::new(
            MoranConfig {
                process_options: ProcessOptions::default(),
                saving_options: SavingOptions {
                    filename: PathBuf::default(),
                    save_population: false,
                },
                stats: StatsConfig::default(),
            },
            subclones,
            0f32,
            Distributions::new(Probs::new(1., 1., 0.999, 1., cells)),
            Proliferation::default(),
            rate_sampler,
            &mut rng,
        )
    }

    fn assert_pop_const(moran: &Moran, cells: u64) {
        assert_eq!(moran.subclones.compute_tot_cells(), cells);
    }

    fn assert_not_all_cells_in_wild_type(moran: &Moran, cells: u64) {
        assert_ne!(moran.subclones.get_clone(0).unwrap().cell_count(), cells);
    }

    fn assert_all_cells_in_wild_type(moran: &Moran, cells: u64) {
        assert_eq!(moran.subclones.get_clone(0).unwrap().cell_count(), cells);
        assert_eq!(the_only_one_subclone_present(&moran.subclones).unwrap(), 0);
    }

    #[quickcheck]
    fn advance_moran_no_fit_variant_test(cells: NonZeroU8, seed: u64) {
        let mut moran = create_moran_no_fit_variants(NonZeroU64::new(cells.get() as u64).unwrap());
        let reaction = NextReaction {
            time: 12.2,
            event: 0,
        };
        let burden_before = genotype::MutationalBurden::from_moran(&moran.process).unwrap();
        assert_pop_const(&moran.process, moran.tot_cells);
        assert_all_cells_in_wild_type(&moran.process, moran.tot_cells);
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);

        moran.process.advance_step(reaction, rng);

        assert_pop_const(&moran.process, moran.tot_cells);
        assert_all_cells_in_wild_type(&moran.process, moran.tot_cells);
        let burden_after = genotype::MutationalBurden::from_moran(&moran.process).unwrap();
        assert!(burden_before.0.keys().sum::<u16>() < burden_after.0.keys().sum::<u16>());
    }

    #[quickcheck]
    fn simulate_moran_uses_self_rates_and_preserves_population(cells: NonZeroU8, seed: u64) {
        // End-to-end: run sosa::simulate against a Moran process, passing
        // DUMMY_RATES as the rates argument. The override on Moran::next_reaction
        // ignores that argument and reads `self.rates` instead; if it didn't,
        // every waiting time would be infinite (DUMMY_RATES is all zeros) and
        // the assert in next_reaction would fire. We also confirm the Moran
        // invariant — population stays constant. Uses the no-fit-variants
        // Moran so simulation time can grow freely without the p>1 panic in
        // next_clone.
        use sosa::{Options, simulate};
        let cells = NonZeroU64::new(cells.get() as u64).unwrap();
        let mut moran = create_moran_no_fit_variants(cells);
        let initial = moran.process.subclones.compute_tot_cells();
        let state = &mut CurrentState {
            population: Variants::variant_counts(&moran.process.subclones),
        };
        let possible_reactions = moran.process.subclones.array_of_gillespie_reactions();
        let options = Options {
            init_iter: 0,
            max_iter_time: IterTime {
                iter: 8,
                time: f32::INFINITY,
            },
            max_cells: u64::MAX,
        };
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);

        let _ = simulate(
            state,
            &DUMMY_RATES,
            &possible_reactions,
            &mut moran.process,
            &options,
            rng,
        );

        assert_eq!(initial, moran.process.subclones.compute_tot_cells());
    }

    #[quickcheck]
    fn advance_step_static_keeps_rates_unchanged(seed: u64) {
        // Static mode: `should_resample` is false for every HitCount, so
        // `self.rates` must stay byte-equal across many advance_step calls
        // — births may or may not fire, rates do not move either way.
        // Fixed-size population (8 cells) and small reaction time keep the
        // Bernoulli `p` inside `next_clone` below 1.
        let mut moran = create_moran_fit_variants(NonZeroU64::new(8).unwrap()).process;
        let initial = moran.rates.0;
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        for _ in 0..16 {
            if moran.subclones.get_clone(0).unwrap().cell_count() == 0 {
                break;
            }
            let reaction = NextReaction {
                time: 0.5,
                event: 0,
            };
            moran.advance_step(reaction, &mut rng);
        }
        assert_eq!(initial, moran.rates.0);
    }

    #[quickcheck]
    fn advance_step_multihits_first_hit_does_not_resample(seed: u64) {
        // Birth from a wild-type parent produces a child with HitCount::First.
        // `should_resample(First)` is false (the slot's pre-sampled rate from
        // tier 0 is already what a first hit expects), so rates must stay
        // unchanged whether or not a birth fires.
        let cells = 8_u64;
        let mut rng_init = ChaCha8Rng::seed_from_u64(0);
        let mut moran = Moran::new(
            MoranConfig {
                process_options: ProcessOptions::default(),
                saving_options: SavingOptions {
                    filename: PathBuf::default(),
                    save_population: false,
                },
                stats: StatsConfig::default(),
            },
            SubClones::new(vec![StemCell::new(); cells as usize], 3),
            0f32,
            Distributions::new(Probs::new(1., 1., 0.999, 1., cells)),
            Proliferation::default(),
            RateSampler::from_config(1.0, &FitnessConfig::Multihits),
            &mut rng_init,
        );
        let initial = moran.rates.0;
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let reaction = NextReaction {
            time: 0.5,
            event: 0,
        };
        moran.advance_step(reaction, &mut rng);
        assert_eq!(initial, moran.rates.0);
    }

    #[test]
    fn advance_step_multihits_second_hit_resamples_target_slot() {
        // Parent clone 1 has HitCount::First; a fit birth produces a child
        // with HitCount::Second, which `should_resample` flags as true. The
        // new clone's slot must be overwritten with a fresh draw from the
        // matching tier; every other slot must be byte-equal to the initial
        // table.
        //
        // The Bernoulli inside `next_clone` may roll false for some seeds, in
        // which case no birth happens. We sweep a few seeds and assert that
        // at least one triggers a resample with the expected invariants.
        let mut saw_resample = false;
        for seed in 0..32_u64 {
            let mut moran = create_moran_multihits_with_first_hit_seed(4);
            let initial = moran.rates.0;
            let mut rng = ChaCha8Rng::seed_from_u64(seed);
            // u = mu / cells = 0.999/8 ≈ 0.125; with reaction.time = 4 the
            // Bernoulli p = u * time ≈ 0.5 — well under 1 and high enough that
            // most seeds fire a birth.
            let reaction = NextReaction {
                time: 4.0,
                event: 1,
            };
            moran.advance_step(reaction, &mut rng);
            let changed: Vec<usize> = (0..MAX_SUBCLONES)
                .filter(|&i| initial[i] != moran.rates.0[i])
                .collect();
            if changed.is_empty() {
                continue;
            }
            assert_eq!(
                changed.len(),
                1,
                "expected at most one slot to be resampled, got {changed:?}",
            );
            let slot = changed[0];
            assert!(
                slot >= 2,
                "neither wild-type nor First-hit slot should change"
            );
            assert!(
                moran.rates.0[slot].is_finite() && moran.rates.0[slot] > 0.0,
                "resampled rate must be finite positive, got {}",
                moran.rates.0[slot],
            );
            assert_eq!(
                moran
                    .subclones
                    .get_clone(slot as CloneId)
                    .unwrap()
                    .hit_count(),
                HitCount::Second,
            );
            saw_resample = true;
            break;
        }
        assert!(saw_resample, "no seed triggered a resample-eligible birth");
    }

    #[quickcheck]
    fn moran_new_rates_match_gillespie_rates_for_static(cells: NonZeroU8, seed: u64) {
        // Moran::new(rate_sampler, rng) must initialise self.rates identically
        // to the old setup of `subclones.gillespie_rates(&fitness, b0, rng)`
        // for the Static case, given the same RNG sequence.
        let cells = NonZeroU64::new(cells.get() as u64).unwrap();
        let b0 = 1.5_f32;
        let fitness = Fitness::GammaSampled {
            shape: 2.0,
            scale: 0.5,
        };
        let make_subclones = || SubClones::new(vec![StemCell::new(); cells.get() as usize], 3);
        let subclones = make_subclones();

        let mut rng_old = ChaCha8Rng::seed_from_u64(seed);
        let expected = subclones.gillespie_rates(&fitness, b0, &mut rng_old);

        let mut rng_new = ChaCha8Rng::seed_from_u64(seed);
        let moran = Moran::new(
            MoranConfig {
                process_options: ProcessOptions::default(),
                saving_options: SavingOptions {
                    filename: PathBuf::default(),
                    save_population: false,
                },
                stats: StatsConfig::default(),
            },
            make_subclones(),
            0.,
            Distributions::default(),
            Proliferation::default(),
            RateSampler::from_config(b0, &FitnessConfig::Static(fitness)),
            &mut rng_new,
        );
        assert_eq!(moran.rates.0, expected.0);
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
        let burden_before = genotype::MutationalBurden::from_moran(&moran.process).unwrap();
        assert_pop_const(&moran.process, moran.tot_cells);
        assert_all_cells_in_wild_type(&moran.process, moran.tot_cells);
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);

        moran.process.advance_step(reaction, rng);

        assert_pop_const(&moran.process, moran.tot_cells);
        if new_clone {
            assert_not_all_cells_in_wild_type(&moran.process, moran.tot_cells);
        }
        let burden_after = genotype::MutationalBurden::from_moran(&moran.process).unwrap();
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

    fn create_exp_fit_variants(cells: NonZeroU64, is_regrowth: bool) -> ExpFitVariant {
        ExpFitVariant {
            tot_cells: cells.get(),
            process: create_exp(true, cells, is_regrowth),
        }
    }

    fn create_exp_no_fit_variants(cells: NonZeroU64, is_regrowth: bool) -> ExpNoFitVariant {
        ExpNoFitVariant {
            tot_cells: cells.get(),
            process: create_exp(false, cells, is_regrowth),
        }
    }

    fn create_exp(fit_variants: bool, cells: NonZeroU64, is_regrowht: bool) -> Exponential {
        let mu = if fit_variants { 0.999 } else { 0. };
        let phase = if is_regrowht {
            ExponentialPhase::Regrowth
        } else {
            ExponentialPhase::Development
        };
        Exponential::new(
            SubClones::new(vec![StemCell::new(); cells.get() as usize], 3),
            Distributions::new(Probs::new(10., 1., mu, 1., cells.get())),
            Proliferation::default(),
            0.,
            phase,
        )
    }

    fn assert_not_all_cells_in_wild_type_exp(exp: &Exponential, cells: u64) {
        assert_ne!(exp.subclones.get_clone(0).unwrap().cell_count(), cells);
    }

    fn assert_all_cells_in_wild_type_exp(exp: &Exponential, cells: u64) {
        assert_eq!(exp.subclones.get_clone(0).unwrap().cell_count(), cells);
        assert_eq!(the_only_one_subclone_present(&exp.subclones).unwrap(), 0);
    }

    #[quickcheck]
    fn advance_exp_no_fit_variant_test(cells: NonZeroU8, seed: u64) {
        let mut exp =
            create_exp_no_fit_variants(NonZeroU64::new(cells.get() as u64).unwrap(), false);
        let reaction = NextReaction {
            time: 12.2,
            event: 0,
        };
        let burden_before = genotype::MutationalBurden::from_exp(&exp.process).unwrap();
        let old_cells = exp.process.subclones.compute_tot_cells();
        assert_all_cells_in_wild_type_exp(&exp.process, exp.tot_cells);
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);

        exp.process.advance_step(reaction, rng);

        assert_eq!(exp.process.subclones.compute_tot_cells(), old_cells + 1);
        let burden_after = genotype::MutationalBurden::from_exp(&exp.process).unwrap();
        assert!(burden_before.0.keys().sum::<u16>() <= burden_after.0.keys().sum::<u16>());
    }

    #[quickcheck]
    fn advance_exp_test(cells: NonZeroU8, seed: u64, new_clone: bool) {
        let mut exp = create_exp_fit_variants(NonZeroU64::new(cells.get() as u64).unwrap(), false);
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
        let burden_before = genotype::MutationalBurden::from_exp(&exp.process).unwrap();
        let old_cells = exp.process.subclones.compute_tot_cells();
        assert_all_cells_in_wild_type_exp(&exp.process, exp.tot_cells);
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);

        exp.process.advance_step(reaction, rng);

        assert_eq!(exp.process.subclones.compute_tot_cells(), old_cells + 1);
        if new_clone {
            assert_not_all_cells_in_wild_type_exp(&exp.process, exp.tot_cells);
        }
        let burden_after = genotype::MutationalBurden::from_exp(&exp.process).unwrap();
        assert!(burden_before.0.keys().sum::<u16>() <= burden_after.0.keys().sum::<u16>());
    }

    #[quickcheck]
    fn switch_to_moran_development_test(cells: NonZeroU8, seed: u64) {
        let exp = create_exp_fit_variants(NonZeroU64::new(cells.get() as u64).unwrap(), false);
        let burden_before = genotype::MutationalBurden::from_exp(&exp.process).unwrap();
        let old_cells = exp.process.subclones.compute_tot_cells();
        assert_all_cells_in_wild_type_exp(&exp.process, exp.tot_cells);
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);
        let moran = switch_to_moran(
            exp.process,
            MoranConfig {
                process_options: ProcessOptions::default(),
                saving_options: SavingOptions {
                    filename: PathBuf::default(),
                    save_population: false,
                },
                stats: StatsConfig::default(),
            },
            Distributions::default(),
            RateSampler::from_config(1.0, &FitnessConfig::Static(Fitness::Neutral)),
            rng,
        );

        let burden_after = genotype::MutationalBurden::from_moran(&moran).unwrap();
        assert!(burden_before.0.keys().sum::<u16>() <= burden_after.0.keys().sum::<u16>());
        assert_eq!(old_cells, moran.subclones.compute_tot_cells());
    }

    #[quickcheck]
    fn switch_to_moran_regrowth_test(cells: NonZeroU8, seed: u64) {
        let exp = create_exp_fit_variants(NonZeroU64::new(cells.get() as u64).unwrap(), true);
        let burden_before = genotype::MutationalBurden::from_exp(&exp.process).unwrap();
        let old_cells = exp.process.subclones.compute_tot_cells();
        assert_all_cells_in_wild_type_exp(&exp.process, exp.tot_cells);
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);
        let moran = switch_to_moran(
            exp.process,
            MoranConfig {
                process_options: ProcessOptions::default(),
                saving_options: SavingOptions {
                    filename: PathBuf::default(),
                    save_population: false,
                },
                stats: StatsConfig::default(),
            },
            Distributions::default(),
            RateSampler::from_config(1.0, &FitnessConfig::Static(Fitness::Neutral)),
            rng,
        );

        let burden_after = genotype::MutationalBurden::from_moran(&moran).unwrap();
        assert_eq!(
            burden_before.0.keys().sum::<u16>(),
            burden_after.0.keys().sum::<u16>()
        );
        assert_eq!(old_cells, moran.subclones.compute_tot_cells());
    }

    fn the_only_one_subclone_present(subclones: &SubClones) -> Option<CloneId> {
        //! returns `None` if more than one subclone is present else the clone
        //! id.
        //!
        //! **warning** this is very slow.
        let tot_cells = subclones.compute_tot_cells();
        (0..MAX_SUBCLONES as CloneId)
            .find(|&id| subclones.get_clone(id).unwrap().cell_count() == tot_cells)
    }

    #[quickcheck]
    fn subsample_moran_to_same_size(seed: u64) {
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);
        let moran = Moran::default();
        let nb_cells_before = moran.subclones.get_cells().len();
        let cells = nb_cells_before;
        let subsampled = moran.into_subsampled(NonZeroUsize::new(cells).unwrap(), rng);
        assert!(the_only_one_subclone_present(&subsampled.subclones).is_none());
        assert!(nb_cells_before == subsampled.subclones.get_cells().len());
    }

    #[quickcheck]
    fn subsample_moran(cells: NonZeroU8, seed: u64) {
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);
        let moran = Moran::default();
        let nb_cells_before = moran.subclones.get_cells().len();
        let subsampled =
            moran.into_subsampled(NonZeroUsize::new(cells.get() as usize).unwrap(), rng);
        if cells.get() == 1 {
            assert!(the_only_one_subclone_present(&subsampled.subclones).is_some());
        } else {
            assert!(the_only_one_subclone_present(&subsampled.subclones).is_none());
        }
        // this is always true as here we subsample to max 255 cells
        assert!(nb_cells_before > subsampled.subclones.get_cells().len());
    }
}
