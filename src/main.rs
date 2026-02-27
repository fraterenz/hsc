use crate::clap_app::Cli;
use anyhow::Context;
use clap_app::Parallel;
use hsc::{
    process::{switch_to_moran, Exponential, Moran, ProcessOptions},
    proliferation::{Division, NeutralMutations, Proliferation},
    snapshots::{SavingCells, SavingOptions, Snapshot, StatsConfig},
    stemcell::StemCell,
    subclone::{Distributions, Fitness, SubClones, Variants},
    write2file, Probs, ProbsPerYear,
};
use indicatif::ParallelProgressIterator;
use log::{debug, info};
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rand_distr::Bernoulli;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use sosa::{simulate, CurrentState, Options, StopReason};
use std::{collections::VecDeque, num::NonZeroUsize, path::PathBuf};
use uuid::Uuid;

pub mod clap_app;

#[derive(Clone, Debug)]
pub struct OptionsTreatment {
    before_treatment: Options,
    after_treatment: Options,
    cells_left: NonZeroUsize,
    regrowth_options: SimulationOptionsExp,
}

#[derive(Clone, Debug)]
pub enum GillespieOptions {
    NoTreatment(Options),
    Treatment(OptionsTreatment),
}

#[derive(Clone, Debug)]
pub struct SimulationOptionsMoran {
    tau: f32,
    probs_per_year: ProbsPerYear,
    process_options: ProcessOptions,
    gillespie_options: GillespieOptions,
    stats: StatsConfig,
    save_population: bool,
    asymmetric: f32,
}

#[derive(Clone, Debug)]
pub struct SimulationOptionsExp {
    tau: f32,
    gillespie_options: Options,
    probs_per_year: ProbsPerYear,
    asymmetric: f32,
}

#[derive(Clone, Debug)]
pub struct AppOptions {
    fitness: Fitness,
    runs: usize,
    seed: u64,
    parallel: Parallel,
    options_moran: SimulationOptionsMoran,
    options_exponential: Option<SimulationOptionsExp>,
    pub snapshots: VecDeque<Snapshot>,
    background: bool,
}

fn main() {
    env_logger::init();
    let app = Cli::build()
        .with_context(|| "Cannot construct the app".to_string())
        .unwrap();

    debug!("app: {app:#?}");

    info!("saving process at timepoints: {:#?}", app.snapshots);
    info!("starting simulations (tot runs: {})", app.runs);

    let run_simulations = |idx| {
        let rng = &mut ChaCha8Rng::seed_from_u64(app.seed);
        rng.set_stream(idx as u64);

        let options_moran_gillespie = match &app.options_moran.gillespie_options {
            GillespieOptions::Treatment(op) => &op.before_treatment,
            GillespieOptions::NoTreatment(op) => op,
        };

        // initial state
        let cells = if app.options_exponential.is_some() {
            // add a neutral mutation such that we have a clonal variant
            let mutation = Uuid::new_v4();
            vec![StemCell::with_mutations(vec![mutation])]
        } else {
            vec![StemCell::new(); options_moran_gillespie.max_cells as usize - 1]
        };
        let subclones = SubClones::new(cells, options_moran_gillespie.max_cells as usize - 1);
        let state = &mut CurrentState {
            population: Variants::variant_counts(&subclones),
        };
        let possible_reactions = subclones.array_of_gillespie_reactions();

        let rates_moran = subclones.gillespie_rates(&app.fitness, 1. / app.options_moran.tau, rng);
        let snapshots = app.snapshots.clone();

        let (mean, std) = app.fitness.get_mean_std();
        let filename: PathBuf = format!(
            "{}mu_{}mean_{}std_{}tau_{}cells_{}idx",
            app.options_moran.probs_per_year.mu,
            mean,
            std,
            app.options_moran.tau,
            options_moran_gillespie.max_cells - 1,
            idx
        )
        .replace('.', "dot")
        .into();
        let neutral_mutation = if app.background {
            NeutralMutations::UponDivisionAndBackground
        } else {
            NeutralMutations::UponDivision
        };

        let mut moran = if let Some(options) = app.options_exponential.as_ref() {
            // convert into rates per division
            let probs_exp = Probs::new(
                options.probs_per_year.mu_background,
                options.probs_per_year.mu_division,
                options.probs_per_year.mu,
                options.asymmetric,
                options.gillespie_options.max_cells,
            );
            let distributions_exp = Distributions::new(probs_exp);
            let rates = subclones.gillespie_rates(&app.fitness, 1.0 / options.tau, rng);
            let division = if (options.asymmetric - 0.).abs() < f32::EPSILON {
                Division::Symmetric
            } else {
                Division::Asymmetric(Bernoulli::new(options.asymmetric as f64).unwrap())
            };
            let proliferation = Proliferation::new(neutral_mutation, division);
            let mut exp = Exponential::new(
                subclones,
                distributions_exp,
                proliferation,
                0.,
                hsc::process::ExponentialPhase::Development,
            );
            debug!("start simulating exp. phase");

            let stop = simulate(
                state,
                &rates,
                &possible_reactions,
                &mut exp,
                &options.gillespie_options,
                rng,
            );
            debug!(
                "exponential simulation {} stopped because {:#?}, nb cells {}",
                idx,
                stop,
                exp.subclones.compute_tot_cells()
            );

            // convert into rates per division
            let probs_moran = Probs::new(
                app.options_moran.probs_per_year.mu_background,
                app.options_moran.probs_per_year.mu_division,
                app.options_moran.probs_per_year.mu,
                app.options_moran.asymmetric,
                options_moran_gillespie.max_cells - 1,
            );
            let moran_distributions = Distributions::new(probs_moran);
            debug!("switching to moran");
            // switch_to_moran start with time 0
            let mut moran = switch_to_moran(
                exp,
                ProcessOptions {
                    path: app.options_moran.process_options.path.clone(),
                    snapshots,
                },
                moran_distributions,
                filename.clone(),
                app.options_moran.stats,
                app.options_moran.save_population,
                rng,
            );
            moran
                .save(moran.time, &SavingCells::WholePopulation, rng)
                .expect("cannot save simulation at the end of the exponential phase");
            moran
        } else {
            let division = if (app.options_moran.asymmetric - 0.).abs() < f32::EPSILON {
                Division::Symmetric
            } else {
                Division::Asymmetric(Bernoulli::new(app.options_moran.asymmetric as f64).unwrap())
            };
            // convert into rates per division
            let probs_moran = Probs::new(
                app.options_moran.probs_per_year.mu_background,
                app.options_moran.probs_per_year.mu_division,
                app.options_moran.probs_per_year.mu,
                app.options_moran.asymmetric,
                options_moran_gillespie.max_cells - 1,
            );
            let moran_distributions = Distributions::new(probs_moran);
            let proliferation = Proliferation::new(neutral_mutation, division);
            Moran::new(
                app.options_moran.process_options.clone(),
                subclones,
                0.,
                SavingOptions {
                    filename: filename.clone(),
                    save_population: app.options_moran.save_population,
                },
                app.options_moran.stats,
                moran_distributions,
                proliferation,
            )
        };
        debug!("simulating Moran phase");

        let stop = simulate(
            state,
            &rates_moran,
            &possible_reactions,
            &mut moran,
            options_moran_gillespie,
            rng,
        );
        moran
            .save(moran.time, &SavingCells::WholePopulation, rng)
            .unwrap();

        // treatment
        let (moran, stop) = if let GillespieOptions::Treatment(op) =
            &app.options_moran.gillespie_options
        {
            let (moran_distributions, snapshots) =
                (moran.distributions.clone(), moran.snapshots.clone());

            // treatment subsamples the total population
            debug!(
                "Subsample Moran at {} cells before starting to regrowth",
                op.cells_left
            );
            let mut regrowth: Exponential = moran.into_subsampled(op.cells_left, rng).into();
            debug!(
                "restart exp growth with time {} from {} cells with options {:#?}",
                regrowth.time,
                regrowth.subclones.get_cells().len(),
                &op.regrowth_options.gillespie_options,
            );
            let state = &mut CurrentState {
                population: Variants::variant_counts(&regrowth.subclones),
            };
            let rates =
                regrowth
                    .subclones
                    .gillespie_rates(&app.fitness, 1. / op.regrowth_options.tau, rng);

            let stop = simulate(
                state,
                &rates,
                &possible_reactions,
                &mut regrowth,
                &op.regrowth_options.gillespie_options,
                rng,
            );
            debug!(
                "regrowing simulation {} stopped at time {} because {:#?} with nb cells {}",
                idx,
                regrowth.time,
                stop,
                regrowth.subclones.compute_tot_cells()
            );
            let mut moran = switch_to_moran(
                regrowth,
                ProcessOptions {
                    path: app.options_moran.process_options.path.clone(),
                    snapshots,
                },
                moran_distributions,
                filename,
                app.options_moran.stats,
                app.options_moran.save_population,
                rng,
            );
            // hack required because simulate always starts at time 0, this is
            // time relative to where the simulation is
            let mut after_treatment_relative = op.after_treatment.clone();
            after_treatment_relative.max_iter_time.time =
                op.after_treatment.max_iter_time.time - moran.time;
            debug!(
                "simulating Moran phase at time {} with {:#?}",
                moran.time, after_treatment_relative,
            );

            let stop = simulate(
                state,
                &rates_moran,
                &possible_reactions,
                &mut moran,
                &after_treatment_relative,
                rng,
            );
            debug!(
                "end of Moran phase due to {:#?} at time {}, saving now",
                stop, moran.time
            );
            moran
                .save(moran.time, &SavingCells::WholePopulation, rng)
                .unwrap();
            (moran, stop)
        } else {
            (moran, stop)
        };
        match stop {
            StopReason::MaxTimeReached => {
                    debug!("Moran simulation {idx} stopped because {stop:#?}");
            },
            StopReason::MaxItersReached => debug!("the simulation stopped earlier than expected because the max number of iterations has been reached"),
            _ => unreachable!("the simulation shouldnt have stopped")
        }

        debug!("saving the SFS for all timepoints");
        write2file(
            &rates_moran.0,
            &moran.path2dir.join("rates").join(format!("{idx}.csv")),
            None,
            false,
        )
        .expect("cannot save the fitness of the subclones");
    };

    std::process::exit({
        // start from seed for the array job
        let start = (app.seed * 10) as usize;
        let start_end = start..app.runs + start;

        match app.parallel {
            Parallel::Debug | Parallel::False => (start_end).for_each(run_simulations),
            Parallel::True => (start_end)
                .into_par_iter()
                .progress_count(app.runs as u64)
                .for_each(run_simulations),
        }
        info!("end simulation");
        0
    });
}
