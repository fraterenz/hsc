use crate::clap_app::Cli;
use anyhow::Context;
use chrono::Utc;
use clap_app::Parallel;
use hsc::{
    process::{Exponential, Moran, ProcessOptions, SavingCells, SavingOptions, Snapshot},
    proliferation::{Division, NeutralMutations, Proliferation},
    stemcell::StemCell,
    subclone::{Distributions, Fitness, SubClones, Variants},
    write2file,
};
use indicatif::ParallelProgressIterator;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rand_distr::Bernoulli;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use sosa::{simulate, CurrentState, Options, StopReason};
use std::collections::VecDeque;
use uuid::Uuid;

pub mod clap_app;

#[derive(Clone, Debug)]
pub struct SimulationOptionsMoran {
    tau: f32,
    mu_background: f32,
    mu_division: f32,
    asymmetric_prob: f32,
    process_options: ProcessOptions,
    gillespie_options: Options,
    save_sfs_only: bool,
    save_population: bool,
}

#[derive(Clone, Debug)]
pub struct SimulationOptionsExp {
    tau: f32,
    gillespie_options: Options,
    mu_background: f32,
    mu_division: f32,
    asymmetric_prob: f32,
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
    pub mu0: f32,
    pub verbosity: u8,
    background: bool,
}

fn main() {
    let app = Cli::build()
        .with_context(|| "Cannot construct the app".to_string())
        .unwrap();

    if app.verbosity > 1 {
        println!("app: {:#?}", app);
    }

    println!(
        "saving variant fraction at timepoints: {:#?}",
        app.snapshots
    );

    println!("{} starting simulation", Utc::now());

    let run_simulations = |idx| {
        let rng = &mut ChaCha8Rng::seed_from_u64(app.seed);
        rng.set_stream(idx as u64);

        // initial state
        let cells = if app.options_exponential.is_some() {
            // add a neutral mutation such that we have a clonal variant
            let mutation = Uuid::new_v4();
            vec![StemCell::with_mutations(vec![mutation])]
        } else {
            vec![StemCell::new(); app.options_moran.gillespie_options.max_cells as usize - 1]
        };
        let subclones = SubClones::new(
            cells,
            app.options_moran.gillespie_options.max_cells as usize,
        );
        let state = &mut CurrentState {
            population: Variants::variant_counts(&subclones),
        };
        let possible_reactions = subclones.gillespie_set_of_reactions();

        // convert into rates per division
        let u = (app.options_moran.tau * app.mu0
            / (app.options_moran.gillespie_options.max_cells as f32)) as f64;
        let distributions = Distributions::new(
            u,
            app.options_moran.mu_background,
            app.options_moran.mu_division,
            app.options_moran.gillespie_options.verbosity,
        );

        let rates = subclones.gillespie_rates(&app.fitness, 1. / app.options_moran.tau, rng);
        let snapshots = app.snapshots.clone();

        let (mean, std) = app.fitness.get_mean_std();
        let filename = format!(
            "{}mu0_{}u_{}mean_{}std_{}tau_{}cells_{}idx",
            app.mu0,
            u,
            mean,
            std,
            app.options_moran.tau,
            app.options_moran.gillespie_options.max_cells - 1,
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
            // use 1 as we dont care about time during the exp growth phase
            let rates = subclones.gillespie_rates(&app.fitness, 1.0 / options.tau, rng);
            let division = if (options.asymmetric_prob - 0f32).abs() > f32::EPSILON {
                Division::Asymmetric(Bernoulli::new(options.asymmetric_prob as f64).unwrap())
            } else {
                Division::Symmetric
            };
            let proliferation = Proliferation::new(neutral_mutation, division);
            let mut exp = Exponential::new(
                subclones,
                Distributions::new(
                    u, // technically we should use another u norm. by options.tau
                    // but we dont care too much about clone exp. in the exp.
                    // phase for now
                    options.mu_background,
                    options.mu_division,
                    app.options_moran.gillespie_options.verbosity,
                ),
                proliferation,
                options.gillespie_options.verbosity,
            );

            let stop = simulate(
                state,
                &rates,
                &possible_reactions,
                &mut exp,
                &options.gillespie_options,
                rng,
            );
            if options.gillespie_options.verbosity > 1 {
                println!(
                    "exponential simulation {} stopped because {:#?}, nb cells {}",
                    idx,
                    stop,
                    exp.subclones.compute_tot_cells()
                );
            }

            // switch_to_moran start with time 0
            exp.switch_to_moran(
                ProcessOptions {
                    path: app.options_moran.process_options.path.clone(),
                    snapshots,
                },
                distributions,
                filename,
                app.options_moran.save_sfs_only,
                app.options_moran.save_population,
                rng,
            )
        } else {
            let division = if (app.options_moran.asymmetric_prob - 0f32).abs() > f32::EPSILON {
                Division::Asymmetric(
                    Bernoulli::new(app.options_moran.asymmetric_prob as f64).unwrap(),
                )
            } else {
                Division::Symmetric
            };
            let proliferation = Proliferation::new(neutral_mutation, division);
            Moran::new(
                app.options_moran.process_options.clone(),
                subclones,
                0.,
                SavingOptions {
                    filename,
                    save_sfs_only: app.options_moran.save_sfs_only,
                    save_population: app.options_moran.save_population,
                },
                distributions,
                proliferation,
                app.options_moran.gillespie_options.verbosity,
            )
        };

        let stop = simulate(
            state,
            &rates,
            &possible_reactions,
            &mut moran,
            &app.options_moran.gillespie_options,
            rng,
        );
        moran
            .save(
                moran.time,
                &SavingCells::WholePopulation,
                moran.save_sfs_only,
                rng,
            )
            .unwrap();
        match stop {
            StopReason::MaxTimeReached => {
                if app.options_moran.gillespie_options.verbosity > 1 {
                    println!("Moran simulation {} stopped because {:#?}", idx, stop);
                }
            },
            StopReason::MaxItersReached => println!("the simulation stopped earlier than expected because the max number of iterations has been reached"),
            _ => unreachable!("the simulation shouldnt have stopped")
        }

        if moran.verbosity > 1 {
            println!("saving the SFS for all timepoints");
        }
        write2file(
            &rates.0,
            &moran.path2dir.join("rates").join(format!("{}.csv", idx)),
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
        println!("{} End simulation", Utc::now());
        0
    });
}
