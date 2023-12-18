use crate::clap_app::Cli;
use anyhow::Context;
use chrono::Utc;
use clap_app::{NeutralMutationRate, Parallel};
use hsc::{
    process::{Exponential, Moran, ProcessOptions, SavingOptions, Snapshot},
    stemcell::StemCell,
    subclone::{Distributions, Fitness, SubClones, Variants},
    write2file,
};
use indicatif::ParallelProgressIterator;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use sosa::{simulate, CurrentState, Options, StopReason};
use std::collections::VecDeque;

pub mod clap_app;

#[derive(Clone, Debug)]
pub struct SimulationOptions {
    process_options: ProcessOptions,
    gillespie_options: Options,
    save_sfs_only: bool,
}

#[derive(Clone, Debug)]
pub struct AppOptions {
    fitness: Fitness,
    runs: usize,
    /// intertime division rate for the wild-type measured in years
    pub tau: f32,
    seed: u64,
    parallel: Parallel,
    options_moran: SimulationOptions,
    options_exponential: Option<SimulationOptions>,
    pub snapshots: VecDeque<Snapshot>,
    pub mu0: f32,
    pub verbosity: u8,
    pub neutral_rate: NeutralMutationRate,
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

        // units: [mut/year]
        let m_background = app.neutral_rate.mu_background;
        // units: [mut/div]
        let m_division = app.neutral_rate.mu_division;

        // initial state
        let cells = if app.options_exponential.is_some() {
            vec![StemCell::new()]
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
        let u = (app.tau * app.mu0 / (app.options_moran.gillespie_options.max_cells as f32)) as f64;
        let distributions = Distributions::new(
            u,
            m_background,
            m_division,
            app.options_moran.gillespie_options.verbosity,
        );

        let rates = subclones.gillespie_rates(&app.fitness, 1. / app.tau, rng);
        let snapshots = app.snapshots.clone();

        let (mean, std) = app.fitness.get_mean_std();
        let filename = format!(
            "{}mu0_{}u_{}mean_{}std_{}tau_{}cells_{}idx",
            app.mu0,
            u,
            mean,
            std,
            app.tau,
            app.options_moran.gillespie_options.max_cells - 1,
            idx
        )
        .replace('.', "dot")
        .into();

        let mut moran = if let Some(options) = app.options_exponential.as_ref() {
            let rate = app
                .neutral_rate
                .mu_exp
                .expect("find no exp neutral rate but options_exp");
            // use 1 as we dont care about time during the exp growth phase
            let rates = subclones.gillespie_rates(&app.fitness, 1.0, rng);
            // we assume no background mutation for the exponential growing phase
            let mut exp = Exponential::new(
                subclones,
                Distributions::new(u, rate, rate, app.options_moran.gillespie_options.verbosity),
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
            )
        } else {
            Moran::new(
                app.options_moran.process_options.clone(),
                subclones,
                0.,
                SavingOptions {
                    filename,
                    save_sfs_only: app.options_moran.save_sfs_only,
                },
                distributions,
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
                moran.subclones.get_cells().len(),
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
