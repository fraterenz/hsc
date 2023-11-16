use crate::clap_app::{Cli, ProcessType};
use chrono::Utc;
use clap_app::{FitnessArg, NeutralMutationRate, Parallel};
use hsc::{
    process::{Exponential, Moran, ProcessOptions, SavingOptions},
    stemcell::StemCell,
    subclone::{from_mean_std_to_shape_scale, Distributions, Fitness, SubClones, Variants},
    write2file,
};
use indicatif::ParallelProgressIterator;
use rand::{Rng, SeedableRng};
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
    fitness: Option<FitnessArg>,
    runs: usize,
    /// intertime division rate for the wild-type measured in years
    pub tau: Option<f32>,
    seed: u64,
    parallel: Parallel,
    options_moran: SimulationOptions,
    options_exponential: Option<SimulationOptions>,
    pub snapshots: Vec<f32>,
    pub p_asymmetric: f64,
    pub mu0: Option<f32>,
    pub verbosity: u8,
    pub neutral_rate: NeutralMutationRate,
}

fn main() {
    let app = Cli::build();

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

        let tau = if let Some(tau) = app.tau {
            tau
        } else {
            rng.gen_range(0.1f32..10f32)
        };
        let process_type = ProcessType::new(app.p_asymmetric);

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

        let (fitness, mean, std) = if let Some(fitness) = app.fitness.as_ref() {
            if let Some(s) = fitness.s {
                assert!(s < 1., "s should be smaller than 1");
                (Fitness::Fixed { s }, s, 0.)
            } else if let Some(mean_std) = fitness.mean_std.as_ref() {
                assert!(mean_std[0] < 1., "the mean should be smaller than 1");
                assert!(mean_std[1] < 0.1, "the std should be smaller than 0.1");
                let (shape, scale) = from_mean_std_to_shape_scale(mean_std[0], mean_std[1]);
                (
                    Fitness::GammaSampled { shape, scale },
                    mean_std[0],
                    mean_std[1],
                )
            } else {
                (Fitness::Neutral, 0., 0.)
            }
        } else {
            // If not present, draw values of use mean_std between (mean_min=0.01, std_min=0.01)
            // and (mean_max=0.4, std_max=0.1)
            let (mean, std) = (rng.gen_range(0.01..0.4), rng.gen_range(0.005..0.1));
            let (shape, scale) = from_mean_std_to_shape_scale(mean, std);
            (Fitness::GammaSampled { shape, scale }, mean, std)
        };
        let mu0 = if let Some(mu0) = app.mu0 {
            mu0
        } else {
            rng.gen_range(0.1f32..20f32)
        };

        // convert into rates per division
        let u = Cli::normalise_mutation_rate(
            (tau * mu0 / (app.options_moran.gillespie_options.max_cells as f32)) as f64,
            process_type,
        );
        let distributions = Distributions::new(
            app.p_asymmetric,
            u,
            m_background,
            m_division,
            app.options_moran.gillespie_options.verbosity,
        );

        let rates = subclones.gillespie_rates(&fitness, 1. / tau, rng);
        let mut snapshots = app.snapshots.clone();
        snapshots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut snapshots = VecDeque::from(snapshots);

        let filename = format!(
            "{}mu0_{}u_{}mean_{}std_{}tau_{}cells_{}idx",
            mu0,
            u,
            mean,
            std,
            tau,
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
            let rates = subclones.gillespie_rates(&fitness, 1.0, rng);
            // since we are using 1, we dont divide the rate by r
            let m = Cli::normalise_mutation_rate(rate, process_type);
            // we assume no background mutation for the exponential growing phase
            let mut exp = Exponential::new(
                subclones,
                Distributions::new(
                    app.p_asymmetric,
                    u,
                    m,
                    m,
                    app.options_moran.gillespie_options.verbosity,
                ),
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

            snapshots.pop_front();
            // switch_to_moran start with time 0
            let moran = exp.switch_to_moran(
                app.options_moran.process_options.clone(),
                snapshots,
                distributions,
                filename,
                app.options_moran.save_sfs_only,
            );
            moran
                .save(
                    moran.time,
                    moran.subclones.compute_tot_cells() as usize,
                    app.options_moran.save_sfs_only,
                    rng,
                )
                .unwrap();
            if let Some(subsample) = moran.cells2subsample.as_ref() {
                for cells in subsample {
                    moran
                        .save(moran.time, *cells, app.options_moran.save_sfs_only, rng)
                        .unwrap();
                }
            }
            moran
        } else {
            Moran::new(
                app.options_moran.process_options.clone(),
                subclones,
                app.snapshots.clone(),
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
