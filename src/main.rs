use crate::clap_app::Cli;
use chrono::Utc;
use clap_app::Parallel;
use hsc::{
    process::{Exponential, Moran, ProcessOptions},
    stemcell::StemCell,
    subclone::{Fitness, SubClones, Variants},
    write2file,
};
use indicatif::ParallelProgressIterator;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use sosa::{simulate, CurrentState, Options};
use std::collections::VecDeque;

pub mod clap_app;

#[derive(Clone, Debug)]
pub struct SimulationOptions {
    process_options: ProcessOptions,
    gillespie_options: Options,
}

#[derive(Clone, Debug)]
pub struct AppOptions {
    fitness: Fitness,
    runs: usize,
    /// division rate for the wild-type
    b0: f32,
    seed: u64,
    parallel: Parallel,
    options_moran: SimulationOptions,
    options_exponential: Option<SimulationOptions>,
    pub snapshots: Vec<f32>,
}

fn main() {
    let app = Cli::build();

    if app.options_moran.gillespie_options.verbosity > 1 {
        println!("app: {:#?}", app);
    }

    println!(
        "saving variant fraction at timepoints: {:#?}",
        app.snapshots
    );

    println!("{} starting simulation", Utc::now());

    let run_simulations = |idx| {
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

        let rng = &mut ChaCha8Rng::seed_from_u64(app.seed);
        rng.set_stream(idx as u64);

        let rates = subclones.gillespie_rates(&app.fitness, app.b0, rng);
        let mut snapshots = app.snapshots.clone();
        snapshots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut snapshots = VecDeque::from(snapshots);

        let mut moran = if let Some(options) = app.options_exponential.as_ref() {
            let mut exp = Exponential::new(
                options.process_options.clone(),
                subclones,
                idx,
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

            let timepoint = snapshots.len();
            snapshots.pop_front();
            let moran = exp.switch_to_moran(app.options_moran.process_options.clone(), snapshots);
            moran
                .save(timepoint, moran.subclones.compute_tot_cells() as usize, rng)
                .unwrap();
            moran
        } else {
            Moran::new(
                app.options_moran.process_options.clone(),
                subclones,
                app.snapshots.clone(),
                idx,
                0.,
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
        if app.options_moran.gillespie_options.verbosity > 1 {
            println!("Moran simulation {} stopped because {:#?}", idx, stop);
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
        let start_end = app.seed as usize..app.runs + app.seed as usize;

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
