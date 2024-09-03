use crate::clap_app::Cli;
use anyhow::Context;
use chrono::Utc;
use clap_app::Parallel;
use hsc::{
    process::{HscDynamics, SaveOptions, SimulationOptions, Snapshot},
    subclone::Fitness,
};
use indicatif::ParallelProgressIterator;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use std::collections::VecDeque;

pub mod clap_app;

#[derive(Clone, Debug)]
pub struct AppOptions {
    fitness: Fitness,
    runs: usize,
    seed: u64,
    parallel: Parallel,
    options_moran: SimulationOptions,
    options_exponential: SimulationOptions,
    options_save: SaveOptions,
    pub snapshots: VecDeque<Snapshot>,
    pub verbosity: u8,
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
        let mut hsc = HscDynamics::new(
            app.verbosity,
            app.options_moran.gillespie_options.max_cells,
            idx,
        );
        hsc.simulate(
            app.options_exponential.clone(),
            app.options_moran.clone(),
            app.options_save.clone(),
            app.fitness.clone(),
            rng,
        )
        .unwrap();
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
