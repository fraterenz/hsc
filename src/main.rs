use std::path::PathBuf;

use chrono::Utc;
use clap_app::Parallel;
use hsc::{
    neutral::StemCell,
    process::{CellDivisionProbabilities, HSCProcess},
    sfs::SubClone,
    MAX_SUBCLONES,
};
use indicatif::ParallelProgressIterator;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
// use rand_distr::{Distribution, Exp};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use sosa::{simulate, CurrentState, Options, ReactionRates};

use crate::clap_app::Cli;

pub mod clap_app;

pub struct SimulationOptions {
    s: f32,
    runs: usize,
    /// division rate for the wild-type
    b0: f32,
    seed: u64,
    parallel: Parallel,
    probabilities: CellDivisionProbabilities,
    path: PathBuf,
    snapshots: Vec<f32>,
    options: Options,
}

// const NITERS: usize = 20;
// const NCELLS: usize = 20;
// const NITERS: usize = 2_000_000;
// const NCELLS: usize = 200_000;
// const RUNS: usize = 1;
//
// // division rate for the wild-type
// const B0: f32 = 1.;
// // avg fit mutations arising in 1 year
// const MU0: f32 = 40.;
// // avg number of neutral mutations per each proliferative event
// const LAMBDA_POISSON: f32 = 1.;
//
// // proliferative advantage conferred by fit mutations
// const S: f32 = 0.15;
// // const LAMBDA_EXP: f32 = 1.;
//
// // probability of getting an asymmetric division per each proliferate event
// const P_ASYMMETRIC: f64 = 0.;
// const SEED: u64 = 26;
// const VERBOSITY: u8 = 1;

fn main() {
    // let mut rng = ChaCha8Rng::seed_from_u64(SEED);
    // let exp = Exp::new(LAMBDA_EXP).unwrap();
    let app = Cli::build();
    let rates = ReactionRates(core::array::from_fn(|i| {
        if i == 0 {
            app.b0
        } else {
            app.b0 * (1. + app.s)
        }
    }));
    if app.options.verbosity > 1 {
        println!("rates: {:#?}", rates);
    }

    if app.options.verbosity > 1 {
        println!(
            "rate of fit variants per cell division {}",
            app.probabilities.p
        );
    }
    // initial state
    let mut subclones: [SubClone; MAX_SUBCLONES] =
        std::array::from_fn(|i| SubClone::new(i, app.options.max_cells as usize));
    for _ in 0..app.options.max_cells - 1 {
        subclones[0].assign_cell(StemCell::new());
    }
    let init_population = core::array::from_fn(|i| subclones[i].cell_count());
    let state = &mut CurrentState {
        population: init_population,
    };
    let possible_reactions = core::array::from_fn(|i| i);

    println!(
        "saving variant fraction at timepoints: {:#?}",
        app.snapshots
    );
    println!("{} starting simulation", Utc::now(),);
    let run_simulations = |idx| {
        let mut rng = ChaCha8Rng::seed_from_u64(app.seed);
        rng.set_stream(idx as u64);

        let mut process = HSCProcess::new(
            app.probabilities.clone(),
            subclones.clone(),
            app.snapshots.clone(),
            app.path.clone(),
            idx,
            0.,
            app.options.verbosity,
        );
        let stop = simulate(
            &mut state.clone(),
            &rates,
            &possible_reactions,
            &mut process,
            &app.options,
            &mut rng,
        );
        if app.options.verbosity > 1 {
            println!("simulation {} stopped because {:#?}", idx, stop);
        }
        process.save(&mut rng).unwrap();
    };

    std::process::exit({
        match app.parallel {
            Parallel::Debug | Parallel::False => (0..app.runs).for_each(run_simulations),
            Parallel::True => (0..app.runs)
                .into_par_iter()
                .progress_count(app.runs as u64)
                .for_each(run_simulations),
        }
        println!("{} End simulation", Utc::now());
        0
    });
}
