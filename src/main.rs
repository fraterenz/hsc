use std::path::PathBuf;

use chrono::Utc;
use hsc::{neutral::StemCell, process::HSCProcess, sfs::SubClone, MAX_SUBCLONES};
use indicatif::ParallelProgressIterator;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
// use rand_distr::{Distribution, Exp};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use sosa::{simulate, CurrentState, Options, ReactionRates};

// const NITERS: usize = 20;
// const NCELLS: usize = 20;
const NITERS: usize = 2_000_000;
const NCELLS: usize = 200_000;
const RUNS: usize = 1;

// division rate for the wild-type
const B0: f32 = 1.;
// avg fit mutations arising in 1 year
const MU0: f32 = 40.;
// avg number of neutral mutations per each proliferative event
const LAMBDA_POISSON: f32 = 1.;

// proliferative advantage conferred by fit mutations
const S: f32 = 0.15;
// const LAMBDA_EXP: f32 = 1.;

// probability of getting an asymmetric division per each proliferate event
const P_ASYMMETRIC: f64 = 0.;
const SEED: u64 = 26;
const VERBOSITY: u8 = 1;

fn main() {
    // let mut rng = ChaCha8Rng::seed_from_u64(SEED);
    // let exp = Exp::new(LAMBDA_EXP).unwrap();
    let rates = ReactionRates(core::array::from_fn(|i| {
        if i == 0 {
            B0
        } else {
            B0 * (1. + S)
        }
    }));
    if VERBOSITY > 1 {
        println!("rates: {:#?}", rates);
    }
    // rate of fit variant per cell division
    let u = if (P_ASYMMETRIC - 0.).abs() <= f64::EPSILON {
        // in the symmetric case we need to divide by 2 because 1 division
        // corresponds to two cells, hence two draws (i.e. two mutational
        // events) per division
        (MU0 / (2. * B0 * NCELLS as f32)) as f64
    } else {
        (MU0 / (B0 * NCELLS as f32)) as f64
    };
    if VERBOSITY > 1 {
        println!("rate of fit variants per cell division {}", u);
    }

    // initial state
    let mut subclones: [SubClone; MAX_SUBCLONES] =
        std::array::from_fn(|i| SubClone::new(i, NCELLS));
    for _ in 0..NCELLS {
        subclones[0].assign_cell(StemCell::new());
    }
    let init_population = core::array::from_fn(|i| subclones[i].cell_count());
    let state = &mut CurrentState {
        population: init_population,
    };

    let options = Options {
        max_iter: NITERS,
        max_cells: (NCELLS + 1) as u64,
        init_iter: 0,
        verbosity: VERBOSITY,
    };

    let possible_reactions = core::array::from_fn(|i| i);
    println!("{} starting simulation", Utc::now(),);

    (0..RUNS)
        // .into_iter()
        .into_par_iter()
        .progress_count(RUNS as u64)
        .for_each(|idx| {
            let mut rng = ChaCha8Rng::seed_from_u64(SEED);
            rng.set_stream(idx as u64);

            let mut process = HSCProcess::new(
                P_ASYMMETRIC,
                LAMBDA_POISSON,
                u,
                subclones.clone(),
                Some(vec![100, 500, NITERS - 1]),
                PathBuf::from("./test"),
                idx,
                vec![0.],
                VERBOSITY,
            );
            let stop = simulate(
                &mut state.clone(),
                &rates,
                &possible_reactions,
                &mut process,
                &options,
                &mut rng,
            );
            if VERBOSITY > 1 {
                println!("simulation {} stopped because {:#?}", idx, stop);
            }
            process.save(&mut rng).unwrap();
        });
    println!("{} end simulation", Utc::now(),);
}
