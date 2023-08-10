use std::{
    fs,
    path::{Path, PathBuf},
};

use anyhow::Context;
use chrono::Utc;
use clap_app::Parallel;
use hsc::{
    genotype::{Sfs, StatisticsMutations},
    process::{CellDivisionProbabilities, HSCProcess},
    stemcell::{load_cells, StemCell},
    subclone::SubClone,
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

fn find_timepoints(path: &Path) -> Vec<u8> {
    let dir2read = path.join("variant_fraction");
    let mut timepoints: Vec<u8> = fs::read_dir(dir2read)
        .unwrap()
        .map(|f| {
            f.unwrap()
                .file_name()
                .into_string()
                .unwrap()
                .parse::<u8>()
                .unwrap()
        })
        .collect();
    timepoints.sort_unstable();
    timepoints
}

fn main() {
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

        if process.verbosity > 1 {
            println!("saving the SFS for all timepoints");
        }
        let timepoints = find_timepoints(&app.path);
        // last timepoint first to create the stats
        let last_t = timepoints
            .first()
            .expect("found empty timepoints")
            .to_owned() as usize;
        let path2last_t = process
            .make_path(hsc::process::Stats2Save::Genotypes, last_t)
            .unwrap()
            .with_extension("json");
        let path2sfs_last_t = process
            .make_path(hsc::process::Stats2Save::Sfs, last_t)
            .unwrap()
            .with_extension("json");
        let cells = load_cells(&path2last_t).expect("cannot load cells");
        if process.verbosity > 0 {
            println!(
                "computing the sfs for the last timepoint {} with {} cells",
                last_t,
                cells.len()
            );
            if process.verbosity > 1 {
                println!("cells: {:#?}", cells);
            }
        }
        // these stats are the most complete and will be used later on for the
        // other timepoints as well
        let mut stats = StatisticsMutations::from_cells(
            &cells,
            &process.distributions.poisson,
            &mut rng,
            process.verbosity,
        )
        .with_context(|| "cannot construct the stats for the sfs")
        .unwrap();
        Sfs::from_stats(&stats, process.verbosity)
            .expect("cannot create SFS from stats")
            .save(&path2sfs_last_t)
            .unwrap();

        for t in timepoints.into_iter().skip(1) {
            let path2t = process
                .make_path(hsc::process::Stats2Save::Genotypes, t as usize)
                .unwrap()
                .with_extension("json");
            let path2sfs_t = process
                .make_path(hsc::process::Stats2Save::Sfs, t as usize)
                .unwrap()
                .with_extension("json");
            let cells = load_cells(&path2t).expect("cannot load cells");
            // save the SFS for the timepoint loading its cells but using the
            // stats from the oldest timepoint
            Sfs::from_cells(
                &cells,
                &mut stats,
                &process.distributions.poisson,
                &mut rng,
                process.verbosity,
            )
            .expect("cannot create SFS from stats")
            .save(&path2sfs_t)
            .unwrap();
        }
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
