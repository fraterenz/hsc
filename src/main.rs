use std::{fs, path::Path};

use anyhow::Context;
use chrono::Utc;
use clap_app::Parallel;
use hsc::{
    genotype::{MutationalBurden, Sfs, StatisticsMutations},
    process::{HSCProcess, ProcessOptions, Stats2Save},
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
    options: Options,
    process_options: ProcessOptions,
    pub snapshots: Vec<f32>,
}

fn find_timepoints(path: &Path, cells: usize) -> Vec<u8> {
    let dir2read = path
        .join(format!("{}cells", cells))
        .join("variant_fraction");
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

fn save_measurements(process: &HSCProcess, rng: &mut ChaCha8Rng) -> anyhow::Result<()> {
    let mut cells2save = vec![process.cells];
    if let Some(subsampling) = process.cells2subsample {
        cells2save.push(subsampling);
    };

    for cell2save in cells2save {
        let timepoints = find_timepoints(&process.path2dir, cell2save);
        // last timepoint first to create the stats
        let last_t = timepoints
            .first()
            .expect("found empty timepoints")
            .to_owned() as usize;
        if process.verbosity > 1 {
            println!("found the most recent timepoint {}", last_t);
        }
        let path2last_t = process
            .make_path(Stats2Save::Genotypes, cell2save, last_t)?
            .with_extension("json");
        let path2burden_last_t = process
            .make_path(Stats2Save::Burden, cell2save, last_t)?
            .with_extension("json");
        let path2sfs_last_t = process.make_path(Stats2Save::Sfs, cell2save, last_t)?;
        let path2sfs_entropy_last_t =
            process.make_path(Stats2Save::SfsEntropy, cell2save, last_t)?;
        let cells = load_cells(&path2last_t)
            .unwrap_or_else(|_| panic!("cannot load cells from {:#?}", path2burden_last_t));
        if process.verbosity > 0 {
            println!(
                "computing the burden for the last timepoint {} with {} cells",
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
            rng,
            process.verbosity,
        )
        .with_context(|| "cannot construct the stats for the burden")?;
        MutationalBurden::from_cells(
            &cells,
            &mut stats,
            &process.distributions.poisson,
            rng,
            process.verbosity,
        )
        .expect("cannot create burden from stats")
        .save(&path2burden_last_t)?;
        Sfs::from_cells(&cells, &stats, process.verbosity)
            .expect("cannot create SFS from stats")
            .save(&path2sfs_last_t)?;
        // remove all entries that do not satisfy the condition, i.e.
        // proliferation events that occurred after `process.snapshot_entropy`
        let stats4entropy = stats.clone().from_stats_removing_recent_entries(
            process
                .proliferation_event_entropy
                .expect("No proliferation_event_entropy found"),
        );
        Sfs::from_cells(&cells, &stats4entropy, process.verbosity)
            .expect("cannot create SFS entropy from stats")
            .save(&path2sfs_entropy_last_t)?;

        for t in timepoints.into_iter().skip(1) {
            let path2t = process
                .make_path(Stats2Save::Genotypes, cell2save, t as usize)?
                .with_extension("json");
            if let Ok(cells) = load_cells(&path2t) {
                let path2burden_t = process
                    .make_path(Stats2Save::Burden, cell2save, t as usize)?
                    .with_extension("json");
                // save the burden for the timepoint loading its cells but
                // using the stats from the oldest timepoint
                MutationalBurden::from_cells(
                    &cells,
                    &mut stats,
                    &process.distributions.poisson,
                    rng,
                    process.verbosity,
                )
                .expect("cannot create SFS from stats")
                .save(&path2burden_t)?;
                // save the sfs now that the stats have been updated
                let path2sfs_t = process.make_path(Stats2Save::Sfs, cell2save, t as usize)?;
                let path2sfs_entropy_t =
                    process.make_path(Stats2Save::SfsEntropy, cell2save, t as usize)?;
                Sfs::from_cells(&cells, &stats, process.verbosity)
                    .expect("cannot create SFS from stats")
                    .save(&path2sfs_t)?;
                Sfs::from_cells(&cells, &stats4entropy, process.verbosity)
                    .expect("cannot create SFS entropy from stats")
                    .save(&path2sfs_entropy_t)?;
            }
        }
        stats.save(&process.make_path(Stats2Save::Stats, cell2save, last_t)?)?;
    }
    Ok(())
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
            app.process_options.probabilities.p
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
            app.process_options.clone(),
            subclones.clone(),
            app.snapshots.clone(),
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
        save_measurements(&process, &mut rng).expect("cannot save");
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
