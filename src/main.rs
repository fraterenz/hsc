use std::{
    fs,
    path::{Path, PathBuf},
};

use anyhow::Context;
use chrono::Utc;
use clap_app::Parallel;
use hsc::{
    genotype::{MutationalBurden, Sfs, StatisticsMutations},
    process::{Exponential, Moran, ProcessOptions, Stats2Save},
    stemcell::{load_cells, StemCell},
    subclone::{Fitness, SubClones, Variants},
    write2file,
};
use indicatif::ParallelProgressIterator;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
// use rand_distr::{Distribution, Exp};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use sosa::{simulate, CurrentState, Options};

use crate::clap_app::Cli;

pub mod clap_app;

struct Paths2Stats {
    burden: PathBuf,
    burden_entropy: PathBuf,
    genotypes: PathBuf,
    sfs: PathBuf,
    sfs_entropy: PathBuf,
    // variant_fraction: PathBuf,
}

impl Paths2Stats {
    fn make_paths(process: &Moran, cells: usize, timepoint: usize) -> anyhow::Result<Self> {
        let genotypes = process
            .make_path(Stats2Save::Genotypes, cells, timepoint)?
            .with_extension("json");
        let burden = process
            .make_path(Stats2Save::Burden, cells, timepoint)?
            .with_extension("json");
        let burden_entropy = process
            .make_path(Stats2Save::BurdenEntropy, cells, timepoint)?
            .with_extension("json");
        let sfs = process.make_path(Stats2Save::Sfs, cells, timepoint)?;
        let sfs_entropy = process.make_path(Stats2Save::SfsEntropy, cells, timepoint)?;
        Ok(Paths2Stats {
            burden,
            burden_entropy,
            genotypes,
            sfs,
            sfs_entropy,
            // variant_fraction,
        })
    }
}

#[derive(Clone, Debug)]
pub struct SimulationOptions {
    fitness: Fitness,
    runs: usize,
    /// division rate for the wild-type
    b0: f32,
    seed: u64,
    parallel: Parallel,
    options_moran: Options,
    options_exponential: Options,
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

fn save_measurements(process: &Moran, rng: &mut ChaCha8Rng) -> anyhow::Result<()> {
    let mut cells2save = vec![process.subclones.compute_tot_cells() as usize];
    if let Some(subsampling) = process.cells2subsample.as_ref() {
        for cell in subsampling {
            cells2save.push(*cell);
        }
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
        let paths_last_t_genotypes = Paths2Stats::make_paths(process, cell2save, last_t)?.genotypes;
        let cells = load_cells(&paths_last_t_genotypes)
            .unwrap_or_else(|_| panic!("cannot load cells from {:#?}", paths_last_t_genotypes));
        if process.verbosity > 0 {
            println!(
                "computing the stats for the last timepoint {} with {} cells",
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

        for t in timepoints.iter() {
            let paths_t = Paths2Stats::make_paths(process, cell2save, *t as usize)?;
            if let Ok(cells) = load_cells(&paths_t.genotypes) {
                // save the burden for the timepoint loading its cells but
                // using the stats from the oldest timepoint
                MutationalBurden::from_cells_update_stats(
                    &cells,
                    &mut stats,
                    &process.distributions.poisson,
                    rng,
                    process.verbosity,
                )
                .unwrap_or_else(|_| panic!("cannot create burden from stats for timepoint {}", t))
                .save(&paths_t.burden)?;
                // save the sfs now that the stats have been updated
                Sfs::from_cells(&cells, &stats, process.verbosity)
                    .unwrap_or_else(|_| panic!("cannot create SFS from stats for timepoint {}", t))
                    .save(&paths_t.sfs)?;
            }
        }
        stats.save(&process.make_path(
            Stats2Save::Stats,
            cell2save,
            *timepoints.first().unwrap() as usize,
        )?)?;

        // DEAL WITH ENTROPY
        // remove all entries that do not satisfy the condition, i.e.
        // proliferation events that occurred after `process.snapshot_entropy`
        let stats4entropy = stats.from_stats_removing_recent_entries(
            process
                .proliferation_event_entropy
                .expect("No proliferation_event_entropy found"),
            process.verbosity,
        );
        for t in timepoints.into_iter() {
            let paths_t = Paths2Stats::make_paths(process, cell2save, t as usize)?;
            if let Ok(cells) = load_cells(&paths_t.genotypes) {
                // dont use `from_cells_update_stats` since we want to keep only the
                // mutations that were present at early ages for the entropy
                MutationalBurden::from_cells(&cells, &stats4entropy, process.verbosity)
                    .expect("cannot create burden from stats for the entropy")
                    .save(&paths_t.burden_entropy)?;
                if process.verbosity > 1 {
                    println!("saving the sfs entropy");
                }
                Sfs::from_cells(&cells, &stats4entropy, process.verbosity)
                    .expect("cannot create SFS entropy from stats")
                    .save(&paths_t.sfs_entropy)?;
            }
        }
    }
    Ok(())
}

fn main() {
    let app = Cli::build();

    if app.options_moran.verbosity > 1 {
        println!(
            "rate of fit variants per cell division {}",
            app.process_options.probabilities.p
        );
        println!("app: {:#?}", app);
    }

    println!(
        "saving variant fraction at timepoints: {:#?}",
        app.snapshots
    );

    println!("{} starting simulation", Utc::now());

    let run_simulations = |idx| {
        // initial state
        let cells = if app.process_options.exponential {
            vec![StemCell::new()]
        } else {
            vec![StemCell::new(); app.options_moran.max_cells as usize - 1]
        };
        let subclones = SubClones::new(cells, app.options_moran.max_cells as usize);
        let state = &mut CurrentState {
            population: Variants::variant_counts(&subclones),
        };
        let possible_reactions = subclones.gillespie_set_of_reactions();

        let rng = &mut ChaCha8Rng::seed_from_u64(app.seed);
        rng.set_stream(idx as u64);

        // let rates = subclones.gillespie_rates(app.fitness, app.b0);
        let rates = subclones.gillespie_rates(&app.fitness, app.b0, rng);

        let mut moran = if app.process_options.exponential {
            let mut exp = Exponential::new(
                app.process_options.clone(),
                subclones,
                idx,
                app.options_exponential.verbosity,
            );

            let stop = simulate(
                state,
                &rates,
                &possible_reactions,
                &mut exp,
                &app.options_exponential,
                rng,
            );
            if app.options_exponential.verbosity > 1 {
                println!(
                    "exponential simulation {} stopped because {:#?}, nb cells {}",
                    idx,
                    stop,
                    exp.subclones.compute_tot_cells()
                );
            }

            exp.switch_to_moran(app.process_options.clone(), app.snapshots.clone())
        } else {
            Moran::new(
                app.process_options.clone(),
                subclones,
                app.snapshots.clone(),
                idx,
                0.,
                app.options_moran.verbosity,
            )
        };

        let stop = simulate(
            state,
            &rates,
            &possible_reactions,
            &mut moran,
            &app.options_moran,
            rng,
        );
        if app.options_moran.verbosity > 1 {
            println!("Moran simulation {} stopped because {:#?}", idx, stop);
        }

        if moran.verbosity > 1 {
            println!("saving the SFS for all timepoints");
        }
        save_measurements(&moran, rng).expect("cannot save");
        write2file(
            &rates.0,
            &moran.path2dir.join("rates").join(format!("{}.csv", idx)),
            None,
            false,
        )
        .expect("cannot save the fitness of the subclones");
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
