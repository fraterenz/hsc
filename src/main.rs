pub fn main() {
    println!("H");
}
//
//use std::{
//    collections::VecDeque,
//    fs,
//    path::{Path, PathBuf},
//};
//
//use anyhow::Context;
//use chrono::Utc;
//use clap_app::Parallel;
//use hsc::{
//    genotype::{MutationalBurden, NeutralMutationPoisson, Sfs, StatisticsMutations},
//    process::{Exponential, Moran, ProcessOptions, Stats2Save},
//    stemcell::{load_cells, StemCell},
//    subclone::{Fitness, SubClones, Variants},
//    write2file,
//};
//use indicatif::ParallelProgressIterator;
//use rand::SeedableRng;
//use rand_chacha::ChaCha8Rng;
//use rayon::prelude::{IntoParallelIterator, ParallelIterator};
//use sosa::{simulate, CurrentState, Options};
//
//use crate::clap_app::Cli;
//
//pub mod clap_app;
//
//struct Paths2Stats {
//    burden: PathBuf,
//    burden_entropy: PathBuf,
//    genotypes: PathBuf,
//    sfs: PathBuf,
//    sfs_entropy: PathBuf,
//    // variant_fraction: PathBuf,
//}
//
//impl Paths2Stats {
//    fn make_paths(process: &Moran, cells: usize, timepoint: usize) -> anyhow::Result<Self> {
//        let genotypes = process
//            .make_path(Stats2Save::Genotypes, cells, timepoint)?
//            .with_extension("json");
//        let burden = process
//            .make_path(Stats2Save::Burden, cells, timepoint)?
//            .with_extension("json");
//        let burden_entropy = process
//            .make_path(Stats2Save::BurdenEntropy, cells, timepoint)?
//            .with_extension("json");
//        let sfs = process.make_path(Stats2Save::Sfs, cells, timepoint)?;
//        let sfs_entropy = process.make_path(Stats2Save::SfsEntropy, cells, timepoint)?;
//        Ok(Paths2Stats {
//            burden,
//            burden_entropy,
//            genotypes,
//            sfs,
//            sfs_entropy,
//            // variant_fraction,
//        })
//    }
//}
//
//#[derive(Clone, Debug)]
//pub struct SimulationOptions {
//    process_options: ProcessOptions,
//    gillespie_options: Options,
//}
//
//#[derive(Clone, Debug)]
//pub struct AppOptions {
//    fitness: Fitness,
//    runs: usize,
//    /// division rate for the wild-type
//    b0: f32,
//    seed: u64,
//    parallel: Parallel,
//    options_moran: SimulationOptions,
//    options_exponential: Option<SimulationOptions>,
//    pub snapshots: Vec<f32>,
//}
//
//fn find_timepoints(path: &Path, cells: usize) -> Vec<u8> {
//    let dir2read = path
//        .join(format!("{}cells", cells))
//        .join("variant_fraction");
//    let mut timepoints: Vec<u8> = fs::read_dir(dir2read)
//        .unwrap()
//        .map(|f| {
//            f.unwrap()
//                .file_name()
//                .into_string()
//                .unwrap()
//                .parse::<u8>()
//                .unwrap()
//        })
//        .collect();
//    timepoints.sort_unstable();
//    timepoints
//}
//
//fn save_measurements(
//    process: &Moran,
//    poisson: &NeutralMutationPoisson,
//    exp_poisson: Option<&NeutralMutationPoisson>,
//    rng: &mut ChaCha8Rng,
//) -> anyhow::Result<()> {
//    let mut cells2save = vec![process.subclones.compute_tot_cells() as usize];
//    if let Some(subsampling) = process.cells2subsample.as_ref() {
//        for cell in subsampling {
//            cells2save.push(*cell);
//        }
//    };
//
//    // CREATE STATS
//    let timepoints = find_timepoints(&process.path2dir, cells2save[0]);
//    // first timepoint first to create the stats
//    let first_t = timepoints
//        .last()
//        .expect("found empty timepoints")
//        .to_owned() as usize;
//    if process.verbosity > 1 {
//        println!("found the most recent timepoint {}", first_t);
//    }
//    let paths_first_t_genotypes =
//        Paths2Stats::make_paths(process, cells2save[0], first_t)?.genotypes;
//    let cells = load_cells(&paths_first_t_genotypes)
//        .unwrap_or_else(|_| panic!("cannot load cells from {:#?}", paths_first_t_genotypes));
//    if process.verbosity > 0 {
//        println!(
//            "computing the stats for the first timepoint {} with {} cells using Poisson distribution {:#?}",
//            first_t,
//            cells.len(),
//            exp_poisson.unwrap_or(poisson)
//        );
//        if process.verbosity > 1 {
//            println!("cells: {:#?}", cells);
//        }
//    }
//    // create cells with the exponential if present
//    let stats = &mut StatisticsMutations::from_cells(
//        &cells,
//        exp_poisson.unwrap_or_else(|| {
//            if process.verbosity > 1 {
//                println!(
//                    "using Poisson distribution from the Moran process {:#?}",
//                    poisson
//                );
//            }
//            poisson
//        }),
//        rng,
//        process.verbosity,
//    )
//    .with_context(|| "cannot construct the stats for the burden")?;
//    // DEAL WITH ENTROPY
//    // remove all entries that do not satisfy the condition, i.e.
//    // proliferation events that occurred after `process.snapshot_entropy`
//    let stats4entropy = process
//        .proliferation_event_entropy
//        .map(|proliferation_event_entropy| {
//            stats
//                .clone()
//                .from_stats_removing_recent_entries(proliferation_event_entropy, process.verbosity)
//        });
//
//    for cell2save in cells2save {
//        for t in timepoints.iter() {
//            let paths_t = Paths2Stats::make_paths(process, cell2save, *t as usize)?;
//            if let Ok(cells) = load_cells(&paths_t.genotypes) {
//                // save the burden for the timepoint loading its cells but
//                // using the stats from the oldest timepoint
//                MutationalBurden::from_cells_update_stats(
//                    &cells,
//                    stats,
//                    poisson,
//                    rng,
//                    process.verbosity,
//                )
//                .unwrap_or_else(|_| panic!("cannot create burden from stats for timepoint {}", t))
//                .save(&paths_t.burden)?;
//                // save the sfs now that the stats have been updated
//                Sfs::from_cells(&cells, stats, process.verbosity)
//                    .unwrap_or_else(|_| panic!("cannot create SFS from stats for timepoint {}", t))
//                    .save(&paths_t.sfs)?;
//                if let Some(stats4entr) = &stats4entropy {
//                    // DEAL WITH ENTROPY
//                    // dont use `from_cells_update_stats` since we want to keep only the
//                    // mutations that were present at early ages for the entropy
//                    MutationalBurden::from_cells(&cells, stats4entr, process.verbosity)
//                        .expect("cannot create burden from stats for the entropy")
//                        .save(&paths_t.burden_entropy)?;
//                    if process.verbosity > 1 {
//                        println!("saving the sfs entropy");
//                    }
//                    Sfs::from_cells(&cells, stats4entr, process.verbosity)
//                        .expect("cannot create SFS entropy from stats")
//                        .save(&paths_t.sfs_entropy)?;
//                }
//            }
//            stats.save(&process.make_path(
//                Stats2Save::Stats,
//                cell2save,
//                *timepoints.first().unwrap() as usize,
//            )?)?;
//        }
//    }
//    Ok(())
//}
//
//fn main() {
//    let app = Cli::build();
//
//    if app.options_moran.gillespie_options.verbosity > 1 {
//        println!("app: {:#?}", app);
//    }
//
//    println!(
//        "saving variant fraction at timepoints: {:#?}",
//        app.snapshots
//    );
//
//    println!("{} starting simulation", Utc::now());
//
//    let run_simulations = |idx| {
//        // initial state
//        let cells = if app.options_exponential.is_some() {
//            vec![StemCell::new()]
//        } else {
//            vec![StemCell::new(); app.options_moran.gillespie_options.max_cells as usize - 1]
//        };
//        let subclones = SubClones::new(
//            cells,
//            app.options_moran.gillespie_options.max_cells as usize,
//        );
//        let state = &mut CurrentState {
//            population: Variants::variant_counts(&subclones),
//        };
//        let possible_reactions = subclones.gillespie_set_of_reactions();
//
//        let rng = &mut ChaCha8Rng::seed_from_u64(app.seed);
//        rng.set_stream(idx as u64);
//
//        let rates = subclones.gillespie_rates(&app.fitness, app.b0, rng);
//        let mut snapshots = app.snapshots.clone();
//        snapshots.sort_by(|a, b| a.partial_cmp(b).unwrap());
//        let mut snapshots = VecDeque::from(snapshots);
//
//        let mut moran = if let Some(options) = app.options_exponential.as_ref() {
//            let mut exp = Exponential::new(
//                options.process_options.clone(),
//                subclones,
//                idx,
//                options.gillespie_options.verbosity,
//            );
//
//            let stop = simulate(
//                state,
//                &rates,
//                &possible_reactions,
//                &mut exp,
//                &options.gillespie_options,
//                rng,
//            );
//            if options.gillespie_options.verbosity > 1 {
//                println!(
//                    "exponential simulation {} stopped because {:#?}, nb cells {}",
//                    idx,
//                    stop,
//                    exp.subclones.compute_tot_cells()
//                );
//            }
//
//            let timepoint = snapshots.len();
//            snapshots.pop_front();
//            let moran = exp.switch_to_moran(app.options_moran.process_options.clone(), snapshots);
//            moran
//                .save(timepoint, moran.subclones.compute_tot_cells() as usize, rng)
//                .unwrap();
//            moran
//        } else {
//            Moran::new(
//                app.options_moran.process_options.clone(),
//                subclones,
//                app.snapshots.clone(),
//                idx,
//                0.,
//                app.options_moran.gillespie_options.verbosity,
//            )
//        };
//
//        let stop = simulate(
//            state,
//            &rates,
//            &possible_reactions,
//            &mut moran,
//            &app.options_moran.gillespie_options,
//            rng,
//        );
//        if app.options_moran.gillespie_options.verbosity > 1 {
//            println!("Moran simulation {} stopped because {:#?}", idx, stop);
//        }
//
//        if moran.verbosity > 1 {
//            println!("saving the SFS for all timepoints");
//        }
//        save_measurements(
//            &moran,
//            &app.options_moran.process_options.neutral_poisson,
//            app.options_exponential
//                .as_ref()
//                .map(|options| &options.process_options.neutral_poisson),
//            rng,
//        )
//        .expect("cannot save");
//        write2file(
//            &rates.0,
//            &moran.path2dir.join("rates").join(format!("{}.csv", idx)),
//            None,
//            false,
//        )
//        .expect("cannot save the fitness of the subclones");
//    };
//
//    std::process::exit({
//        match app.parallel {
//            Parallel::Debug | Parallel::False => (0..app.runs).for_each(run_simulations),
//            Parallel::True => (0..app.runs)
//                .into_par_iter()
//                .progress_count(app.runs as u64)
//                .for_each(run_simulations),
//        }
//        println!("{} End simulation", Utc::now());
//        0
//    });
//}
