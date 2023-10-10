use crate::clap_app::Cli;
use chrono::Utc;
use clap_app::Parallel;
use hsc::{
    process::{Exponential, Moran, ProcessOptions},
    stemcell::StemCell,
    subclone::{from_shape_scale_to_mean_std, Fitness, SubClones, Variants},
    write2file,
};
use indicatif::ParallelProgressIterator;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use sosa::{simulate, CurrentState, Options};
use std::{collections::VecDeque, path::PathBuf};

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

fn create_filename(u: f64, fitness: Fitness, idx: usize) -> PathBuf {
    let (mean, std) = match fitness {
        Fitness::Neutral => (0., 0.),
        Fitness::Fixed { s } => (s, 0.),
        Fitness::GammaSampled { shape, scale } => from_shape_scale_to_mean_std(shape, scale),
    };
    format!("{:.2}u_{:.2}mean_{:.2}std_{}id", u, mean, std, idx)
        .replace('.', "dot")
        .into()
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

        let filename = create_filename(
            app.options_moran.process_options.distributions.u,
            app.fitness,
            idx,
        );
        dbg!(&filename);

        let mut moran = if let Some(options) = app.options_exponential.as_ref() {
            let mut exp = Exponential::new(
                options.process_options.clone(),
                subclones,
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
            let moran = exp.switch_to_moran(
                app.options_moran.process_options.clone(),
                snapshots,
                filename,
            );
            moran
                .save(timepoint, moran.subclones.compute_tot_cells() as usize, rng)
                .unwrap();
            if let Some(subsample) = moran.cells2subsample.as_ref() {
                for cells in subsample {
                    moran.save(timepoint, *cells, rng).unwrap();
                }
            }
            moran
        } else {
            Moran::new(
                app.options_moran.process_options.clone(),
                subclones,
                app.snapshots.clone(),
                0.,
                filename,
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
