use crate::clap_app::{Cli, ProcessType};
use chrono::Utc;
use clap_app::{FitnessArg, Parallel};
use hsc::{
    process::{Exponential, Moran, ProcessOptions},
    stemcell::StemCell,
    subclone::{
        from_mean_std_to_shape_scale, from_shape_scale_to_mean_std, Distributions, Fitness,
        SubClones, Variants,
    },
    write2file,
};
use indicatif::ParallelProgressIterator;
use rand::{Rng, SeedableRng};
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
    fitness: Option<FitnessArg>,
    runs: usize,
    /// division rate for the wild-type
    b0: f32,
    seed: u64,
    parallel: Parallel,
    options_moran: SimulationOptions,
    options_exponential: Option<SimulationOptions>,
    pub snapshots: Vec<f32>,
    pub p_asymmetric: f64,
    pub mu0: Option<f32>,
}

fn create_filename(u: f64, fitness: Fitness, cells: u64, b0: f32, idx: usize) -> PathBuf {
    let (mean, std) = match fitness {
        Fitness::Neutral => (0., 0.),
        Fitness::Fixed { s } => (s, 0.),
        Fitness::GammaSampled { shape, scale } => from_shape_scale_to_mean_std(shape, scale),
    };
    format!(
        "{}u_{}mean_{}std_{}b0_{}cells_{}id",
        u,
        mean,
        std,
        b0,
        cells - 1,
        idx
    )
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

        let fitness = if let Some(fitness) = app.fitness.as_ref() {
            if let Some(s) = fitness.s {
                Fitness::Fixed { s }
            } else if let Some(mean_std) = fitness.mean_std.as_ref() {
                let (shape, scale) = from_mean_std_to_shape_scale(mean_std[0], mean_std[1]);
                Fitness::GammaSampled { shape, scale }
            } else {
                Fitness::Neutral
            }
        } else {
            // If not present, draw values of use mean_std between (mean_min=0.01, std_min=0.01)
            // and (mean_max=0.4, std_max=0.1)
            let (mean, std) = (rng.gen_range(0.01..0.4), rng.gen_range(0.01..0.1));
            let (shape, scale) = from_mean_std_to_shape_scale(mean, std);
            Fitness::GammaSampled { shape, scale }
        };
        let mu0 = if let Some(mu0) = app.mu0 {
            mu0
        } else {
            rng.gen_range(1..12) as f32
        };

        // convert into rates per division
        let process_type = ProcessType::new(app.p_asymmetric);
        let u = Cli::normalise_mutation_rate(
            (mu0 / (app.b0 * app.options_moran.gillespie_options.max_cells as f32)) as f64,
            process_type,
        );
        let distributions = Distributions::new(
            app.p_asymmetric,
            u,
            app.options_moran.gillespie_options.verbosity,
        );

        let rates = subclones.gillespie_rates(&fitness, app.b0, rng);
        let mut snapshots = app.snapshots.clone();
        snapshots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut snapshots = VecDeque::from(snapshots);

        let filename = create_filename(
            u,
            fitness,
            app.options_moran.gillespie_options.max_cells,
            app.b0,
            idx,
        );

        let mut moran = if let Some(options) = app.options_exponential.as_ref() {
            let mut exp = Exponential::new(
                options.process_options.clone(),
                subclones,
                distributions.clone(),
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
                distributions,
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
