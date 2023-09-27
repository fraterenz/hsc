use clap::{ArgAction, Args, Parser};
use hsc::{genotype::NeutralMutationPoisson, process::ProcessOptions, subclone::Distributions};
use num_traits::{Float, NumCast};
use sosa::{IterTime, NbIndividuals, Options};
use std::path::PathBuf;

use crate::{AppOptions, Fitness, SimulationOptions};

#[derive(Clone, Copy, Debug)]
enum ProcessType {
    MoranSymmetric,
    MoranAsymmetric,
}

impl ProcessType {
    pub fn new(p_asymmetric: f64) -> Self {
        if (p_asymmetric - 0.).abs() < f64::EPSILON {
            ProcessType::MoranSymmetric
        } else {
            ProcessType::MoranAsymmetric
        }
    }
}

#[derive(Clone, Debug)]
pub enum Parallel {
    False,
    True,
    Debug,
}

#[derive(Args, Debug, Clone)]
#[group(required = true, multiple = false)]
struct FitnessArg {
    /// Neutral scenario with all clones having the same fitness coefficient of
    /// 0
    #[arg(long, action = ArgAction::SetTrue)]
    neutral: Option<bool>,
    /// proliferative advantage conferred by fit mutations assuming all clones
    /// have the same advantange, units: mutation / cell
    #[arg(long)]
    s: Option<f32>,
    // /// The Gamma distribution used to sample the fitness coefficients
    // /// representing the proliferative advantage conferred by fit mutations
    // #[command(flatten)]
    // gamma: Option<GammaFitness>,
    /// The mean and the standard deviation of the Gamma distribution used to
    /// sample the fitness coefficients representing the proliferative
    /// advantage conferred by fit mutations
    #[arg(long, num_args = 2)]
    mean_std: Option<Vec<f32>>,
}

#[derive(Args, Debug, Clone)]
#[group(required = true, multiple = true)]
struct NeutralMutationRate {
    /// Poisson rate for the exponential growing phase, leave empty to simulate
    /// just a Moran process with fixed population size
    #[arg(long)]
    mu_exp: Option<f32>,
    /// Background mutation rate (all neutral mutations **not** occuring in the
    /// mitotic phase) for for the constant population phase
    #[arg(long)]
    mu_background: f32,
    /// Division mutation rate (all neutral mutations occuring upon
    /// cell-division, mitotic phase) for for the constant population phase
    #[arg(long)]
    mu_division: f32,
}

#[derive(Debug, Parser)] // requires `derive` feature
#[command(name = "hsc", version, about = "TODO")]
pub struct Cli {
    /// The years simulated will be `years + 1`, that is simulations start at
    /// year zero and the end at year `year + 1`
    #[arg(short, long, default_value_t = 1)]
    years: usize,
    /// Number of cells to simulate
    #[arg(short, long, default_value_t = 100)]
    cells: NbIndividuals,
    #[arg(short, long, default_value_t = 1)]
    runs: usize,
    /// division rate for the wild-type in 1 year, units: division / (year * cell)
    #[arg(long, default_value_t = 1.)]
    b0: f32,
    /// avg fit mutations arising in 1 year, units: division / year
    #[arg(long, default_value_t = 2.)]
    mu0: f32,
    /// avg number of neutral mutations per each proliferative event assuming a
    /// Poisson distribution, units: mutation / (cell * year)
    #[command(flatten)]
    neutral_rate: NeutralMutationRate,
    #[command(flatten)]
    fitness: FitnessArg,
    /// probability of getting an asymmetric division per each proliferate event
    #[arg(long, default_value_t = 0.)]
    p_asymmetric: f64,
    #[arg(long, default_value_t = 26)]
    seed: u64,
    /// Triggers debug mode: max verbosity, 1 sequential simulation, 10 cells,
    /// 20 iterations with high mutation rate
    #[arg(short, long, action = ArgAction::SetTrue, default_value_t = false)]
    debug: bool,
    #[arg(short, long, action = ArgAction::Count, conflicts_with = "debug", default_value_t = 0)]
    verbosity: u8,
    /// Path to store the results of the simulations
    #[arg(value_name = "DIR")]
    path: PathBuf,
    /// Run sequentially each run instead of using rayon for parallelisation
    #[arg(long, action = ArgAction::SetTrue, default_value_t = false, conflicts_with = "debug")]
    sequential: bool,
    #[arg(long, default_value_t = 10)]
    /// Number of snapshots to take to save the simulation. Those timepoints
    /// will be linespaced, starting from 1 to `years`.
    snapshots: u8,
    /// Number of cells to subsample before saving the measurements.
    /// If not specified, do not subsample. If subsampling is performed, the
    /// measurements of the whole population will also be saved.
    #[arg(long, num_args = 0..)]
    subsample: Option<Vec<usize>>,
}

impl Cli {
    fn build_snapshots_from_time(n_snapshots: usize, time: f32) -> Vec<f32> {
        let dx = time / ((n_snapshots - 1) as f32);
        let mut x = vec![0.; n_snapshots];
        for i in 1..n_snapshots - 1 {
            x[i] = x[i - 1] + dx;
        }

        x.shrink_to_fit();
        x[n_snapshots - 1] = time;
        x
    }

    fn normalise_mutation_rate<F: Float>(rate: F, process_type: ProcessType) -> F {
        match process_type {
            ProcessType::MoranAsymmetric => rate,
            ProcessType::MoranSymmetric => rate / NumCast::from(2).unwrap(),
        }
    }

    pub fn build() -> AppOptions {
        let cli = Cli::parse();

        let fitness = if let Some(s) = cli.fitness.s {
            Fitness::Fixed { s }
        } else if let Some(mean_std) = cli.fitness.mean_std {
            Fitness::GammaSampled {
                shape: mean_std[0].powf(2.) / mean_std[1].powf(2.),
                scale: mean_std[1].powf(2.) / mean_std[0],
            }
        } else {
            Fitness::Neutral
        };

        let (parallel, runs) = if cli.debug {
            (Parallel::Debug, 1)
        } else if cli.sequential {
            (Parallel::False, cli.runs)
        } else {
            (Parallel::True, cli.runs)
        };

        let years = cli.years;
        let (max_cells, years, b0, mu0, verbosity) = if cli.debug {
            (11, 1, 1., 4., u8::MAX)
        } else {
            (cli.cells + 1, years, cli.b0, cli.mu0, cli.verbosity)
        };

        let max_iter = 2 * max_cells as usize * b0 as usize * years;
        let snapshots = Cli::build_snapshots_from_time(cli.snapshots as usize, years as f32);

        let process_type = ProcessType::new(cli.p_asymmetric);

        // convert into rates per division
        let u = Cli::normalise_mutation_rate((mu0 / (b0 * max_cells as f32)) as f64, process_type);
        let m_background = cli.neutral_rate.mu_background / b0;
        let m_division = cli.neutral_rate.mu_division / b0;

        let distributions = Distributions::new(cli.p_asymmetric, u, verbosity);

        // Moran
        let process_options = ProcessOptions {
            distributions,
            path: cli.path.clone(),
            cells2subsample: cli.subsample.clone(),
            neutral_poisson: NeutralMutationPoisson::new(m_division, m_background)
                .expect("wrong lambda"), // TODO
        };
        let options_moran = SimulationOptions {
            process_options,
            gillespie_options: Options {
                max_iter_time: IterTime {
                    iter: max_iter,
                    time: years as f32,
                },
                max_cells,
                init_iter: 0,
                verbosity,
            },
        };

        // Exp
        let options_exponential = cli.neutral_rate.mu_exp.map(|rate| {
            let m = Cli::normalise_mutation_rate(rate / b0, process_type);
            let u = (mu0 / (b0 * max_cells as f32)) as f64;
            let distributions = Distributions::new(cli.p_asymmetric, u, verbosity);
            let process_options = ProcessOptions {
                distributions,
                path: cli.path,
                cells2subsample: cli.subsample,
                // we assume no background mutation for the exponential growing phase
                neutral_poisson: NeutralMutationPoisson::new(m, m).expect("wrong lambda"),
            };
            SimulationOptions {
                process_options,
                gillespie_options: Options {
                    max_iter_time: IterTime {
                        iter: usize::MAX,
                        time: f32::INFINITY,
                    },
                    max_cells: max_cells - 1,
                    init_iter: 0,
                    verbosity,
                },
            }
        });

        AppOptions {
            parallel,
            fitness,
            runs,
            b0,
            seed: cli.seed,
            snapshots,
            options_moran,
            options_exponential,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verify_cli() {
        use clap::CommandFactory;
        Cli::command().debug_assert()
    }
}
