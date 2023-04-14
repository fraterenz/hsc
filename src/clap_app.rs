use std::path::PathBuf;

use clap::{ArgAction, Parser};
use hsc::process::CellDivisionProbabilities;
use sosa::{NbIndividuals, Options};

use crate::SimulationOptions;

pub enum Parallel {
    False,
    True,
    Debug,
}

#[derive(Debug, Parser)] // requires `derive` feature
#[command(name = "hsc")]
#[command(version, about = "TODO")]
pub struct Cli {
    #[arg(long, default_value_t = 100)]
    niters: usize,
    #[arg(long, default_value_t = 100)]
    ncells: NbIndividuals,
    #[arg(long, default_value_t = 1)]
    runs: usize,
    /// division rate for the wild-type
    #[arg(long, default_value_t = 1.)]
    b0: f32,
    /// avg fit mutations arising in 1 year
    #[arg(long, default_value_t = 2.)]
    mu0: f32,
    /// avg number of neutral mutations per each proliferative event
    #[arg(long, default_value_t = 1.)]
    lambda_poisson: f32,
    #[arg(long, default_value_t = 0.15)]
    /// proliferative advantage conferred by fit mutations
    s: f32,
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
    #[arg(short, long, action = ArgAction::SetTrue, default_value_t = false, conflicts_with = "debug")]
    sequential: bool,
}

impl Cli {
    pub fn build() -> SimulationOptions {
        let cli = Cli::parse();

        let (parallel, runs) = if cli.debug {
            (Parallel::Debug, 1)
        } else if cli.sequential {
            (Parallel::False, cli.runs)
        } else {
            (Parallel::True, cli.runs)
        };

        let (max_cells, max_iter, b0, mu0, verbosity) = if cli.debug {
            (11, 20usize, 1., 4., u8::MAX)
        } else {
            (cli.ncells + 1, cli.niters, cli.b0, cli.mu0, cli.verbosity)
        };

        // rate of fit variant per cell division
        let u = if (cli.p_asymmetric - 0.).abs() <= f64::EPSILON {
            // in the symmetric case we need to divide by 2 because 1 division
            // corresponds to two cells, hence two draws (i.e. two mutational
            // events) per division
            (mu0 / (2. * b0 * max_cells as f32)) as f64
        } else {
            (mu0 / (b0 * max_cells as f32)) as f64
        };
        let probabilities = CellDivisionProbabilities {
            p_asymmetric: cli.p_asymmetric,
            lambda_poisson: cli.lambda_poisson,
            p: u,
        };

        let options = Options {
            max_iter,
            max_cells,
            init_iter: 0,
            verbosity,
        };

        SimulationOptions {
            parallel,
            s: cli.s,
            runs,
            b0,
            seed: cli.seed,
            probabilities,
            path: cli.path,
            options,
        }
    }
}

#[test]
fn verify_cli() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
