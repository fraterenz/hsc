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
    /// Number of years to simulate
    #[arg(short, long, default_value_t = 1)]
    years: usize,
    /// Number of cells to simulate
    #[arg(short, long, default_value_t = 100)]
    cells: NbIndividuals,
    #[arg(short, long, default_value_t = 1)]
    runs: usize,
    /// division rate for the wild-type in 1 year
    #[arg(long, default_value_t = 1.)]
    b0: f32,
    /// avg fit mutations arising in 1 year
    #[arg(long, default_value_t = 2.)]
    mu0: f32,
    /// avg number of neutral mutations per each proliferative event assuming a
    /// Poisson distribution
    #[arg(long, default_value_t = 1.)]
    neutral_rate: f32,
    #[arg(short, default_value_t = 0.15)]
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
    #[arg(long, action = ArgAction::SetTrue, default_value_t = false, conflicts_with = "debug")]
    sequential: bool,
}

impl Cli {
    fn build_snapshots_from_time(time: f32) -> Vec<f32> {
        let n = 10;
        let dx = time / ((n - 1) as f32);
        let mut x = vec![0.1; n];
        for i in 1..n - 1 {
            x[i] = x[i - 1] + dx;
        }

        x.shrink_to_fit();
        x[n - 1] = time;
        x
    }

    pub fn build() -> SimulationOptions {
        let cli = Cli::parse();

        let (parallel, runs) = if cli.debug {
            (Parallel::Debug, 1)
        } else if cli.sequential {
            (Parallel::False, cli.runs)
        } else {
            (Parallel::True, cli.runs)
        };

        let (max_cells, years, b0, mu0, verbosity) = if cli.debug {
            (11, 1, 1., 4., u8::MAX)
        } else {
            (cli.cells + 1, cli.years + 1, cli.b0, cli.mu0, cli.verbosity)
        };

        let max_iter = max_cells as usize * b0 as usize * years;
        let snapshots = Cli::build_snapshots_from_time(cli.years as f32);

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
            lambda_poisson: cli.neutral_rate,
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
            snapshots,
            probabilities,
            path: cli.path,
            options,
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