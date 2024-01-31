use anyhow::{bail, ensure};
use clap::{ArgAction, Args, Parser};
use hsc::{
    process::{ProcessOptions, Snapshot},
    subclone::{from_mean_std_to_shape_scale, Fitness},
};
use num_traits::{Float, NumCast};
use sosa::{IterTime, NbIndividuals, Options};
use std::{collections::VecDeque, ops::RangeInclusive, path::PathBuf};

use crate::{AppOptions, SimulationOptions};

#[derive(Clone, Copy, Debug)]
pub enum ProcessType {
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
#[group(required = false, multiple = false)]
pub struct FitnessArg {
    /// Neutral scenario with all clones having the same fitness coefficient of
    /// 0
    #[arg(long, action = ArgAction::SetTrue, default_value_t = false)]
    pub neutral: bool,
    /// proliferative advantage conferred by fit mutations assuming all clones
    /// have the same advantange, units: mutation / cell
    ///
    /// s cannot be greater than 1.
    #[arg(long, default_value_t = 0.11, value_parser = fitness_in_range)]
    pub s: f32,
    /// The mean and the standard deviation of the Gamma distribution used to
    /// sample the fitness coefficients representing the proliferative
    /// advantage conferred by fit mutations.
    ///
    /// The mean and the std cannot be greater than 1 and 0.1 respectively
    #[arg(long, num_args = 2)]
    pub mean_std: Option<Vec<f32>>,
}

const MEAN_RANGE: RangeInclusive<f32> = 0.01..=2.;
const STD_RANGE: RangeInclusive<f32> = 0.001..=0.5;

fn fitness_in_range(s: &str) -> Result<f32, String> {
    let fitness: f32 = s
        .parse()
        .map_err(|_| format!("`{s}` isn't a fitness number"))?;
    if MEAN_RANGE.contains(&fitness) {
        Ok(fitness)
    } else {
        Err(format!(
            "mean not in range {}-{}",
            MEAN_RANGE.start(),
            MEAN_RANGE.end()
        ))
    }
}

#[derive(Args, Debug, Clone)]
#[group(required = false, multiple = true)]
pub struct NeutralMutationRate {
    /// Poisson rate for the exponential growing phase, leave empty to simulate
    /// just a Moran process with fixed population size
    #[arg(long)]
    pub mu_exp: Option<f32>,
    /// Background mutation rate (all neutral mutations **not** occuring in the
    /// mitotic phase) for for the constant population phase
    /// units: [mut/year]
    #[arg(long, default_value_t = 12.)]
    pub mu_background: f32,
    /// Division mutation rate (all neutral mutations occuring upon
    /// cell-division, mitotic phase) for for the constant population phase
    /// units: [mut/division]
    #[arg(long, default_value_t = 1.2)]
    pub mu_division: f32,
}

#[derive(Debug, Parser)] // requires `derive` feature
#[command(name = "hsc", version, about = "TODO")]
pub struct Cli {
    /// The years simulated will be `years + 1`, that is simulations start at
    /// year zero and the end at year `year + 1`
    #[arg(short, long, default_value_t = 9)]
    years: usize,
    /// Number of cells to simulate
    #[arg(short, long, default_value_t = 100)]
    cells: NbIndividuals,
    #[arg(short, long, default_value_t = 1)]
    runs: usize,
    /// time between wild-type stem cells divisions in years, units: year / (division * cell)
    #[arg(long, default_value_t = 1.)]
    tau: f32,
    /// avg fit mutations arising in 1 year, units: division / year. If not passed
    /// generated a random value between 1 and 14
    #[arg(long, default_value_t = 4.)]
    mu0: f32,
    /// avg number of neutral mutations per each proliferative event assuming a
    /// Poisson distribution, units: mutation / (cell * year)
    #[command(flatten)]
    neutral_rate: NeutralMutationRate,
    #[command(flatten)]
    /// If not present, use a constant fitness of 0.11
    fitness: FitnessArg,
    #[arg(long, default_value_t = 26)]
    seed: u64,
    /// Save the whole population
    #[arg(long, action = ArgAction::SetTrue, default_value_t = false)]
    save_population: bool,
    /// Save only the SFS (when you dont want too many files to be saved)
    #[arg(long, action = ArgAction::SetTrue, default_value_t = false)]
    save_sfs_only: bool,
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
    #[arg(long, requires = "subsamples", num_args = 0..)]
    /// Snapshots to take to save the simulation, requires `subsamples`.
    ///
    /// The combination of `snapshots` with `subsamples` gives four different
    /// behaviours:
    ///
    /// 1. when `snapshots.len() = 1` and `subsamples.len() = 1`: subsample once with the number of cells corresponding to `snapshots[0]`
    ///
    /// 2. when `snapshots.len() > 1` and `subsamples.len() = 1`: for every `s` in `snapshots`, subsample with the number of cells corresponding to `snapshots[0]`
    ///
    /// 3. when `snapshots.len() = 1` and `subsamples.len() > 1`: for every `c` in `subsamples`, subsample once with the number of cells corresponding to `c`
    ///
    /// 4. when `snapshots.len() > 1` and `subsamples.len() > 1`: for every pair `(s, c)` in `snapshots.zip(subsamples)`, subsample at time `s` with `c` cells
    snapshots: Option<Vec<f32>>,
    /// Number of cells to subsample before saving the measurements, leave
    /// empty when no subsample is needed.
    /// If subsampling is performed, the measurements of the whole population
    /// will also be saved.
    ///
    /// See help for `snapshots` for more details.
    #[arg(long, requires = "snapshots", num_args = 0..)]
    subsamples: Option<Vec<usize>>,
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

    pub fn normalise_mutation_rate<F: Float>(rate: F, process_type: ProcessType) -> F {
        match process_type {
            ProcessType::MoranAsymmetric => rate,
            ProcessType::MoranSymmetric => rate / NumCast::from(2).unwrap(),
        }
    }

    pub fn build() -> anyhow::Result<AppOptions> {
        let cli = Cli::parse();

        let (parallel, runs) = if cli.debug {
            (Parallel::Debug, 1)
        } else if cli.sequential {
            (Parallel::False, cli.runs)
        } else {
            (Parallel::True, cli.runs)
        };

        let years = cli.years;
        let (max_cells, years, verbosity) = if cli.debug {
            (11, 1, u8::MAX)
        } else {
            (cli.cells + 1, years, cli.verbosity)
        };

        let max_iter = 10 * max_cells as usize * 10 * years;

        let mut snapshots = match (cli.subsamples, cli.snapshots) {
            (Some(sub), Some(snap)) => {
                match (&sub[..], &snap[..]) {
                    // subsample `unique_sub` once at `unique_snap` time
                    ([unique_sub], [unique_snap]) => VecDeque::from_iter(
                        [(unique_sub, unique_snap)]
                            .into_iter()
                            .map(|(&cells2sample, &time)| Snapshot { cells2sample, time }),
                    ),
                    // subsample `unique_sub` at several `snaps` times
                    ([unique_sub], snaps) => VecDeque::from_iter(
                        snaps
                            .iter()
                            .zip(std::iter::repeat(unique_sub))
                            .map(|(&time, &cells2sample)| Snapshot { cells2sample, time }),
                    ),
                    // subsample with different cells `unique_sub` once at `unique_snap` time
                    (subs, [unique_snap]) => VecDeque::from_iter(
                        subs.iter()
                            .zip(std::iter::repeat(unique_snap))
                            .map(|(&cells2sample, &time)| Snapshot { cells2sample, time }),
                    ),
                    // subsample with many cells `subs` at specific `snaps`
                    (subs, snaps) => {
                        ensure!(
                            subs.len() == snaps.len(),
                            "the lenght of snapshots do not match the lenght of subsamples"
                        );
                        VecDeque::from_iter(
                            subs.iter()
                                .zip(snaps)
                                .map(|(&cells2sample, &time)| Snapshot { cells2sample, time }),
                        )
                    }
                }
            }
            (None, None) => VecDeque::from_iter(
                Cli::build_snapshots_from_time(10usize, years as f32 - 1.)
                    .into_iter()
                    .map(|t| Snapshot {
                        cells2sample: cli.cells as usize,
                        time: t,
                    }),
            ),
            _ => unreachable!(),
        };

        snapshots.make_contiguous();
        snapshots
            .as_mut_slices()
            .0
            .sort_by(|s1, s2| s1.time.partial_cmp(&s2.time).unwrap());

        ensure!(
            snapshots.iter().all(|s| s.time < cli.years as f32),
            "times to take `snapshots` must be smaller than total `years`"
        );
        ensure!(
            snapshots
                .iter()
                .all(|s| s.cells2sample <= cli.cells as usize),
            "the cells to subsample must be smaller or equal than the total population size"
        );

        // Moran
        let process_options = ProcessOptions {
            path: cli.path.clone(),
            snapshots: snapshots.clone(),
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
            save_sfs_only: cli.save_sfs_only,
            save_population: cli.save_population,
        };

        // Exp
        let options_exponential = cli.neutral_rate.mu_exp.map(|_| {
            assert!(cli.neutral_rate.mu_exp.is_some());
            let process_options = ProcessOptions {
                path: cli.path,
                snapshots: snapshots.clone(),
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
                save_sfs_only: cli.save_sfs_only,
                save_population: cli.save_population,
            }
        });

        let fitness = if let Some(mean_std) = cli.fitness.mean_std {
            if !MEAN_RANGE.contains(&mean_std[0]) {
                bail!(format!(
                    "fitness not in range {}-{}",
                    MEAN_RANGE.start(),
                    MEAN_RANGE.end()
                ));
            }
            if !STD_RANGE.contains(&mean_std[1]) {
                bail!(format!(
                    "std not in range {}-{}",
                    STD_RANGE.start(),
                    STD_RANGE.end()
                ));
            }
            let (shape, scale) = from_mean_std_to_shape_scale(mean_std[0], mean_std[1]);
            Fitness::GammaSampled { shape, scale }
        } else if cli.fitness.neutral {
            Fitness::Neutral
        } else {
            Fitness::Fixed { s: cli.fitness.s }
        };

        Ok(AppOptions {
            parallel,
            fitness,
            runs,
            tau: cli.tau,
            seed: cli.seed,
            snapshots,
            options_moran,
            options_exponential,
            mu0: cli.mu0,
            neutral_rate: cli.neutral_rate,
            verbosity: cli.verbosity,
        })
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
