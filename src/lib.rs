//! Simulate the dynamics of a stem cell population undergoing proliferation
//! and differentiation according to a Moran process (fixed population size).
//!
//! Whenever a cell divides, it acquires a Poisson number of passenger
//! mutations, which do not give any proliferative advantage.
//!
//! On top of that, cells can also acquire one proliferative advantageous
//! mutation upon division.
//! In this case, the cell creates a new clone, which has a birth-rate of
//! `lambda_i = lambda_0 ( 1 + s_i)`, where `lambda_0` is the birth-rate of the
//! wild-type, the clone without any proliferative mutations.

use std::{
    fs,
    io::{BufWriter, Write},
    path::Path,
};

use anyhow::Context;

/// The neutral mutations representing the genotype of the stem cells.
pub mod genotype;
/// The events to simulate for this Markov process.
pub mod process;
/// The proliferation of cells with the simulation of neutral and fit mutations.
pub mod proliferation;
/// The agents whose state defines the system simulated by the process.
pub mod stemcell;
/// The classes defining the proliferative advantage.
pub mod subclone;

/// Maximal number of fit clones that can arise during the simulation.
///
/// If the parameters of the simulation provided by the user (e.g. fit mutation
/// rate) are too extreme, the program will exit with error. To avoid this,
/// increase here the number of clones.
pub const MAX_SUBCLONES: usize = 1200;

/// The time at birth measured in years used for background mutations in the
/// exponential growing phase.
pub const TIME_AT_BIRTH: f32 = 9. / 12.;

#[derive(Clone, Debug)]
pub struct ProbsPerYear {
    /// Arrival rate of neutral background mutations per year
    pub mu_background: f32,
    /// Arrival rate of neutral "divisional" mutations per year
    pub mu_division: f32,
    /// Arrival rate of fit mutants per year
    pub mu: f32,
}

/// Probabilities used in the simulations.
#[derive(Clone, Debug)]
pub enum Probs {
    Asymmetric {
        /// Arrival rate of fit mutants per cell per division
        u: f32,
        /// Probabilities per year
        probs_per_year: ProbsPerYear,
        /// Probability of asymmetric division per cell per division
        asymmetric: f32,
    },
    Symmetric {
        /// Arrival rate of fit mutants per cell per division
        u: f32,
        /// Probabilities per year
        probs_per_year: ProbsPerYear,
    },
}

impl Probs {
    pub fn new(
        mu_background: f32,
        mu_division: f32,
        mu: f32,
        asymmetric: f32,
        cells: u64,
        verbosity: u8,
    ) -> Probs {
        //! ## Panics
        //! Panics when `mu` is greater than `cells` or when asymmetric is not
        //! within interval of 0 and 1.
        assert!(mu <= cells as f32);
        let probs_per_year = ProbsPerYear {
            mu_background,
            mu_division,
            mu,
        };
        let u = mu / (cells as f32);
        let probs = if (asymmetric - 0.).abs() > f32::EPSILON {
            Probs::Asymmetric {
                u,
                probs_per_year,
                asymmetric,
            }
        } else {
            Probs::Symmetric { u, probs_per_year }
        };
        if verbosity > 0 {
            println!("probs {probs:#?}");
        }
        assert!((0f32..1.).contains(&u), "Invalid u: u>=0 and u<1");
        assert!(
            (0f32..=1.).contains(&asymmetric),
            "Invalid asymmetric: asymmetric>=0 and asymmetric<=1"
        );
        probs
    }

    pub fn is_asymmetric(&self) -> bool {
        match self {
            Probs::Symmetric { .. } => false,
            Probs::Asymmetric { .. } => true,
        }
    }
}

pub fn write2file<T: std::fmt::Display>(
    data: &[T],
    path: &Path,
    header: Option<&str>,
    endline: bool,
) -> anyhow::Result<()> {
    //! Write vector of float into new file with a precision of 6 decimals.
    //! Write NAN if the slice to write to file is empty.
    fs::create_dir_all(path.parent().unwrap()).expect("Cannot create dir");
    let f = fs::OpenOptions::new()
        .read(true)
        .append(true)
        .create(true)
        .open(path)
        .with_context(|| "Cannot open stream")?;

    let mut buffer = BufWriter::new(f);

    if !data.is_empty() {
        if let Some(h) = header {
            writeln!(buffer, "{h}")?;
        }

        for ele in data.iter() {
            write!(buffer, "{ele:.6},")?;
        }

        if endline {
            writeln!(buffer)?;
        }
    } else {
        write!(buffer, "{},", f32::NAN)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck::{Arbitrary, Gen};
    use std::num::NonZeroU8;

    #[derive(Clone, Debug)]
    pub struct LambdaFromNonZeroU8(pub f32);

    impl Arbitrary for LambdaFromNonZeroU8 {
        fn arbitrary(g: &mut Gen) -> LambdaFromNonZeroU8 {
            let lambda: NonZeroU8 = NonZeroU8::arbitrary(g);
            LambdaFromNonZeroU8(lambda.get() as f32)
        }
    }

    #[test]
    #[should_panic]
    fn panic_asymmetric_neg_cells_test() {
        Probs::new(1.1, 1.1, 0.1, -0.1, 10, 0);
    }

    #[test]
    #[should_panic]
    fn panic_asymmetric_inf_cells_test() {
        Probs::new(1.1, 1.1, 0.1, f32::INFINITY, 10, 0);
    }

    #[test]
    #[should_panic]
    fn panic_asymmetric_nan_cells_test() {
        Probs::new(1.1, 1.1, 0.1, f32::NAN, 10, 0);
    }

    #[test]
    #[should_panic]
    fn panic_mu_gr_cells_test() {
        Probs::new(1.1, 1.1, 12., 0., 10, 0);
    }

    #[test]
    #[should_panic]
    fn panic_mu_neg_cells_test() {
        Probs::new(1.1, 1.1, -0.1, 0., 10, 0);
    }

    #[test]
    #[should_panic]
    fn panic_mu_inf_cells_test() {
        Probs::new(1.1, 1.1, f32::INFINITY, 0., 10, 0);
    }

    #[test]
    #[should_panic]
    fn panic_mu_nan_cells_test() {
        Probs::new(1.1, 1.1, f32::NAN, 0., 10, 0);
    }
}
