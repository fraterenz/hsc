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
            writeln!(buffer, "{}", h)?;
        }

        for ele in data.iter() {
            write!(buffer, "{:.6},", ele)?;
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
}
