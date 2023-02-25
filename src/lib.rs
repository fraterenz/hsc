//! Simulate the dynamics of a stem cell population undergoing proliferation
//! and differentiation according to a Moran process (fixed population size).
//!
//! At every cell division, a stem cell can undergo one of the following events:
//!
//! 1. a symmetric division (see [`process::symmetric_division`]),
//! 2. an asymmetric division (see [`process::asymmetric_division`])
//!
//! Whenever a cell divides, it acquires a Poisson number of passenger
//! mutations, which do not give any proliferative advantage.
//!
//! On top of that, cells can also acquire one proliferative advantageous
//! mutation upon divison.
//! In this case, the cell creates a new clone, which has a birth-rate of
//! `lambda_i = lambda_0 ( 1 + s_i)`, where `lambda_0` is the birth-rate of the
//! wild-type, the clone without any proliferative mutations.

/// The events to simulate for this Markov process.
pub mod process;
/// The site frequency spectrum (also known as the variant allele frequency).
pub mod sfs;
