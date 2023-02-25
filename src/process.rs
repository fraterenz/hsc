use std::rc::Rc;

use rand_distr::{Bernoulli, Distribution, Poisson};
use rand_pcg::Pcg64Mcg;

use crate::sfs::{assign_cell, CloneId, NeutralMutation, StemCell, SubClone};

/// Distribution probabilities used to model acquisition of neutral and
/// conferring-fitness mutations upon cell division.
pub struct Distributions {
    poisson: Poisson<f32>,
    bern: Bernoulli,
}

impl Distributions {
    /// # Arguments
    ///
    /// * `lambda` - The parameter of the Poisson distribution modelling the
    /// acquisition of neutral mutations (the average number of neutral
    /// mutations acquired by cells upon cell division).
    /// * `p` - The parameter of a Bernouilli trial modelling the probability
    /// of acquiring a mutation conferrin a positive advantage to the cell.
    ///
    /// # Panics
    /// When p<0 or p>1 and when lambda < 0 or nan.
    pub fn new(lambda: f32, p: f64) -> Self {
        Self {
            poisson: Poisson::new(lambda).expect("Invalid lambda found: lambda <= 0 or nan"),
            bern: Bernoulli::new(p).expect("Invalid p: p<0 or p>1"),
        }
    }

    pub fn nb_neutral_mutations(&self, rng: &mut Pcg64Mcg) -> u8 {
        //! The number of neutral mutations acquired upon cell division.
        self.poisson.sample(rng) as u8
    }

    pub fn acquire_p_mutation(&self, rng: &mut Pcg64Mcg) -> bool {
        self.bern.sample(rng)
    }
}

pub fn asymmetric_division(cell: &mut StemCell, distr: &Distributions, rng: &mut Pcg64Mcg) {
    //! An asymmetric division generates a differentiated cell and a stem cell.
    //!
    //! Upon asymmetric division, a cell gets a Poisson number of neutral
    //! mutations (with rate `lambda`, see [`Distributions::new`]) and one or
    //! zero mutations conferring a fitness advantage (with probability `p`,
    //! see [`Distributions::new`]).
    for _ in 0..distr.nb_neutral_mutations(rng) {
        cell.mutations.push(NeutralMutation::new_v4());
    }

    if distr.acquire_p_mutation(rng) {
        let new_clone = SubClone {
            id: CloneId::new_v4(),
            s: 2.9,
        };
        assign_cell(cell, &Rc::new(new_clone));
    }
}

pub fn symmetric_division(cell: &mut StemCell, distr: &Distributions, rng: &mut Pcg64Mcg) {
    //! A symmetric division generates two stem cells.
    //!
    //! Since we are simulating a Moran process, another stem cell must be
    //! removed from the system each time a symmetric division occurs.
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck::{Arbitrary, Gen};
    use quickcheck_macros::quickcheck;
    use rand::SeedableRng;
    use std::num::NonZeroU8;

    #[should_panic]
    #[test]
    fn new_distribution_wrong_lambda_test() {
        Distributions::new(f32::NAN, 0.9);
    }

    #[should_panic]
    #[test]
    fn new_distribution_wrong_lambda_neg_test() {
        Distributions::new(-1.0, 0.9);
    }

    #[should_panic]
    #[test]
    fn new_distribution_wrong_p_test() {
        Distributions::new(1., f64::NAN);
    }

    #[should_panic]
    #[test]
    fn new_distribution_wrong_p_inf_test() {
        Distributions::new(1., f64::INFINITY);
    }

    #[should_panic]
    #[test]
    fn new_distribution_wrong_p_neg_test() {
        Distributions::new(1.0, -0.9);
    }
    #[derive(Clone, Debug)]
    struct LambdaFromNonZeroU8(f32);

    impl Arbitrary for LambdaFromNonZeroU8 {
        fn arbitrary(g: &mut Gen) -> LambdaFromNonZeroU8 {
            let lambda: NonZeroU8 = NonZeroU8::arbitrary(g);
            LambdaFromNonZeroU8(lambda.get() as f32)
        }
    }

    #[quickcheck]
    fn asymmetric_division_biased_one_test(lambda: LambdaFromNonZeroU8, seed: u64) -> bool {
        let mut rng = Pcg64Mcg::seed_from_u64(seed);
        let distr = Distributions::new(lambda.0, 1f64);

        let neutral_clone = SubClone {
            id: CloneId::new_v4(),
            s: 0.,
        };
        let mut cell = StemCell {
            id: 0,
            subclone: Rc::new(neutral_clone),
            mutations: Vec::new(),
        };
        asymmetric_division(&mut cell, &distr, &mut rng);
        cell.subclone.id != neutral_clone.id
    }

    #[quickcheck]
    fn asymmetric_division_biased_zero_test(lambda: LambdaFromNonZeroU8, seed: u64) -> bool {
        let mut rng = Pcg64Mcg::seed_from_u64(seed);
        let distr = Distributions::new(lambda.0, 0f64);

        let neutral_clone = SubClone {
            id: CloneId::new_v4(),
            s: 0.,
        };
        let mut cell = StemCell {
            id: 0,
            subclone: Rc::new(neutral_clone),
            mutations: Vec::new(),
        };
        asymmetric_division(&mut cell, &distr, &mut rng);
        cell.subclone.id == neutral_clone.id
    }
}
