use std::{
    collections::{HashMap, HashSet},
    num::{NonZeroU64, NonZeroU8},
};

use rand::Rng;

use crate::process::NeutralMutationPoisson;

/// We assume that each division generates a different set of neutral
/// mutations (inifinte-site assumption?) that defines a genotype.
pub struct Genotype {}

/// The neutral mutations are not implemented individually, but a set of
/// mutations [`GenotypeId`] is instead assigned to each cell upon division.
///
/// At the end of the simulation, we simulate the individual mutations from the
/// sets assigned to each cell by generating a Poisson number for each
/// [`GenotypeId`] to recreate the neutral SFS.
///
/// A [`GenotypeId`] of `0` indicates no mutations.
type GenotypeId = usize;

/// The number of neutral mutations that are produced by a division event.
/// We assume that a maximal number of 255 neutral mutations can be generated
/// upon one proliferative event.
pub type NbPoissonNeutralMutations = NonZeroU8;

impl Genotype {
    pub fn mutate(stem_cell: &mut StemCell, proliferative_events: usize) {
        //! Mutate a cell by assigning it a new neutral mutation.
        //!
        //! Since we store the divisions performed by each cell, instead of the
        //! neutral mutations (to avoid generating a random number at every
        //! division?), we mutate a cell by storing the iterations at which the
        //! cell proliferates.
        //! The mapping between genotypes (i.e. proliferation id) and neutral
        //! mutations is performed by [`Genotype::sfs_neutral`].
        stem_cell.mutation_set.insert(proliferative_events);
    }

    pub fn sfs_neutral(
        cells: &[StemCell],
        poisson_dist: &NeutralMutationPoisson,
        verbosity: u8,
        rng: &mut impl Rng,
    ) -> HashMap<u64, u64> {
        //! The neutral SFS's keys are the number of mutations present in j
        //! cells (X_j) and values are the number of j cells.
        //!
        //! To plot the neutral SFS, plot on the y-axis the keys (X_j) and on
        //! the x-axis the values (j cells).
        //!
        //! Calling `sfs_neutral` transforms the set of genotypes present in
        //! `cells` [`StemCell`] into neutral unique mutations.
        //!
        //! In a Moran process with neutral selection, the expected number of
        //! variants present in j cells E(X_j) is equal to Nu(1/j) (in a sample
        //! of n cells?), where N is the total population size and u is the
        //! neutral mutation rate per unit of time, j are the cells.
        //!
        //! # Example
        //! ```
        //! use hsc::neutral::{Genotype, StemCell};
        //! use hsc::process::NeutralMutationPoisson;
        //! use rand_distr::Poisson;
        //! use rand_chacha::ChaChaRng;
        //! use rand_chacha::rand_core::SeedableRng;
        //! use std::num::NonZeroU64;
        //!
        //! let cells = [
        //!     StemCell::with_set_of_mutations(vec![1, 2]),
        //!     StemCell::with_set_of_mutations(vec![1, 10, 11])
        //! ];
        //! let verbosity = 0;
        //! // neutral mutations per cell-division
        //! let rate_neutral = 1.;
        //! let poisson = NeutralMutationPoisson(Poisson::new(rate_neutral).unwrap());
        //! let mut rng = ChaChaRng::seed_from_u64(26);
        //!
        //! let mut counts = Genotype::sfs_neutral(&cells, &poisson, verbosity, &mut rng);
        //! // removes duplicates since the vec is sorted
        //! counts.dedup();
        //! // testing requires to create a vec of unique values, since the
        //! // number of mutations is a Poisson random point process.
        //! assert_eq!(
        //!     counts,
        //!     [
        //!         NonZeroU64::new(1).unwrap(),
        //!         NonZeroU64::new(2).unwrap(),
        //!     ]
        //! );
        //! ```
        if verbosity > 0 {
            println!("Computing the sfs neutral for {} cells", cells.len());
        }
        assert!(!cells.is_empty(), "found empty cells");
        let mut total_burden = 0usize;
        let mut counts = HashMap::with_capacity(cells.len());
        // key is the genotype id and value is the number of mutations
        // it's a "database" of mutations
        let mut genotypes_poisson = HashMap::new();
        for cell in cells {
            // retrieve divisions assigned to this cell during the simulation
            for division_id in &cell.mutation_set {
                // construct the database iteratively
                if genotypes_poisson.get(&division_id).is_none() {
                    let poisson_nb = poisson_dist.nb_neutral_mutations(rng);
                    genotypes_poisson.insert(division_id, poisson_nb);
                    total_burden += poisson_nb.get() as usize;
                }
                counts
                    .entry(division_id)
                    .and_modify(|counter| *counter += 1)
                    .or_insert(1);
            }
        }
        let mut sfs = HashMap::with_capacity(total_burden);
        for (id, count) in counts.into_iter() {
            for _ in 0..genotypes_poisson[&id].get() {
                sfs.entry(count)
                    .and_modify(|counter| *counter += 1)
                    .or_insert(1);
            }
        }
        sfs
    }
}

#[derive(Debug, Clone)]
/// Hematopoietic stem and progenitor cells (HSPCs) are a rare population of
/// precursor cells that possess the capacity for self-renewal and multilineage
/// differentiation.
///
/// They carry a set of neutral mutations and are assigned to [`crate::sfs::SubClone`].
pub struct StemCell {
    /// Mutations ids which are just the id of the divisions performed by this
    /// cell.
    mutation_set: HashSet<GenotypeId>,
}

impl Default for StemCell {
    fn default() -> Self {
        //! Creates a stem cell without any neutral mutations.
        Self::new()
    }
}

impl StemCell {
    pub fn new() -> StemCell {
        //! Construct a new cell without any neutral mutations.
        StemCell {
            mutation_set: HashSet::new(),
        }
    }

    pub fn with_set_of_mutations(mutation_set: Vec<usize>) -> StemCell {
        //! Create a stem cell with a set of neutral mutations.
        //!
        //! The mutation set is not a collection of mutations but a collection
        //! of genotypes that will be converted later on into mutations, see
        //! [`Genotype::sfs_neutral`].
        let mut cell = StemCell::new();
        let mutation_set = HashSet::from_iter(mutation_set.into_iter());
        cell.mutation_set = mutation_set;
        cell
    }

    pub fn has_mutations(&self) -> bool {
        !self.mutation_set.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck::{Arbitrary, Gen};
    use quickcheck_macros::quickcheck;
    use rand::SeedableRng;
    use rand_chacha::ChaChaRng;
    use rand_distr::Poisson;
    use std::num::NonZeroU8;

    #[quickcheck]
    fn genotype_mutate_test(divisions: usize) -> bool {
        let mut cell = StemCell::new();
        Genotype::mutate(&mut cell, divisions);
        cell.has_mutations() && cell.mutation_set.len() == 1
    }

    #[derive(Debug, Clone)]
    struct DistinctMutationsId(
        NbPoissonNeutralMutations,
        NbPoissonNeutralMutations,
        NbPoissonNeutralMutations,
    );

    impl Arbitrary for DistinctMutationsId {
        fn arbitrary(g: &mut Gen) -> DistinctMutationsId {
            let mut trials = 0;
            let mut1: NonZeroU8 = NonZeroU8::arbitrary(g);
            let mut mut2 = NonZeroU8::arbitrary(g);
            let mut mut3 = NonZeroU8::arbitrary(g);
            while mut2 <= mut1 {
                mut2 = NonZeroU8::arbitrary(g);
                trials += 1;
                if trials == u8::MAX {
                    mut2 = unsafe { NonZeroU8::new_unchecked(1) };
                    break;
                }
            }
            trials = 0;

            while mut3 <= mut1 || mut3 <= mut2 {
                mut3 = NonZeroU8::arbitrary(g);
                trials += 1;
                if trials == u8::MAX {
                    mut3 = unsafe { NonZeroU8::new_unchecked(2) };
                    break;
                }
            }

            DistinctMutationsId(mut1, mut2, mut3)
        }
    }

    #[test]
    #[should_panic]
    fn test_sfs_neutral_no_cells() {
        let cells = [];
        let verbosity = 1;
        let distributions = NeutralMutationPoisson(Poisson::new(10.).unwrap());
        let mut rng = ChaChaRng::seed_from_u64(26);
        Genotype::sfs_neutral(&cells, &distributions, verbosity, &mut rng);
    }

    #[quickcheck]
    fn test_sfs_neutral_cell_no_mutations(seed: u64, verbosity: u8) -> bool {
        let cells = [StemCell::new()];
        let distributions = NeutralMutationPoisson(Poisson::new(10.).unwrap());
        let mut rng = ChaChaRng::seed_from_u64(seed);
        Genotype::sfs_neutral(&cells, &distributions, verbosity, &mut rng).is_empty()
    }

    #[quickcheck]
    fn test_sfs_neutral_one_cell(mutations: DistinctMutationsId, seed: u64) -> bool {
        let mutations = vec![
            mutations.0.get() as usize,
            mutations.1.get() as usize,
            mutations.2.get() as usize,
        ];
        let cells = [StemCell::with_set_of_mutations(mutations)];
        let verbosity = 0;
        let distributions = NeutralMutationPoisson(Poisson::new(10.).unwrap());
        let mut rng = ChaChaRng::seed_from_u64(seed);

        let mut counts = Genotype::sfs_neutral(&cells, &distributions, verbosity, &mut rng);
        // removes duplicates since the vec is sorted
        counts.dedup();
        let expected = [unsafe { NonZeroU64::new_unchecked(1) }];
        counts == expected
    }

    #[quickcheck]
    fn test_sfs_neutral_two_cells_nothing_in_common(
        mutations: DistinctMutationsId,
        seed: u64,
    ) -> bool {
        let mutations1 = vec![mutations.0.get() as usize];
        let mutations2 = vec![mutations.1.get() as usize, mutations.2.get() as usize];
        let cells = [
            StemCell::with_set_of_mutations(mutations1),
            StemCell::with_set_of_mutations(mutations2),
        ];
        let verbosity = 0;
        let distributions = NeutralMutationPoisson(Poisson::new(10.).unwrap());
        let mut rng = ChaChaRng::seed_from_u64(seed);

        let mut counts = Genotype::sfs_neutral(&cells, &distributions, verbosity, &mut rng);
        // removes duplicates since the vec is sorted
        counts.dedup();
        let expected = [unsafe { NonZeroU64::new_unchecked(1) }];
        counts == expected
    }

    #[quickcheck]
    fn test_sfs_neutral_two_cells_some_in_common(
        mutations: DistinctMutationsId,
        seed: u64,
    ) -> bool {
        let mutations1 = vec![mutations.0.get() as usize];
        let mutations2 = vec![
            mutations.0.get() as usize,
            mutations.1.get() as usize,
            mutations.2.get() as usize,
        ];
        let cells = [
            StemCell::with_set_of_mutations(mutations1),
            StemCell::with_set_of_mutations(mutations2),
        ];
        let verbosity = 0;
        let distributions = NeutralMutationPoisson(Poisson::new(10.).unwrap());
        let mut rng = ChaChaRng::seed_from_u64(seed);

        let mut counts = Genotype::sfs_neutral(&cells, &distributions, verbosity, &mut rng);
        // removes duplicates since the vec is sorted
        counts.dedup();
        let expected = [unsafe { NonZeroU64::new_unchecked(1) }, unsafe {
            NonZeroU64::new_unchecked(2)
        }];
        // we dont know how many entries there will be, since it's random (
        // Possion point process).
        counts == expected
    }

    #[quickcheck]
    fn test_sfs_neutral_two_cells_all_in_common(mutations: DistinctMutationsId, seed: u64) -> bool {
        let mutations = vec![
            mutations.0.get() as usize,
            mutations.1.get() as usize,
            mutations.2.get() as usize,
        ];
        let cells = [
            StemCell::with_set_of_mutations(mutations.clone()),
            StemCell::with_set_of_mutations(mutations),
        ];
        let verbosity = 0;
        let distributions = NeutralMutationPoisson(Poisson::new(10.).unwrap());
        let mut rng = ChaChaRng::seed_from_u64(seed);

        let mut counts = Genotype::sfs_neutral(&cells, &distributions, verbosity, &mut rng);
        // removes duplicates since the vec is sorted
        counts.dedup();
        let expected = [unsafe { NonZeroU64::new_unchecked(2) }];
        counts == expected
    }
}
