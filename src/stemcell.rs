use std::{
    collections::{HashMap, HashSet},
    num::NonZeroU8,
};

use rand::Rng;

use crate::process::NeutralMutationPoisson;

/// The genotype of a [`StemCell`].
///
/// We assume that each division generates a different set of mutations
/// (inifinte-site assumption?) that defines a genotype.
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
        //! Mutate a cell by assigning it a new mutation.
        //!
        //! Since we store the divisions performed by each cell, instead of the
        //! mutations (to avoid generating a random number at every division?),
        //! we mutate a cell by storing the iterations at which the cell
        //! proliferates.
        //! The mapping between genotypes (i.e. proliferation id) and mutations
        //! is performed by [`Genotype::sfs`].
        stem_cell.mutation_set.insert(proliferative_events);
    }

    pub fn sfs(
        cells: &[StemCell],
        poisson_dist: &NeutralMutationPoisson,
        verbosity: u8,
        rng: &mut impl Rng,
    ) -> HashMap<u64, u64> {
        //! Compute the site frequency spectrum (SFS) of mutations found in the
        //! stem cell population `cells`.
        //! The SFS's keys are the j cells and values are the number of
        //! mutations present in j cells (X_j).
        //! To plot the SFS, plot on the x-axis the keys (j cells) and on the
        //! y-axis the values (X_j mutations in j cells).
        //!
        //! Calling `sfs` transforms the set of genotypes present in `cells`
        //! [`StemCell`] into unique mutations.
        //!
        //! In a Moran process with neutral selection, the expected number of
        //! variants present in j cells E\[X_j\] is equal to Nu(1/j) (in a
        //! sample of n cells?), where N is the total population size and u is
        //! the neutral mutation rate per unit of time, j are the cells.
        //!
        //! # Example
        //! ```
        //! use hsc::stemcell::{Genotype, StemCell};
        //! use hsc::process::NeutralMutationPoisson;
        //! use rand_distr::Poisson;
        //! use rand_chacha::ChaChaRng;
        //! use rand_chacha::rand_core::SeedableRng;
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
        //! let counts = Genotype::sfs(&cells, &poisson, verbosity, &mut rng);
        //! let expected_jcells = [1u64, 2u64];
        //! let mut jcells = counts.into_keys().collect::<Vec<u64>>();
        //! jcells.sort_unstable();
        //! assert_eq!(jcells, expected_jcells);
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
/// They carry a set of neutral mutations and are assigned to [`crate::subclone::SubClone`].
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
        //! [`Genotype::sfs`].
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
        Genotype::sfs(&cells, &distributions, verbosity, &mut rng);
    }

    #[quickcheck]
    fn test_sfs_neutral_cell_no_mutations(seed: u64, verbosity: u8) -> bool {
        let cells = [StemCell::new()];
        let distributions = NeutralMutationPoisson(Poisson::new(10.).unwrap());
        let mut rng = ChaChaRng::seed_from_u64(seed);
        Genotype::sfs(&cells, &distributions, verbosity, &mut rng).is_empty()
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

        let counts = Genotype::sfs(&cells, &distributions, verbosity, &mut rng);
        let expected_jcells = [1u64];
        counts.into_keys().collect::<Vec<u64>>() == expected_jcells
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

        let counts = Genotype::sfs(&cells, &distributions, verbosity, &mut rng);
        let expected_jcells = [1u64];
        counts.into_keys().collect::<Vec<u64>>() == expected_jcells
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

        let counts = Genotype::sfs(&cells, &distributions, verbosity, &mut rng);
        let expected_jcells = [1u64, 2u64];
        let mut jcells = counts.into_keys().collect::<Vec<u64>>();
        jcells.sort_unstable();
        jcells == expected_jcells
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

        let counts = Genotype::sfs(&cells, &distributions, verbosity, &mut rng);
        let expected_jcells = [2u64];
        counts.into_keys().collect::<Vec<u64>>() == expected_jcells
    }
}
