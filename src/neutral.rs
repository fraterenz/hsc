use std::{
    collections::{HashMap, HashSet},
    num::NonZeroU8,
};

use rand::Rng;

use crate::process::Distributions;

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
pub type GenotypeId = usize;

/// The number of neutral mutations that are produced by a division event.
/// We assume that a maximal number of 255 neutral mutations can be generated
/// upon one proliferative event.
pub type NbPoissonNeutralMutations = NonZeroU8;

impl Genotype {
    pub fn mutate(stem_cell: &mut StemCell, proliferative_events: usize) {
        //! Mutate a cell by assigning it new neutral mutation based on the
        //! current iteration `proliferative_events`.
        stem_cell.mutation_set.insert(proliferative_events);
    }

    pub fn sfs_neutral(
        cells: &[StemCell],
        distributions: &Distributions,
        verbosity: u8,
        rng: &mut impl Rng,
    ) -> HashMap<NbPoissonNeutralMutations, usize> {
        //! The key is the number of mutations, the value if the number of
        //! cells with that number of mutation.
        // TODO: maybe reduce the capacity
        if verbosity > 0 {
            println!("Computing the sfs neutral for {} cells", cells.len());
        }
        let mut sfs = HashMap::with_capacity(cells.len());
        // key is the genotype id and value is the number of mutations
        let mut genotypes_poisson = HashMap::new();
        for cell in cells {
            for mutation in &cell.mutation_set {
                // mutation 0 means no neutral mutations
                if mutation == &0 {
                    continue;
                }
                if genotypes_poisson.get(&mutation).is_none() {
                    let poisson_nb = distributions.nb_neutral_mutations(rng);
                    genotypes_poisson.insert(mutation, poisson_nb);
                }
                sfs.entry(genotypes_poisson[&mutation])
                    .and_modify(|counter| *counter += 1)
                    .or_insert(1);
            }
        }
        sfs
    }
}

#[derive(Debug, Clone)]
pub struct StemCell {
    /// Mutations ids
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

    pub fn has_mutations(&self) -> bool {
        !self.mutation_set.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck::{Arbitrary, Gen};
    use quickcheck_macros::quickcheck;
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
            let mut1: NonZeroU8 = NonZeroU8::arbitrary(g);
            let mut mut2 = NonZeroU8::arbitrary(g);
            let mut mut3 = NonZeroU8::arbitrary(g);
            while mut2 <= mut1 {
                mut2 = NonZeroU8::arbitrary(g);
            }

            while mut3 <= mut1 || mut3 <= mut2 {
                mut3 = NonZeroU8::arbitrary(g);
            }

            DistinctMutationsId(mut1, mut2, mut3)
        }
    }
}
