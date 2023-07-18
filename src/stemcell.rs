use rand::Rng;
use std::{
    collections::{hash_map::Entry, HashMap, HashSet},
    num::{NonZeroU64, NonZeroU8, NonZeroUsize},
};

use crate::process::NeutralMutationPoisson;

/// The mutations are not implemented individually, but a set of mutations
/// [`GenotypeId`] is instead assigned to each cell upon division.
///
/// At the end of the simulation, we simulate the individual mutations from the
/// sets assigned to each cell by generating a Poisson number for each
/// [`GenotypeId`] to recreate the SFS.
type GenotypeId = usize;

/// The number of mutations that are produced by a division event.
/// We assume that a maximal number of 255 neutral mutations can be generated
/// upon one proliferative event.
pub type NbPoissonMutations = NonZeroU8;

/// Site frequency spectrum ([SFS]()) is the distribution of the allele
/// frequencies of a given set of loci (often SNPs) in a population or sample.
/// The SFS's here is implemeted as a mapping, whose keys are the j cells and
/// values are the number of mutations present in j cells (X_j).
///
/// To plot the SFS, plot on the x-axis the keys (j cells) and on the
/// y-axis the values (X_j mutations in j cells), transform axes in log scale.
///
/// In a Moran process with neutral selection, the expected number of
/// variants present in j cells E\[X_j\] is equal to Nu(1/j) (in a
/// sample of n cells?), where N is the total population size and u is
/// the neutral mutation rate per unit of time, j are the cells.
pub struct Sfs(pub HashMap<u64, u64>);

impl Sfs {
    pub fn from_cells(
        cells: &[StemCell],
        poisson_dist: &NeutralMutationPoisson,
        verbosity: u8,
        rng: &mut impl Rng,
    ) -> HashMap<u64, u64> {
        //! Compute the site frequency spectrum (SFS) of mutations found in the
        //! stem cell population `cells`.
        //!
        //! Calling `sfs` transforms the set of proliferation events present in
        //! `cells` into unique mutations, see [`StemCell::proliferation_events_id`].
        //!
        //! # Example
        //! ```
        //! use hsc::stemcell::{Sfs, StemCell};
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
        //! let counts = Sfs::from_cells(&cells, &poisson, verbosity, &mut rng);
        //! let expected_jcells = [1u64, 2u64];
        //! let mut jcells = counts.into_keys().collect::<Vec<u64>>();
        //! jcells.sort_unstable();
        //! assert_eq!(jcells, expected_jcells);
        //! ```
        if verbosity > 0 {
            println!("Computing the sfs for {} cells", cells.len());
        }
        let stats = StatisticsMutations::from_cells(cells, poisson_dist, rng);
        if verbosity > 1 {
            println!("stats for SFS: {:#?}", stats);
        }
        let sfs = Sfs::from_stats(stats);
        if verbosity > 1 {
            println!("sfs: {:#?}", sfs);
        }
        sfs
    }

    fn from_stats(stats: StatisticsMutations) -> HashMap<u64, u64> {
        let mut sfs = HashMap::with_capacity(stats.total_burden.get());
        for stats_mut_id in stats.counts.into_values() {
            for _ in 0..stats_mut_id.poisson_mut_number.get() {
                sfs.entry(stats_mut_id.cell_count.get())
                    .and_modify(|counter| *counter += 1)
                    .or_insert(1);
            }
        }
        sfs
    }
}

/// Mapping between [`GenotypeId`] and the number of mutations and the number
/// of cells.
#[derive(Debug)]
struct StatisticsMutations {
    counts: HashMap<GenotypeId, CountsMutations>,
    total_burden: NonZeroUsize,
}

#[derive(Debug)]
struct CountsMutations {
    // the number of cells carrying a mutation id
    cell_count: NonZeroU64,
    // the number of mutations corresponding to a mutation id
    poisson_mut_number: NonZeroU64,
}

impl StatisticsMutations {
    fn from_cells(
        cells: &[StemCell],
        poisson_dist: &NeutralMutationPoisson,
        rng: &mut impl Rng,
    ) -> Self {
        assert!(!cells.is_empty(), "found empty cells");
        let mut total_burden = 0usize;
        let mut counts: HashMap<GenotypeId, CountsMutations> = HashMap::with_capacity(cells.len());
        for cell in cells {
            // retrieve divisions assigned to this cell during the simulation
            for division_id in &cell.proliferation_events_id {
                let mutation_id = counts.entry(*division_id);
                let mut_number = match &mutation_id {
                    // if it's the first time we encounter this mutation id
                    // in the population, we assign to this mut id some poisson
                    // mutations
                    Entry::Vacant(_) => {
                        let poisson_nb = NonZeroU64::from(poisson_dist.nb_neutral_mutations(rng));
                        total_burden += poisson_nb.get() as usize;
                        poisson_nb
                    }
                    // else we retrieve the poisson mutations already assigned
                    Entry::Occupied(entry) => entry.get().poisson_mut_number,
                };
                // we update the counter of cells by adding this cell to the
                // mutation id. We also update the counter of mutations
                mutation_id
                    .and_modify(|stat| {
                        stat.cell_count =
                            unsafe { NonZeroU64::new_unchecked(1 + stat.cell_count.get()) };
                        stat.poisson_mut_number = mut_number;
                    })
                    .or_insert(CountsMutations {
                        cell_count: unsafe { NonZeroU64::new_unchecked(1) },
                        poisson_mut_number: mut_number,
                    });
            }
        }
        StatisticsMutations {
            counts,
            total_burden: NonZeroUsize::new(total_burden).expect("Found tot burden of 0"),
        }
    }
}

#[derive(Debug, Clone)]
/// Hematopoietic stem and progenitor cells (HSPCs) are a rare population of
/// precursor cells that possess the capacity for self-renewal and multilineage
/// differentiation.
///
/// They carry a set of neutral mutations and are assigned to [`crate::subclone::SubClone`].
pub struct StemCell {
    /// A collection of ids which identify the iterations upon which the cell
    /// prolfierates during the simulation.
    pub proliferation_events_id: HashSet<GenotypeId>,
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
            proliferation_events_id: HashSet::new(),
        }
    }

    pub fn with_set_of_mutations(mutation_set: Vec<usize>) -> StemCell {
        //! Create a stem cell with a set of neutral mutations.
        //!
        //! The mutation set is not a collection of mutations but a collection
        //! of genotypes that will be converted later on into mutations, see
        //! [`Sfs::from_cells`].
        let mut cell = StemCell::new();
        let mutation_set = HashSet::from_iter(mutation_set.into_iter());
        cell.proliferation_events_id = mutation_set;
        cell
    }

    pub fn has_mutations(&self) -> bool {
        !self.proliferation_events_id.is_empty()
    }

    pub fn record_proliferation_event(&mut self, event_id: usize) {
        //! Record the iteration when the cell has undergone cell-division
        //! (proliferation).
        //!
        //! Since we store the divisions performed by each cell, instead of the
        //! mutations (to avoid generating a random number at every division?),
        //! we mutate a cell by storing the iterations at which the cell
        //! proliferates.
        //! The mapping between genotypes (i.e. proliferation id) and mutations
        //! is performed by [`Sfs::from_cells`].
        self.proliferation_events_id.insert(event_id);
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
        cell.record_proliferation_event(divisions);
        cell.has_mutations() && cell.proliferation_events_id.len() == 1
    }

    #[derive(Debug, Clone)]
    struct DistinctMutationsId(NbPoissonMutations, NbPoissonMutations, NbPoissonMutations);

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
        Sfs::from_cells(&cells, &distributions, verbosity, &mut rng);
    }

    #[test]
    #[should_panic]
    fn test_sfs_neutral_cell_no_mutations() {
        let cells = [StemCell::new()];
        let distributions = NeutralMutationPoisson(Poisson::new(10.).unwrap());
        let mut rng = ChaChaRng::seed_from_u64(26);
        Sfs::from_cells(&cells, &distributions, 4, &mut rng).is_empty();
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

        let counts = Sfs::from_cells(&cells, &distributions, verbosity, &mut rng);
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

        let counts = Sfs::from_cells(&cells, &distributions, verbosity, &mut rng);
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

        let counts = Sfs::from_cells(&cells, &distributions, verbosity, &mut rng);
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

        let counts = Sfs::from_cells(&cells, &distributions, verbosity, &mut rng);
        let expected_jcells = [2u64];
        counts.into_keys().collect::<Vec<u64>>() == expected_jcells
    }

    #[quickcheck]
    fn test_sfs_from_stats(idx: DistinctMutationsId) -> bool {
        let random_nbs: [NonZeroU64; 6] = std::array::from_fn(|i| {
            NonZeroU64::new(idx.0.get() as u64 + i as u64).expect("overflow?")
        });
        let counts: HashMap<GenotypeId, CountsMutations> = HashMap::from([
            (
                random_nbs[0].get() as usize,
                CountsMutations {
                    cell_count: random_nbs[1],
                    poisson_mut_number: random_nbs[2],
                },
            ),
            (
                random_nbs[3].get() as usize,
                CountsMutations {
                    cell_count: random_nbs[4],
                    poisson_mut_number: random_nbs[5],
                },
            ),
        ]);
        let stats = StatisticsMutations {
            counts,
            total_burden: unsafe {
                NonZeroUsize::new_unchecked(
                    random_nbs[1].get() as usize * random_nbs[2].get() as usize
                        + random_nbs[4].get() as usize * random_nbs[5].get() as usize,
                )
            },
        };
        let sfs = Sfs::from_stats(stats);
        let mut keys = sfs.clone().into_keys().collect::<Vec<u64>>();
        let mut values = sfs.into_values().collect::<Vec<u64>>();
        keys.sort_unstable();
        values.sort_unstable();

        let mut expected_keys = vec![random_nbs[1].get(), random_nbs[4].get()];
        let mut expected_values = vec![random_nbs[2].get(), random_nbs[5].get()];
        expected_keys.sort_unstable();
        expected_values.sort_unstable();

        keys == expected_keys && values == expected_values
    }
}
