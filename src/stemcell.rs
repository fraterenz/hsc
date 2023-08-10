use anyhow::ensure;
use rand::Rng;
use rustc_hash::{FxHashMap, FxHashSet};
use std::{
    collections::{hash_map::Entry, HashSet},
    num::NonZeroU64,
};

use crate::process::NeutralMutationPoisson;

/// The mutations are not implemented individually, but a set of mutations
/// [`GenotypeId`] is instead assigned to each cell upon division.
///
/// At the end of the simulation, we simulate the individual mutations from the
/// sets assigned to each cell by generating a Poisson number for each
/// [`GenotypeId`] to recreate the SFS.
pub type GenotypeId = usize;

/// The number of mutations that are produced by a division event.
/// We assume that a maximal number of 255 neutral mutations can be generated
/// upon one proliferative event.
pub type NbPoissonMutations = u8;

/// Site frequency spectrum ([SFS](https://en.wikipedia.org/wiki/Allele_frequency_spectrum))
/// is the distribution of the allele frequencies of a given set of loci in a
/// population or sample.
/// The SFS's here is implemeted as a mapping, whose values are the number of
/// cells (X_j) and whose keys are the number of mutations (Y_j) jmuts present
/// in j cells.
///
/// To plot the SFS, plot on the x-axis the values (jcells) and on the y-axis
/// the keys (X_j mutations in j cells, jmuts), transform axes in log scale.
pub struct Sfs(pub FxHashMap<u64, u64>);

impl Sfs {
    pub fn from_stats(stats: &StatisticsMutations, verbosity: u8) -> anyhow::Result<Self> {
        //! Compute the site frequency spectrum (SFS) of mutations found in the
        //! stem cell population.
        //!
        //! The SFS is implemented with keys being the jmuts and the values
        //! being the jcells with with jmuts.
        // let mut sfs = HashMap::with_capacity(stats.nb_variants);
        ensure!(!stats.counts.is_empty(), "found empty stats");
        let mut sfs: FxHashMap<u64, u64> = FxHashMap::default();
        sfs.shrink_to(stats.nb_variants);
        for stats_mut_id in stats.counts.values() {
            sfs.entry(stats_mut_id.poisson_mut_number)
                .and_modify(|jcells| *jcells += stats_mut_id.cell_count.get())
                .or_insert(stats_mut_id.cell_count.get());
        }
        if verbosity > 1 {
            println!("sfs: {:#?}", sfs);
        }
        Ok(Sfs(sfs))
    }
}

#[derive(Debug, Clone)]
struct CountsMutations {
    // the number of cells carrying a mutation id
    cell_count: NonZeroU64,
    // the number of mutations corresponding to a mutation id
    poisson_mut_number: u64,
}

/// Mapping between [`GenotypeId`] and the number of mutations and the number
/// of cells.
#[derive(Debug, Default, Clone)]
pub struct StatisticsMutations {
    counts: FxHashMap<GenotypeId, CountsMutations>,
    nb_variants: usize,
}

impl StatisticsMutations {
    pub fn from_cells(
        cells: &[StemCell],
        poisson_dist: &NeutralMutationPoisson,
        rng: &mut impl Rng,
        verbosity: u8,
    ) -> anyhow::Result<Self> {
        //! Create the mapping between proliferation events and the number of
        //! Poisson neutral mutations for this population of cells.
        ensure!(!cells.is_empty());
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.shrink_to(cells.len());
        let mut nb_variants = 0usize;
        // self.counts.shrink_to(cells.len());
        for (i, cell) in cells.iter().enumerate() {
            let printiter = verbosity > 2 && i % 1000 == 0;
            if printiter {
                println!("{}th cell", i);
            }
            // retrieve divisions assigned to this cell during the simulation
            for division_id in &cell.proliferation_events_id {
                let mutation_id = counts.entry(*division_id);
                let mut_number = match &mutation_id {
                    // if it's the first time we encounter this mutation id
                    // in the population, we assign to this mut id some poisson
                    // mutations
                    Entry::Vacant(_) => {
                        if printiter {
                            println!("vacant entry");
                        }
                        let poisson_nb = poisson_dist.nb_neutral_mutations(rng);
                        nb_variants += poisson_nb as usize;
                        poisson_nb as u64
                    }

                    // else we retrieve the poisson mutations already assigned
                    Entry::Occupied(entry) => {
                        if printiter {
                            println!("occupied entry {:#?}", &entry);
                        }
                        entry.get().poisson_mut_number
                    }
                };
                // we update the counter of cells by adding this cell to the
                // mutation id. We also update the counter of mutations
                mutation_id
                    .and_modify(|stat| {
                        let new_cell_count = 1 + stat.cell_count.get();
                        stat.cell_count = unsafe { NonZeroU64::new_unchecked(new_cell_count) };
                        stat.poisson_mut_number = mut_number;
                    })
                    .or_insert(CountsMutations {
                        cell_count: unsafe { NonZeroU64::new_unchecked(1) },
                        poisson_mut_number: mut_number,
                    });
            }
        }

        if verbosity > 1 {
            println!("total nb of variants: {}", nb_variants);
        }

        Ok(StatisticsMutations {
            counts,
            nb_variants,
        })
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
    proliferation_events_id: FxHashSet<GenotypeId>,
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
            proliferation_events_id: FxHashSet::default(),
        }
    }

    pub fn with_set_of_mutations(mutation_set: Vec<usize>) -> StemCell {
        //! Create a stem cell with a set of neutral mutations.
        //!
        //! The mutation set is not a collection of mutations but a collection
        //! of genotypes that will be converted later on into mutations, see
        //! [`StatisticsMutations`].
        let mut cell = StemCell::new();
        let mutation_set = HashSet::from_iter(mutation_set.into_iter());
        cell.proliferation_events_id = mutation_set;
        cell
    }

    pub fn has_mutations(&self) -> bool {
        !self.proliferation_events_id.is_empty()
    }

    pub fn record_division(&mut self, event_id: usize) {
        //! Record the iteration when the cell has undergone cell-division
        //! (proliferation).
        //!
        //! Since we store the divisions performed by each cell, instead of the
        //! mutations (to avoid generating a random number at every division?),
        //! we mutate a cell by storing the iterations at which the cell
        //! proliferates.
        //! The mapping between genotypes (i.e. proliferation id) and mutations
        //! is performed by [`StatisticsMutations`].
        self.proliferation_events_id.insert(event_id);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck::{Arbitrary, Gen};
    use quickcheck_macros::quickcheck;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;
    use rand_distr::Poisson;
    use std::num::NonZeroU8;

    #[quickcheck]
    fn genotype_mutate_test(divisions: usize) -> bool {
        let mut cell = StemCell::new();
        cell.record_division(divisions);
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

            DistinctMutationsId(mut1.get(), mut2.get(), mut3.get())
        }
    }

    #[test]
    #[should_panic]
    fn test_update_from_empty_cells() {
        let rng = &mut ChaCha8Rng::seed_from_u64(26);
        let poisson = &NeutralMutationPoisson(Poisson::new(1.).unwrap());
        StatisticsMutations::from_cells(&[], poisson, rng, 0).unwrap();
    }

    #[quickcheck]
    fn test_update_from_cell_without_mutations(lambda: NonZeroU8, seed: u64) -> bool {
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);
        let poisson = &NeutralMutationPoisson(Poisson::new(lambda.get() as f32).unwrap());
        let cells = &[StemCell::default()];
        let stats = StatisticsMutations::from_cells(cells, poisson, rng, 0).unwrap();
        stats.counts.is_empty()
    }

    #[quickcheck]
    fn test_update_from_cell_with_mutations(lambda: NonZeroU8, seed: u64) -> bool {
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);
        let poisson = &NeutralMutationPoisson(Poisson::new(lambda.get() as f32).unwrap());
        let cells = &[StemCell::with_set_of_mutations(vec![1])];
        let stats = StatisticsMutations::from_cells(cells, poisson, rng, 0).unwrap();
        stats
            .counts
            .into_values()
            .map(|ele| ele.cell_count.get())
            .collect::<Vec<u64>>()
            == vec![1]
    }

    #[quickcheck]
    fn test_stats_muts_from_cell_with_mutations_two_cells_same_genotype(
        lambda: NonZeroU8,
        seed: u64,
    ) -> bool {
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);
        let poisson = &NeutralMutationPoisson(Poisson::new(lambda.get() as f32).unwrap());
        let cells = &[
            StemCell::with_set_of_mutations(vec![1]),
            StemCell::with_set_of_mutations(vec![1]),
        ];
        let stats = StatisticsMutations::from_cells(cells, poisson, rng, 0).unwrap();
        stats
            .counts
            .into_values()
            .map(|ele| ele.cell_count.get())
            .collect::<Vec<u64>>()
            == vec![2]
    }

    #[quickcheck]
    fn test_stats_muts_from_cell_with_mutations_two_cells(lambda: NonZeroU8, seed: u64) -> bool {
        // might fail when the random number generates the same number of mutations
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);
        let poisson = &NeutralMutationPoisson(Poisson::new(lambda.get() as f32).unwrap());
        let cells = &[
            StemCell::with_set_of_mutations(vec![1]),
            StemCell::with_set_of_mutations(vec![2]),
        ];
        let stats = StatisticsMutations::from_cells(cells, poisson, rng, 0).unwrap();
        stats
            .counts
            .into_values()
            .map(|ele| ele.cell_count.get())
            .collect::<Vec<u64>>()
            == vec![1, 1]
    }

    #[test]
    #[should_panic]
    fn test_sfs_neutral_no_cells() {
        let verbosity = 1;
        let stats = StatisticsMutations {
            counts: FxHashMap::default(),
            nb_variants: 0,
        };
        Sfs::from_stats(&stats, verbosity).unwrap();
    }

    #[quickcheck]
    fn test_sfs_neutral_cell_no_mutations(cell_count: NonZeroU64) -> bool {
        let verbosity = 1;
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.insert(
            0,
            CountsMutations {
                cell_count,
                poisson_mut_number: 0,
            },
        );
        let stats = StatisticsMutations {
            counts,
            nb_variants: 0,
        };
        Sfs::from_stats(&stats, verbosity)
            .unwrap()
            .0
            .into_iter()
            .all(|(jmut, jcell)| jcell == cell_count.get() && jmut == 0)
    }

    #[quickcheck]
    fn test_sfs_one_proliferation_event(idx: DistinctMutationsId) -> bool {
        let random_nbs: [u64; 6] = std::array::from_fn(|i| idx.0 as u64 + i as u64);
        dbg!(&random_nbs);
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.insert(
            random_nbs[0] as usize,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(random_nbs[1]) },
                poisson_mut_number: random_nbs[2],
            },
        );
        let stats = StatisticsMutations {
            counts,
            nb_variants: random_nbs[2] as usize,
        };
        let sfs = Sfs::from_stats(&stats, 0).unwrap();
        let mut keys = sfs.0.clone().into_keys().collect::<Vec<u64>>();
        let mut values = sfs.0.into_values().collect::<Vec<u64>>();
        keys.sort_unstable();
        values.sort_unstable();

        let mut expected_jcells = vec![random_nbs[1]];
        let mut expected_jmuts = vec![random_nbs[2]];
        expected_jcells.sort_unstable();
        expected_jmuts.sort_unstable();

        keys == expected_jmuts && values == expected_jcells
    }

    #[quickcheck]
    fn test_sfs_two_proliferation_events_all_in_common(idx: DistinctMutationsId) -> bool {
        let random_nbs: [u64; 6] = std::array::from_fn(|i| idx.0 as u64 + i as u64);
        dbg!(&random_nbs);
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.insert(
            random_nbs[0] as usize,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(random_nbs[1]) },
                poisson_mut_number: random_nbs[2],
            },
        );
        counts.insert(
            random_nbs[3] as usize,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(random_nbs[4]) },
                poisson_mut_number: random_nbs[2],
            },
        );
        let stats = StatisticsMutations {
            counts,
            nb_variants: (random_nbs[2] + random_nbs[2]) as usize,
        };
        let sfs = Sfs::from_stats(&stats, 0).unwrap();
        let mut keys = sfs.0.clone().into_keys().collect::<Vec<u64>>();
        let mut values = sfs.0.into_values().collect::<Vec<u64>>();
        keys.sort_unstable();
        values.sort_unstable();

        let mut expected_jcells = vec![random_nbs[1] + random_nbs[4]];
        let mut expected_jmuts = vec![random_nbs[2]];
        expected_jcells.sort_unstable();
        expected_jmuts.sort_unstable();

        keys == expected_jmuts && values == expected_jcells
    }

    #[quickcheck]
    fn test_sfs_two_proliferation_events_nothing_in_common(idx: DistinctMutationsId) -> bool {
        let random_nbs: [u64; 6] = std::array::from_fn(|i| idx.0 as u64 + i as u64);
        dbg!(&random_nbs);
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.insert(
            random_nbs[0] as usize,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(random_nbs[1]) },
                poisson_mut_number: random_nbs[2],
            },
        );
        counts.insert(
            random_nbs[3] as usize,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(random_nbs[4]) },
                poisson_mut_number: random_nbs[5],
            },
        );
        let stats = StatisticsMutations {
            counts,
            nb_variants: (random_nbs[2] + random_nbs[5]) as usize,
        };
        let sfs = Sfs::from_stats(&stats, 0).unwrap();
        let mut keys = sfs.0.clone().into_keys().collect::<Vec<u64>>();
        let mut values = sfs.0.into_values().collect::<Vec<u64>>();
        keys.sort_unstable();
        values.sort_unstable();

        let mut expected_jcells = vec![random_nbs[1], random_nbs[4]];
        let mut expected_jmuts = vec![random_nbs[2], random_nbs[5]];
        expected_jcells.sort_unstable();
        expected_jmuts.sort_unstable();

        keys == expected_jmuts && values == expected_jcells
    }
}
