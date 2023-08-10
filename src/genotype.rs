use std::{
    collections::hash_map::Entry,
    fs,
    num::{NonZeroU16, NonZeroU64},
    path::Path,
};

use anyhow::{ensure, Context};
use rand::Rng;
use rustc_hash::FxHashMap;

use crate::{process::NeutralMutationPoisson, stemcell::StemCell};

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
pub struct Sfs(pub FxHashMap<NonZeroU16, u64>);

impl Sfs {
    pub fn from_stats(stats: &StatisticsMutations, verbosity: u8) -> anyhow::Result<Self> {
        //! Compute the site frequency spectrum (SFS) of mutations found in the
        //! stem cell population.
        //!
        //! The SFS is implemented with keys being the jmuts and the values
        //! being the jcells with with jmuts.
        // let mut sfs = HashMap::with_capacity(stats.nb_variants);
        ensure!(!stats.counts.is_empty(), "found empty stats");
        let mut sfs: FxHashMap<NonZeroU16, u64> = FxHashMap::default();
        sfs.shrink_to(stats.nb_variants);
        for stats_mut_id in stats.counts.values() {
            if stats_mut_id.poisson_mut_number > 0 {
                sfs.entry(unsafe {
                    NonZeroU16::new_unchecked(stats_mut_id.poisson_mut_number as u16)
                })
                .and_modify(|jcells| *jcells += stats_mut_id.cell_count.get())
                .or_insert(stats_mut_id.cell_count.get());
            }
        }
        if verbosity > 1 {
            println!("sfs: {:#?}", sfs);
        }
        Ok(Sfs(sfs))
    }

    pub fn from_cells(
        cells: &[StemCell],
        stats: &mut StatisticsMutations,
        poisson_dist: &NeutralMutationPoisson,
        rng: &mut impl Rng,
        verbosity: u8,
    ) -> anyhow::Result<Self> {
        //! Compute the site frequency spectrum (SFS) of mutations found in the
        //! stem cell population.
        //! Use this function when the cells used to construct the `stats` are
        //! not the same as the `cells` given here as the argument of the
        //! function.
        //!
        //! The SFS is implemented with keys being the jmuts and the values
        //! being the jcells with with jmuts.
        // let mut sfs = HashMap::with_capacity(stats.nb_variants);
        ensure!(!stats.counts.is_empty(), "found empty stats");
        let mut sfs: FxHashMap<NonZeroU16, u64> = FxHashMap::default();
        sfs.shrink_to(stats.nb_variants);
        for cell in cells.iter() {
            for id in cell.proliferation_events_id.iter() {
                if let Some(counts) = stats.counts.get(id) {
                    if counts.poisson_mut_number > 0 {
                        sfs.entry(unsafe {
                            NonZeroU16::new_unchecked(counts.poisson_mut_number as u16)
                        })
                        .and_modify(|jcells| *jcells += 1)
                        .or_insert(1);
                    }
                } else {
                    stats.store_new_genotype_id(*id, poisson_dist, rng, false);
                    if stats.counts[id].poisson_mut_number > 0 {
                        sfs.insert(
                            unsafe {
                                NonZeroU16::new_unchecked(
                                    stats.counts[id].poisson_mut_number as u16,
                                )
                            },
                            stats.counts[id].cell_count.get(),
                        );
                    }
                }
            }
        }
        if verbosity > 1 {
            println!("sfs: {:#?}", sfs);
        }
        Ok(Sfs(sfs))
    }

    pub fn save(&self, path2file: &Path) -> anyhow::Result<()> {
        let path2file = path2file.with_extension("json");
        // some simulations might have only cells in the wild-type clone
        let sfs = serde_json::to_string(&self.0).with_context(|| "cannot serialize the sfs")?;
        fs::write(path2file, sfs).with_context(|| "Cannot save the total SFS".to_string())?;

        Ok(())
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
        let mut stats = StatisticsMutations {
            counts,
            nb_variants: 0,
        };
        // self.counts.shrink_to(cells.len());
        for (i, cell) in cells.iter().enumerate() {
            let printiter = verbosity > 2 && i % 1000 == 0;
            if printiter {
                println!("{}th cell", i);
            }
            // retrieve divisions assigned to this cell during the simulation
            for division_id in &cell.proliferation_events_id {
                stats.store_new_genotype_id(*division_id, poisson_dist, rng, printiter);
            }
        }

        if verbosity > 1 {
            println!("total nb of variants: {}", stats.nb_variants);
        }

        Ok(stats)
    }

    pub fn store_new_genotype_id(
        &mut self,
        division_id: GenotypeId,
        poisson_dist: &NeutralMutationPoisson,
        rng: &mut impl Rng,
        printiter: bool,
    ) {
        let mutation_id = self.counts.entry(division_id);
        let mut_number = match &mutation_id {
            // if it's the first time we encounter this mutation id
            // in the population, we assign to this mut id some poisson
            // mutations
            Entry::Vacant(_) => {
                if printiter {
                    println!("vacant entry");
                }
                let poisson_nb = poisson_dist.nb_neutral_mutations(rng);
                self.nb_variants += poisson_nb as usize;
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

#[cfg(test)]
mod tests {
    use std::num::NonZeroU8;

    use quickcheck::{Arbitrary, Gen};
    use quickcheck_macros::quickcheck;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;
    use rand_distr::Poisson;

    use super::*;

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
    fn test_stats_muts_from_empty_cells() {
        let rng = &mut ChaCha8Rng::seed_from_u64(26);
        let poisson = &NeutralMutationPoisson(Poisson::new(1.).unwrap());
        StatisticsMutations::from_cells(&[], poisson, rng, 0).unwrap();
    }

    #[quickcheck]
    fn test_stats_muts_from_cell_without_mutations(lambda: NonZeroU8, seed: u64) -> bool {
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);
        let poisson = &NeutralMutationPoisson(Poisson::new(lambda.get() as f32).unwrap());
        let cells = &[StemCell::default()];
        let stats = StatisticsMutations::from_cells(cells, poisson, rng, 0).unwrap();
        stats.counts.is_empty()
    }

    #[quickcheck]
    fn test_stats_muts_from_cell_with_mutations(lambda: NonZeroU8, seed: u64) -> bool {
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
        Sfs::from_stats(&stats, verbosity).unwrap().0.is_empty()
    }

    #[quickcheck]
    fn test_sfs_one_proliferation_event(idx: DistinctMutationsId) -> bool {
        let random_nbs: [u64; 6] = std::array::from_fn(|i| idx.0 as u64 + i as u64);
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
        let mut keys = sfs.0.clone().into_keys().collect::<Vec<NonZeroU16>>();
        let mut values = sfs.0.into_values().collect::<Vec<u64>>();
        keys.sort_unstable();
        values.sort_unstable();

        let mut expected_jcells = vec![random_nbs[1]];
        let mut expected_jmuts = vec![random_nbs[2]];
        expected_jcells.sort_unstable();
        expected_jmuts.sort_unstable();

        keys.into_iter()
            .map(|e| e.get() as u64)
            .collect::<Vec<u64>>()
            == expected_jmuts
            && values == expected_jcells
    }

    #[quickcheck]
    fn test_sfs_two_proliferation_events_all_in_common(idx: DistinctMutationsId) -> bool {
        let random_nbs: [u64; 6] = std::array::from_fn(|i| idx.0 as u64 + i as u64);
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
        let mut keys = sfs.0.clone().into_keys().collect::<Vec<NonZeroU16>>();
        let mut values = sfs.0.into_values().collect::<Vec<u64>>();
        keys.sort_unstable();
        values.sort_unstable();

        let mut expected_jcells = vec![random_nbs[1] + random_nbs[4]];
        let mut expected_jmuts = vec![NonZeroU16::new(random_nbs[2] as u16).unwrap()];
        expected_jcells.sort_unstable();
        expected_jmuts.sort_unstable();

        keys == expected_jmuts && values == expected_jcells
    }

    #[quickcheck]
    fn test_sfs_two_proliferation_events_nothing_in_common(idx: DistinctMutationsId) -> bool {
        let random_nbs: [NonZeroU16; 6] = std::array::from_fn(|i| unsafe {
            NonZeroU16::new_unchecked(idx.0 as u16 + i as u16 + 1)
        });
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.insert(
            random_nbs[0].get() as usize,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(random_nbs[1].get() as u64) },
                poisson_mut_number: random_nbs[2].get() as u64,
            },
        );
        counts.insert(
            random_nbs[3].get() as usize,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(random_nbs[4].get() as u64) },
                poisson_mut_number: random_nbs[5].get() as u64,
            },
        );
        let stats = StatisticsMutations {
            counts,
            nb_variants: (random_nbs[2].get() + random_nbs[5].get()) as usize,
        };
        let sfs = Sfs::from_stats(&stats, 0).unwrap();
        let mut keys = sfs.0.clone().into_keys().collect::<Vec<NonZeroU16>>();
        let mut values = sfs.0.into_values().collect::<Vec<u64>>();
        keys.sort_unstable();
        values.sort_unstable();

        let mut expected_jcells = vec![random_nbs[1].get() as u64, random_nbs[4].get() as u64];
        let mut expected_jmuts = vec![random_nbs[2], random_nbs[5]];
        expected_jcells.sort_unstable();
        expected_jmuts.sort_unstable();

        keys == expected_jmuts && values == expected_jcells
    }
}
