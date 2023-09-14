use anyhow::{ensure, Context};
use rand::Rng;
use rand_distr::{Distribution, Poisson};
use rustc_hash::FxHashMap;
use serde::Serialize;
use std::{
    collections::hash_map::Entry,
    fs,
    num::{NonZeroU16, NonZeroU64},
    path::Path,
};

use crate::{stemcell::StemCell, write2file};

/// The Poisson probability distribution modeling the appearance of neutral
/// mutations upon cell-division, aka neutral mutations are assumed to follow
/// a [Poisson point process](https://en.wikipedia.org/wiki/Poisson_point_process).
#[derive(Debug, Clone)]
pub struct NeutralMutationPoisson(pub Poisson<f32>);

impl NeutralMutationPoisson {
    pub fn nb_neutral_mutations(&self, rng: &mut impl Rng) -> NbPoissonMutations {
        //! The number of neutral mutations acquired upon cell division.
        let mut mutations = self.0.sample(rng);
        while mutations >= u8::MAX as f32 || mutations.is_sign_negative() || mutations.is_nan() {
            mutations = self.0.sample(rng);
        }
        mutations as u16
    }
}

// poisson: NeutralMutationPoisson(
//     Poisson::new(lambda_poisson).expect("Invalid lambda found: lambda <= 0 or nan"),
// ),

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
pub type NbPoissonMutations = u16;

/// [Site frequency spectrum](https://en.wikipedia.org/wiki/Allele_frequency_spectrum)
/// implemented as a vector of variants where each entry is the number of stem
/// cells carrying that variant.
pub struct Sfs(Vec<NonZeroU64>);

impl Sfs {
    pub fn from_cells(
        cells: &[StemCell],
        stats: &StatisticsMutations,
        verbosity: u8,
    ) -> anyhow::Result<Self> {
        //! Compute the SFS from the stem cell population.
        //! The SFS is implemented as a vector of variants where each entry is
        //! the number of stem cells carrying that variant.
        //! This vector is returned sorted in ascending way.
        //!
        //! The SFS will be computed using the mutations in cells only if those
        //! genotypes id are present in stats, otherwise they will be filtered
        //! out.
        ensure!(!stats.counts.is_empty(), "found empty stats");
        let mut sfs_genotypes: FxHashMap<usize, NonZeroU64> = FxHashMap::default();
        sfs_genotypes.shrink_to(stats.nb_variants as usize);
        for cell in cells.iter() {
            for id in cell.proliferation_events_id.iter() {
                sfs_genotypes
                    .entry(*id)
                    .and_modify(|cell_count| {
                        *cell_count = unsafe { NonZeroU64::new_unchecked(cell_count.get() + 1) };
                    })
                    .or_insert(unsafe { NonZeroU64::new_unchecked(1) });
            }
        }
        if verbosity > 1 {
            println!(
                "sfs considering all cells all variants: {:#?}",
                sfs_genotypes
            );
        }
        let mut sfs = Vec::new();
        for (id, nb_cells) in sfs_genotypes.into_iter() {
            if let Some(counts) = stats.counts.get(&id) {
                for _ in 0..counts.poisson_mut_number {
                    sfs.push(nb_cells);
                }
            }
        }
        sfs.sort_unstable();
        if verbosity > 1 {
            println!("sfs considering only variants in stats: {:#?}", sfs);
        }
        Ok(Sfs(sfs))
    }

    pub fn save(&self, path2file: &Path) -> anyhow::Result<()> {
        let path2file = path2file.with_extension("csv");
        write2file(&self.0, &path2file, None, false)?;

        Ok(())
    }
}

/// Single-cell mutational burden is a mapping of cells sharing a number of
/// mutations.
///
/// To plot it, plot on the x-axis the keys (number of mutations) and on the
/// y-axis the keys (the number of cells with those mutations).
pub struct MutationalBurden(pub FxHashMap<NonZeroU16, u64>);

impl MutationalBurden {
    pub fn from_cells(
        cells: &[StemCell],
        stats: &StatisticsMutations,
        verbosity: u8,
    ) -> anyhow::Result<Self> {
        //! Compute the single-cell mutational burden from the stem cell
        //! population.
        //! Use this function when the cells used to construct the `stats` are
        //! not the same as the `cells` given here as the argument of the
        //! function.
        // let mut burden = HashMap::with_capacity(stats.nb_variants);
        ensure!(!stats.counts.is_empty(), "found empty stats");
        let mut burden: FxHashMap<NonZeroU16, u64> = FxHashMap::default();
        burden.shrink_to(stats.nb_variants as usize);
        for cell in cells.iter() {
            let mut cell_burden = 0;
            for id in cell.proliferation_events_id.iter() {
                if let Some(counts) = stats.counts.get(id) {
                    cell_burden += counts.poisson_mut_number;
                }
            }

            if cell_burden > 0 {
                burden
                    .entry(unsafe { NonZeroU16::new_unchecked(cell_burden) })
                    .and_modify(|jcells| *jcells += 1)
                    .or_insert(1);
            }
        }
        if verbosity > 1 {
            println!("burden: {:#?}", burden);
        }
        Ok(MutationalBurden(burden))
    }

    pub fn from_cells_update_stats(
        cells: &[StemCell],
        stats: &mut StatisticsMutations,
        poisson_dist: &NeutralMutationPoisson,
        rng: &mut impl Rng,
        verbosity: u8,
    ) -> anyhow::Result<Self> {
        //! Compute the single-cell mutational burden from the stem cell
        //! population.
        //! Use this function when the cells used to construct the `stats` are
        //! not the same as the `cells` given here as the argument of the
        //! function.
        //!
        //! This updates also `stats` inserting new entries if `cells` carry
        //! mutations not seen before (i.e. not stored in `stats`).
        // let mut burden = HashMap::with_capacity(stats.nb_variants);
        ensure!(!stats.counts.is_empty(), "found empty stats");
        let mut burden: FxHashMap<NonZeroU16, u64> = FxHashMap::default();
        burden.shrink_to(stats.nb_variants as usize);
        for cell in cells.iter() {
            let mut cell_burden = 0;
            for id in cell.proliferation_events_id.iter() {
                if let Some(counts) = stats.counts.get(id) {
                    cell_burden += counts.poisson_mut_number;
                } else {
                    stats.store_new_genotype_id(*id, poisson_dist, rng, false);
                    cell_burden += stats.counts[id].poisson_mut_number;
                }
            }

            if cell_burden > 0 {
                burden
                    .entry(unsafe { NonZeroU16::new_unchecked(cell_burden) })
                    .and_modify(|jcells| *jcells += 1)
                    .or_insert(1);
            }
        }
        if verbosity > 1 {
            println!("burden: {:#?}", burden);
        }
        Ok(MutationalBurden(burden))
    }

    pub fn save(&self, path2file: &Path) -> anyhow::Result<()> {
        let path2file = path2file.with_extension("json");
        let burden =
            serde_json::to_string(&self.0).with_context(|| "cannot serialize the burden")?;
        fs::write(path2file, burden)
            .with_context(|| "Cannot save the total single cel burden".to_string())?;

        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Serialize)]
struct CountsMutations {
    // the number of cells carrying a mutation id
    cell_count: NonZeroU64,
    // the number of mutations corresponding to a mutation id
    poisson_mut_number: u16,
}

/// Mapping between [`GenotypeId`] and the number of mutations and the number
/// of cells.
#[derive(Debug, Default, Clone, PartialEq, Serialize)]
pub struct StatisticsMutations {
    counts: FxHashMap<GenotypeId, CountsMutations>,
    nb_variants: u64,
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

    pub fn from_stats_removing_recent_entries(
        self,
        t_considered: GenotypeId,
        verbosity: u8,
    ) -> Self {
        //! Return a copy of the stats with entries that are smaller or equal
        //! to `t_considered`.
        //!
        //! ## Panics
        //! Panics if that stats is empty.
        assert!(!self.counts.is_empty());
        let mut nb_variants = 0;
        let mut counts = FxHashMap::default();
        let total = self.counts.len();
        for (id, count) in self.counts.into_iter() {
            if id <= t_considered {
                nb_variants += count.poisson_mut_number as u64;
                counts.insert(id, count);
            }
        }
        if verbosity > 1 {
            println!(
                "remove {} entries from a total of {} entries with t_considered {}",
                total - counts.len(),
                total,
                t_considered
            );

            println!("nb of variants: {} and coutns: {:#?}", nb_variants, counts);
        }
        counts.shrink_to_fit();
        StatisticsMutations {
            counts,
            nb_variants,
        }
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
                self.nb_variants += poisson_nb as u64;
                poisson_nb
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

    pub fn save(&self, path2file: &Path) -> anyhow::Result<()> {
        fs::create_dir_all(path2file.parent().unwrap()).with_context(|| "Cannot create dir")?;
        let path2file = path2file.with_extension("csv");

        let stats = serde_json::to_string(&self).with_context(|| "cannot serialize stats")?;
        fs::write(path2file.clone(), stats)
            .with_context(|| format!("Cannot save stats to {:#?}", path2file))?;
        Ok(())
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

    use crate::tests::LambdaFromNonZeroU8;

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

            DistinctMutationsId(mut1.get() as u16, mut2.get() as u16, mut3.get() as u16)
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
    fn test_burden_neutral_no_cells() {
        let verbosity = 1;
        let mut stats = StatisticsMutations {
            counts: FxHashMap::default(),
            nb_variants: 0,
        };
        MutationalBurden::from_cells_update_stats(
            &[],
            &mut stats,
            &NeutralMutationPoisson(Poisson::new(1.).unwrap()),
            &mut ChaCha8Rng::seed_from_u64(64),
            verbosity,
        )
        .unwrap();
    }

    #[quickcheck]
    fn test_burden_neutral_cell_no_mutations(cell_count: NonZeroU64) -> bool {
        let verbosity = 1;
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.insert(
            0,
            CountsMutations {
                cell_count,
                poisson_mut_number: 0,
            },
        );
        let mut stats = StatisticsMutations {
            counts,
            nb_variants: 0,
        };
        MutationalBurden::from_cells_update_stats(
            &[StemCell::default()],
            &mut stats,
            &NeutralMutationPoisson(Poisson::new(1.).unwrap()),
            &mut ChaCha8Rng::seed_from_u64(64),
            verbosity,
        )
        .unwrap()
        .0
        .is_empty()
    }

    #[quickcheck]
    fn test_burden_one_proliferation_event(idx: DistinctMutationsId) -> bool {
        let random_nbs: [u64; 6] = std::array::from_fn(|i| idx.0 as u64 + i as u64);
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.insert(
            random_nbs[0] as usize,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(random_nbs[1]) },
                poisson_mut_number: random_nbs[2] as u16,
            },
        );
        let mut stats = StatisticsMutations {
            counts,
            nb_variants: random_nbs[2],
        };
        let mut cells = Vec::new();
        for _ in 0..random_nbs[1] {
            cells.push(StemCell::with_set_of_mutations(vec![
                random_nbs[0] as usize,
            ]));
        }
        let burden = MutationalBurden::from_cells_update_stats(
            &cells,
            &mut stats,
            &NeutralMutationPoisson(Poisson::new(1.).unwrap()),
            &mut ChaCha8Rng::seed_from_u64(64),
            1,
        )
        .unwrap();
        let mut keys = burden.0.clone().into_keys().collect::<Vec<NonZeroU16>>();
        let mut values = burden.0.into_values().collect::<Vec<u64>>();
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
    fn test_burden_two_proliferation_events_all_in_common(idx: DistinctMutationsId) -> bool {
        let random_nbs: [u64; 6] = std::array::from_fn(|i| idx.0 as u64 + i as u64);
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.insert(
            random_nbs[0] as usize,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(random_nbs[1]) },
                poisson_mut_number: random_nbs[2] as u16,
            },
        );
        counts.insert(
            random_nbs[3] as usize,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(random_nbs[4]) },
                poisson_mut_number: random_nbs[2] as u16,
            },
        );
        let mut stats = StatisticsMutations {
            counts,
            nb_variants: random_nbs[2] + random_nbs[2],
        };
        let mut cells = Vec::new();
        for _ in 0..random_nbs[1] {
            cells.push(StemCell::with_set_of_mutations(vec![
                random_nbs[0] as usize,
            ]));
        }
        for _ in 0..random_nbs[4] {
            cells.push(StemCell::with_set_of_mutations(vec![
                random_nbs[3] as usize,
            ]));
        }
        let burden = MutationalBurden::from_cells_update_stats(
            &cells,
            &mut stats,
            &NeutralMutationPoisson(Poisson::new(1.).unwrap()),
            &mut ChaCha8Rng::seed_from_u64(64),
            1,
        )
        .unwrap();
        let mut keys = burden.0.clone().into_keys().collect::<Vec<NonZeroU16>>();
        let mut values = burden.0.into_values().collect::<Vec<u64>>();
        keys.sort_unstable();
        values.sort_unstable();

        let mut expected_jcells = vec![random_nbs[1] + random_nbs[4]];
        let mut expected_jmuts = vec![NonZeroU16::new(random_nbs[2] as u16).unwrap()];
        expected_jcells.sort_unstable();
        expected_jmuts.sort_unstable();

        keys == expected_jmuts && values == expected_jcells
    }

    #[quickcheck]
    fn test_burden_two_proliferation_events_nothing_in_common(idx: DistinctMutationsId) -> bool {
        let random_nbs: [NonZeroU16; 6] =
            std::array::from_fn(|i| unsafe { NonZeroU16::new_unchecked(idx.0 + i as u16 + 1) });
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.insert(
            random_nbs[0].get() as usize,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(random_nbs[1].get() as u64) },
                poisson_mut_number: random_nbs[2].get(),
            },
        );
        counts.insert(
            random_nbs[3].get() as usize,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(random_nbs[4].get() as u64) },
                poisson_mut_number: random_nbs[5].get(),
            },
        );
        let mut stats = StatisticsMutations {
            counts,
            nb_variants: (random_nbs[2].get() + random_nbs[5].get()) as u64,
        };
        let mut cells = Vec::new();
        for _ in 0..random_nbs[1].get() {
            cells.push(StemCell::with_set_of_mutations(vec![
                random_nbs[0].get() as usize
            ]));
        }
        for _ in 0..random_nbs[4].get() {
            cells.push(StemCell::with_set_of_mutations(vec![
                random_nbs[3].get() as usize
            ]));
        }
        let burden = MutationalBurden::from_cells_update_stats(
            &cells,
            &mut stats,
            &NeutralMutationPoisson(Poisson::new(1.).unwrap()),
            &mut ChaCha8Rng::seed_from_u64(64),
            1,
        )
        .unwrap();
        let mut keys = burden.0.clone().into_keys().collect::<Vec<NonZeroU16>>();
        let mut values = burden.0.into_values().collect::<Vec<u64>>();
        keys.sort_unstable();
        values.sort_unstable();

        let mut expected_jcells = vec![random_nbs[1].get() as u64, random_nbs[4].get() as u64];
        let mut expected_jmuts = vec![random_nbs[2], random_nbs[5]];
        expected_jcells.sort_unstable();
        expected_jmuts.sort_unstable();

        keys == expected_jmuts && values == expected_jcells
    }

    #[quickcheck]
    fn test_burden_from_cells_nothing_new(
        idx: DistinctMutationsId,
        state: u64,
        lambda: NonZeroU8,
    ) -> bool {
        let random_nbs: [NonZeroU16; 6] =
            std::array::from_fn(|i| unsafe { NonZeroU16::new_unchecked(idx.0 + i as u16 + 1) });
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.insert(
            random_nbs[0].get() as usize,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(random_nbs[1].get() as u64) },
                poisson_mut_number: random_nbs[2].get(),
            },
        );
        counts.insert(
            random_nbs[3].get() as usize,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(random_nbs[4].get() as u64) },
                poisson_mut_number: random_nbs[5].get(),
            },
        );
        let mut stats = StatisticsMutations {
            counts,
            nb_variants: (random_nbs[2].get() + random_nbs[5].get()) as u64,
        };
        let stats_copy = stats.clone();
        let mut cells = Vec::new();
        for _ in 0..random_nbs[1].get() {
            cells.push(StemCell::with_set_of_mutations(vec![
                random_nbs[0].get() as usize
            ]));
        }
        for _ in 0..random_nbs[4].get() {
            cells.push(StemCell::with_set_of_mutations(vec![
                random_nbs[3].get() as usize
            ]));
        }
        let burden = MutationalBurden::from_cells_update_stats(
            &cells,
            &mut stats,
            &NeutralMutationPoisson(Poisson::new(1.).unwrap()),
            &mut ChaCha8Rng::seed_from_u64(64),
            1,
        )
        .unwrap();
        let mut keys = burden.0.clone().into_keys().collect::<Vec<NonZeroU16>>();
        let mut values = burden.0.into_values().collect::<Vec<u64>>();
        keys.sort_unstable();
        values.sort_unstable();
        // again (same cells, nothing has changed)
        let mut cells = vec![];
        for _ in 0..random_nbs[1].get() {
            cells.push(StemCell::with_set_of_mutations(vec![
                random_nbs[0].get() as usize
            ]));
        }
        for _ in 0..random_nbs[4].get() {
            cells.push(StemCell::with_set_of_mutations(vec![
                random_nbs[3].get() as usize
            ]));
        }
        let burden = MutationalBurden::from_cells_update_stats(
            &cells,
            &mut stats,
            &NeutralMutationPoisson(Poisson::new(lambda.get() as f32).unwrap()),
            &mut ChaCha8Rng::seed_from_u64(state),
            0,
        )
        .unwrap();
        let mut keys_again = burden.0.clone().into_keys().collect::<Vec<NonZeroU16>>();
        let mut values_again = burden.0.into_values().collect::<Vec<u64>>();
        keys_again.sort_unstable();
        values_again.sort_unstable();

        let mut expected_jcells = vec![random_nbs[1].get() as u64, random_nbs[4].get() as u64];
        let mut expected_jmuts = vec![random_nbs[2], random_nbs[5]];
        expected_jcells.sort_unstable();
        expected_jmuts.sort_unstable();

        keys == expected_jmuts
            && expected_jmuts == keys_again
            && values == expected_jcells
            && values == values_again
            && stats == stats_copy
    }

    #[quickcheck]
    fn test_sfs_two_proliferation_events_two_cells_some_in_common_subclonal(
        idx: DistinctMutationsId,
        verbosity: u8,
    ) -> bool {
        let random_nbs: [NonZeroU16; 4] =
            std::array::from_fn(|i| unsafe { NonZeroU16::new_unchecked(idx.0 + i as u16 + 1) });
        let genotype1 = random_nbs[0].get() as usize;
        let genotype2 = random_nbs[2].get() as usize;
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.insert(
            genotype1,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(2) },
                poisson_mut_number: random_nbs[1].get(),
            },
        );
        counts.insert(
            genotype2,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(2) },
                poisson_mut_number: random_nbs[3].get(),
            },
        );
        let stats = StatisticsMutations {
            counts,
            nb_variants: (random_nbs[1].get() + random_nbs[3].get()) as u64,
        };
        let cells = [
            StemCell::with_set_of_mutations(vec![genotype1, genotype2]),
            StemCell::with_set_of_mutations(vec![genotype1]),
        ];

        let sfs = Sfs::from_cells(&cells, &stats, verbosity).unwrap();
        let mut expect_sfs = vec![];
        for _ in 0..random_nbs[1].get() {
            expect_sfs.push(NonZeroU64::new(2).unwrap());
        }
        for _ in 0..random_nbs[3].get() {
            expect_sfs.push(NonZeroU64::new(1).unwrap());
        }
        expect_sfs.sort_unstable();

        sfs.0 == expect_sfs
    }

    #[quickcheck]
    fn test_sfs_two_proliferation_events_two_cells_all_in_common_clonal(
        idx: DistinctMutationsId,
        verbosity: u8,
    ) -> bool {
        let random_nbs: [NonZeroU16; 4] =
            std::array::from_fn(|i| unsafe { NonZeroU16::new_unchecked(idx.0 + i as u16 + 1) });
        let genotype1 = random_nbs[0].get() as usize;
        let genotype2 = random_nbs[2].get() as usize;
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.insert(
            genotype1,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(2) },
                poisson_mut_number: random_nbs[1].get(),
            },
        );
        counts.insert(
            genotype2,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(2) },
                poisson_mut_number: random_nbs[3].get(),
            },
        );
        let stats = StatisticsMutations {
            counts,
            nb_variants: (random_nbs[1].get() + random_nbs[3].get()) as u64,
        };
        let cells = [
            StemCell::with_set_of_mutations(vec![genotype1, genotype2]),
            StemCell::with_set_of_mutations(vec![genotype1, genotype2]),
        ];

        let sfs = Sfs::from_cells(&cells, &stats, verbosity).unwrap();
        let mut expect_sfs = vec![];
        for _ in 0..random_nbs[1].get() {
            expect_sfs.push(NonZeroU64::new(2).unwrap());
        }
        for _ in 0..random_nbs[3].get() {
            expect_sfs.push(NonZeroU64::new(2).unwrap());
        }
        expect_sfs.sort_unstable();

        sfs.0 == expect_sfs
    }

    #[quickcheck]
    fn test_sfs_two_proliferation_events_two_cells_nothing_in_common(
        idx: DistinctMutationsId,
        verbosity: u8,
    ) -> bool {
        let random_nbs: [NonZeroU16; 4] =
            std::array::from_fn(|i| unsafe { NonZeroU16::new_unchecked(idx.0 + i as u16 + 1) });
        let genotype1 = random_nbs[0].get() as usize;
        let genotype2 = random_nbs[2].get() as usize;
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.insert(
            genotype1,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(1) },
                poisson_mut_number: random_nbs[1].get(),
            },
        );
        counts.insert(
            genotype2,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(1) },
                poisson_mut_number: random_nbs[3].get(),
            },
        );
        let stats = StatisticsMutations {
            counts,
            nb_variants: (random_nbs[1].get() + random_nbs[3].get()) as u64,
        };
        let cells = [
            StemCell::with_set_of_mutations(vec![genotype1]),
            StemCell::with_set_of_mutations(vec![genotype2]),
        ];

        let sfs = Sfs::from_cells(&cells, &stats, verbosity).unwrap();
        let mut expect_sfs = vec![];
        for _ in 0..random_nbs[1].get() {
            expect_sfs.push(NonZeroU64::new(1).unwrap());
        }
        for _ in 0..random_nbs[3].get() {
            expect_sfs.push(NonZeroU64::new(1).unwrap());
        }
        expect_sfs.sort_unstable();

        sfs.0 == expect_sfs
    }

    #[quickcheck]
    fn test_from_stats_removing_recent_entries(idx: DistinctMutationsId) -> bool {
        let random_nbs: [NonZeroU16; 4] =
            std::array::from_fn(|i| unsafe { NonZeroU16::new_unchecked(idx.0 + i as u16 + 1) });
        let genotype1 = random_nbs[0].get() as usize;
        let genotype2 = random_nbs[2].get() as usize;
        let use_genotype1_as_t = genotype1 < genotype2;
        let t_considered = if use_genotype1_as_t {
            genotype1
        } else {
            genotype2
        };
        let mut counts: FxHashMap<GenotypeId, CountsMutations> = FxHashMap::default();
        counts.insert(
            genotype1,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(1) },
                poisson_mut_number: random_nbs[1].get(),
            },
        );
        counts.insert(
            genotype2,
            CountsMutations {
                cell_count: unsafe { NonZeroU64::new_unchecked(1) },
                poisson_mut_number: random_nbs[3].get(),
            },
        );
        let stats = StatisticsMutations {
            counts: counts.clone(),
            nb_variants: (random_nbs[1].get() + random_nbs[3].get()) as u64,
        };
        let stats = stats.from_stats_removing_recent_entries(t_considered, 1);
        if use_genotype1_as_t {
            stats.counts.get(&genotype2).is_none()
                && stats.nb_variants == random_nbs[1].get() as u64
                && *stats.counts.get(&genotype1).unwrap() == counts[&genotype1]
        } else {
            stats.counts.get(&genotype1).is_none()
                && stats.nb_variants == random_nbs[3].get() as u64
                && *stats.counts.get(&genotype2).unwrap() == counts[&genotype2]
        }
    }

    #[quickcheck]
    fn test_nb_neutral_mutations(lambda_poisson: LambdaFromNonZeroU8) {
        let mut rng = ChaCha8Rng::seed_from_u64(26);
        NeutralMutationPoisson(Poisson::new(lambda_poisson.0).unwrap())
            .nb_neutral_mutations(&mut rng);
    }

    #[quickcheck]
    fn poisson_neutral_mutations(seed: u64) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let poisson = Poisson::new(0.0001).unwrap();
        NeutralMutationPoisson(poisson).nb_neutral_mutations(&mut rng) < 2
    }
}
