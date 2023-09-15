use anyhow::Context;
use rand::Rng;
use rand_distr::{Distribution, Poisson};
use rustc_hash::FxHashMap;
use std::{fs, path::Path};

use crate::stemcell::StemCell;

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
/// implemented as mapping with keys being j cells (x-axis) and values being
/// the number of variants with j cells (y-axis).
pub struct Sfs(pub FxHashMap<u64, u64>);

impl Sfs {
    pub fn from_cells(cells: &[StemCell], verbosity: u8) -> anyhow::Result<Self> {
        //! Compute the SFS from the stem cell population.
        let mut sfs_variants = FxHashMap::default();
        for cell in cells.iter() {
            for variant in cell.variants.iter() {
                sfs_variants
                    .entry(variant)
                    .and_modify(|counter| *counter += 1u64)
                    .or_insert(1);
            }
        }
        let mut sfs = FxHashMap::default();
        for nb_cells in sfs_variants.values() {
            sfs.entry(*nb_cells)
                .and_modify(|counter| *counter += 1)
                .or_insert(1u64);
        }

        if verbosity > 0 {
            println!("sfs: {:#?}", sfs);
        }
        Ok(Sfs(sfs))
    }

    pub fn save(&self, path2file: &Path) -> anyhow::Result<()> {
        let path2file = path2file.with_extension("json");
        let sfs = serde_json::to_string(&self.0).with_context(|| "cannot serialize the SFS")?;
        fs::write(path2file, sfs).with_context(|| "Cannot save the SFS ".to_string())?;

        Ok(())
    }
}

/// Single-cell mutational burden is a mapping of cells sharing a number of
/// mutations.
///
/// To plot it, plot on the x-axis the keys (number of mutations) and on the
/// y-axis the keys (the number of cells with those mutations).
pub struct MutationalBurden(pub FxHashMap<u16, u64>);

impl MutationalBurden {
    pub fn from_cells(cells: &[StemCell], verbosity: u8) -> anyhow::Result<Self> {
        //! Compute the single-cell mutational burden from the stem cell
        //! population.
        let mut burden: FxHashMap<u16, u64> = FxHashMap::default();
        for cell in cells.iter() {
            burden
                .entry(cell.burden() as u16)
                .and_modify(|counter| *counter += 1)
                .or_insert(1u64);
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
