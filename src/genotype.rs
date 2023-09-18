use anyhow::Context;
use rand::Rng;
use rand_distr::{Distribution, Poisson};
use rustc_hash::FxHashMap;
use std::{fs, path::Path};
use uuid::Uuid;

use crate::stemcell::StemCell;

pub type Variant = Uuid;

fn nb_neutral_mutations(poisson: &Poisson<f32>, rng: &mut impl Rng) -> NbPoissonMutations {
    //! The number of neutral mutations acquired upon cell division.
    let mut mutations = poisson.sample(rng);
    while mutations >= u16::MAX as f32 || mutations.is_sign_negative() || mutations.is_nan() {
        mutations = poisson.sample(rng);
    }
    mutations as NbPoissonMutations
}

/// The Poisson probability distribution modeling the appearance of neutral
/// mutations, which are assumed to follow a
/// [Poisson point process](https://en.wikipedia.org/wiki/Poisson_point_process).
#[derive(Debug, Clone)]
pub struct NeutralMutationPoisson {
    background: Poisson<f32>,
    division: Poisson<f32>,
}

impl NeutralMutationPoisson {
    pub fn new(lambda_division: f32, lambda_background: f32) -> anyhow::Result<Self> {
        //! Create two Poisson distributions, one modelling the neutral
        //! mutations acquired upon cell-division and the other modelling the
        //! acquisition of neutral background mutations, i.e. all mutations
        //! occuring not during cell-division.
        Ok(Self {
            background: Poisson::new(lambda_background)
                .with_context(|| {
                    format!(
                        "invalid value of lambda for the background mutations {}",
                        lambda_background
                    )
                })
                .unwrap(),
            division: Poisson::new(lambda_division)
                .with_context(|| {
                    format!(
                        "invalid value of lambda for the division mutations {}",
                        lambda_division
                    )
                })
                .unwrap(),
        })
    }

    pub fn new_muts_upon_division(&self, rng: &mut impl Rng) -> Option<Vec<Variant>> {
        //! Generate neutral mutations acquired upon cell division sampling the
        //! Poisson
        //!
        //! ## Returns
        //! Returns `None` when the number of mutations have
        let nb_mutations = nb_neutral_mutations(&self.division, rng);
        generate_mutations(nb_mutations)
    }

    pub fn new_muts_background(&self, rng: &mut impl Rng) -> Option<Vec<Variant>> {
        //! The number of neutral mutations acquired upon cell division.
        let nb_mutations = nb_neutral_mutations(&self.background, rng);
        generate_mutations(nb_mutations)
    }
}

/// The number of mutations that are produced by a division event.
/// We assume that a maximal number of 255 neutral mutations can be generated
/// upon one proliferative event.
type NbPoissonMutations = u16;

fn generate_mutations(nb_mutations: NbPoissonMutations) -> Option<Vec<Variant>> {
    if nb_mutations == 0 {
        None
    } else {
        Some((0..nb_mutations).map(|_| Uuid::new_v4()).collect())
    }
}

/// [Site frequency spectrum](https://en.wikipedia.org/wiki/Allele_frequency_spectrum)
/// implemented as mapping with values being j cells (x-axis) and keys being
/// the number of variants with j cells (y-axis).
pub struct Sfs(pub FxHashMap<u64, u64>);

impl Sfs {
    pub fn from_cells(cells: &[&StemCell], verbosity: u8) -> anyhow::Result<Self> {
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
    pub fn from_cells(cells: &[&StemCell], verbosity: u8) -> anyhow::Result<Self> {
        //! Compute the single-cell mutational burden from the stem cell
        //! population.
        let mut burden: FxHashMap<u16, u64> = FxHashMap::default();
        for cell in cells.iter() {
            burden
                .entry(cell.burden() as u16)
                .and_modify(|counter| *counter += 1)
                .or_insert(1u64);
        }

        if verbosity > 0 {
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
    use quickcheck_macros::quickcheck;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;
    use rand_distr::Poisson;

    use super::*;

    #[quickcheck]
    fn poisson_neutral_mutations(seed: u64) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let poisson = Poisson::new(0.0001).unwrap();
        NeutralMutationPoisson {
            background: poisson.clone(),
            division: poisson,
        }
        .new_muts_upon_division(&mut rng)
        .is_none()
    }
}
