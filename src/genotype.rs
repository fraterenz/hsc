use anyhow::{ensure, Context};
use rand::Rng;
use rand_distr::{Distribution, Poisson};
use rustc_hash::FxHashMap;
use std::{fs, path::Path};
use uuid::Uuid;

use crate::{
    process::{Exponential, Moran},
    stemcell::StemCell,
};

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
    pub lambda_background: f32,
    pub lambda_division: f32,
    division: Poisson<f32>,
}

impl PartialEq for NeutralMutationPoisson {
    fn eq(&self, other: &Self) -> bool {
        (self.lambda_division - other.lambda_division).abs() < f32::EPSILON
            && (self.lambda_background - other.lambda_background).abs() < f32::EPSILON
    }
}

impl Default for NeutralMutationPoisson {
    fn default() -> Self {
        NeutralMutationPoisson::new(1., 1.).unwrap()
    }
}

impl NeutralMutationPoisson {
    pub fn new(lambda_division: f32, lambda_background: f32) -> anyhow::Result<Self> {
        //! Create two Poisson distributions, one modelling the neutral
        //! mutations acquired upon cell-division and the other modelling the
        //! acquisition of neutral background mutations, i.e. all mutations
        //! occuring not during cell-division.
        ensure!(lambda_division > 0., "invalid value of lambda_division");
        ensure!(lambda_background > 0., "invalid value of lambda_background");
        Ok(Self {
            lambda_background,
            lambda_division,
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

    pub fn new_muts_background(
        &self,
        interdivison_time: f32,
        rng: &mut impl Rng,
        verbosity: u8,
    ) -> Option<Vec<Variant>> {
        //! The number of neutral mutations acquired upon cell division.
        if verbosity > 1 {
            println!(
                "interdivison_time = {} and background lambda {}",
                interdivison_time, self.lambda_background
            );
        }
        if interdivison_time > 0.001 {
            let background = Poisson::new(self.lambda_background * interdivison_time).unwrap();
            let nb_mutations = nb_neutral_mutations(&background, rng);
            if verbosity > 1 {
                println!("{} background mutations", nb_mutations);
            }
            generate_mutations(nb_mutations)
        } else {
            None
        }
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
/// implemented as mapping with keys being j cells (x-axis) and values being
/// the number of variants with j cells (y-axis).
pub struct Sfs(pub FxHashMap<u64, u64>);

impl Sfs {
    pub fn from_cells(cells: &[&StemCell], verbosity: u8) -> anyhow::Result<Self> {
        //! Compute the SFS from the stem cell population.
        if verbosity > 0 {
            println!("computing the SFS from {} cells", cells.len());
            if verbosity > 2 {
                println!("computing the SFS from {:#?}", &cells);
            }
        }
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

    pub fn save(&self, path2file: &Path, verbosity: u8) -> anyhow::Result<()> {
        let path2file = path2file.with_extension("json");
        let sfs = serde_json::to_string(&self.0).with_context(|| "cannot serialize the SFS")?;
        if verbosity > 0 {
            println!("SFS in {:#?}", path2file)
        }
        fs::write(path2file, sfs).with_context(|| "Cannot save the SFS ".to_string())?;

        Ok(())
    }
}

/// Single-cell mutational burden is a mapping of cells sharing a number of
/// mutations.
///
/// To plot it, plot on the x-axis the keys (number of mutations) and on the
/// y-axis the values (the number of cells with those mutations).
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

    pub fn from_moran(moran: &Moran, verbosity: u8) -> anyhow::Result<Self> {
        //! Create the single-cell mutational burden from all cells in the Moran
        //! process
        let cells: Vec<&StemCell> = moran
            .subclones
            .get_cells_with_clones_idx()
            .iter()
            .map(|ele| ele.0)
            .collect();
        MutationalBurden::from_cells(&cells, verbosity)
    }

    pub fn from_exp(exp: &Exponential, verbosity: u8) -> anyhow::Result<Self> {
        //! Create the single-cell mutational burden from all cells in the Moran
        //! process
        let cells: Vec<&StemCell> = exp
            .subclones
            .get_cells_with_clones_idx()
            .iter()
            .map(|ele| ele.0)
            .collect();
        MutationalBurden::from_cells(&cells, verbosity)
    }

    pub fn save(&self, path2file: &Path, verbosity: u8) -> anyhow::Result<()> {
        let path2file = path2file.with_extension("json");
        let burden =
            serde_json::to_string(&self.0).with_context(|| "cannot serialize the burden")?;
        if verbosity > 0 {
            println!("saving burden in {:#?}", path2file);
        }
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

    use crate::tests::LambdaFromNonZeroU8;

    use super::*;

    #[quickcheck]
    fn poisson_neutral_mutations(seed: u64) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let poisson = Poisson::new(0.0001).unwrap();
        NeutralMutationPoisson {
            lambda_background: 0.0001,
            lambda_division: 0.0001,
            division: poisson,
        }
        .new_muts_upon_division(&mut rng)
        .is_none()
    }

    #[test]
    fn test_sfs_4cells_6mutations() {
        // variants
        let square = Uuid::new_v4();
        let circle = Uuid::new_v4();
        let star = Uuid::new_v4();
        let triangle = Uuid::new_v4();
        let diamond = Uuid::new_v4();
        let thunder = Uuid::new_v4();
        // cells with variants
        let cell1 = StemCell::with_mutations(vec![square, circle, star]);
        let cell2 = StemCell::with_mutations(vec![square, circle, diamond]);
        let cell3 = StemCell::with_mutations(vec![square, triangle, star]);
        let cell4 = StemCell::with_mutations(vec![square, thunder]);

        // compute sfs and sort it by jcells
        let mut sfs = Sfs::from_cells(&[&cell1, &cell2, &cell3, &cell4], 0)
            .unwrap()
            .0
            .into_iter()
            .collect::<Vec<(u64, u64)>>();
        sfs.sort_unstable_by(|&entry1, &entry2| entry1.0.cmp(&entry2.0));
        let jcells = sfs
            .clone()
            .into_iter()
            .map(|ele| ele.0)
            .collect::<Vec<u64>>();
        let jmuts = sfs.into_iter().map(|ele| ele.1).collect::<Vec<u64>>();

        assert_eq!(jcells, [1, 2, 4]);
        assert_eq!(jmuts, [3, 2, 1]);
    }

    #[test]
    fn test_sfs_4cells_6mutations_with_variant_with_3cells() {
        // variants
        let square = Uuid::new_v4();
        let circle = Uuid::new_v4();
        let star = Uuid::new_v4();
        let triangle = Uuid::new_v4();
        let diamond = Uuid::new_v4();
        let thunder = Uuid::new_v4();
        // cells with variants
        let cell1 = StemCell::with_mutations(vec![square, circle, star]);
        let cell2 = StemCell::with_mutations(vec![square, circle, diamond]);
        let cell3 = StemCell::with_mutations(vec![square, triangle, star]);
        let cell4 = StemCell::with_mutations(vec![square, circle, thunder]);

        // compute sfs and sort it by jcells
        let mut sfs = Sfs::from_cells(&[&cell1, &cell2, &cell3, &cell4], 0)
            .unwrap()
            .0
            .into_iter()
            .collect::<Vec<(u64, u64)>>();
        sfs.sort_unstable_by(|&entry1, &entry2| entry1.0.cmp(&entry2.0));
        let jcells = sfs
            .clone()
            .into_iter()
            .map(|ele| ele.0)
            .collect::<Vec<u64>>();
        let jmuts = sfs.into_iter().map(|ele| ele.1).collect::<Vec<u64>>();

        assert_eq!(jcells, [1, 2, 3, 4]);
        assert_eq!(jmuts, [3, 1, 1, 1]);
    }

    #[quickcheck]
    fn partial_eq_neutral_poisson_test(
        lambda_division: LambdaFromNonZeroU8,
        lambda_background: LambdaFromNonZeroU8,
    ) -> bool {
        let poissons = NeutralMutationPoisson::new(lambda_division.0, lambda_background.0).unwrap();
        poissons == NeutralMutationPoisson::new(lambda_division.0, lambda_background.0).unwrap()
    }
}
