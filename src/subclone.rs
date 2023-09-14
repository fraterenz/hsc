use std::path::Path;

use crate::{process::Distributions, stemcell::StemCell, write2file, MAX_SUBCLONES};
use anyhow::{ensure, Context};
use rand::{seq::IteratorRandom, Rng};
use rand_distr::{Distribution, Gamma};
use sosa::ReactionRates;

/// Fitness models implemented so far.
#[derive(Clone, Debug)]
pub enum Fitness {
    /// According to the `Neutral` fitness model, all subclones have the same
    /// birth-rates `b0`, where `b0` is the birth-rate of the wild-type neutral
    /// clone.
    Neutral,
    /// According to the `Fixed` fitness model, the birth-rates for the
    /// non-wildtype clones are `b0(1+s)`, where `b0` is the birth-rate of the
    /// wild-type neutral clone.
    /// That is, all the subclones have the same proliferate rates except the
    /// wild-type one.
    Fixed { s: f32 },
    /// According to the `GammaSampled` fitness model, the birth-rates for the
    /// non-wildtype clones are `b0(1+s_i)`, where `b0` is the birth-rate of the
    /// wild-type neutral clone and `s_i` are sampled from the Gamma
    /// distribution.
    /// In this case, all subclones can have different proliferate rates.
    GammaSampled { shape: f32, scale: f32 },
}

/// Id of the [`SubClone`]s.
pub type CloneId = usize;

#[derive(Debug, Clone)]
/// A group of cells sharing the same genetic background with a specific
/// proliferation rate.
///
/// The main loop of the simulation delegates the proliferation of cells to
/// this structure, meaning that the `SubClone` will randomly pick one of its
/// cells and make it proliferate.
/// Upon proliferation, the cell can be assigned to a new clone with
/// probability `p` (see [`crate::process::CellDivisionProbabilities`]).
pub struct SubClone {
    cells: Vec<StemCell>,
    pub id: CloneId,
}

impl Iterator for SubClone {
    type Item = StemCell;

    fn next(&mut self) -> Option<Self::Item> {
        self.cells.pop()
    }
}

impl SubClone {
    pub fn new(id: CloneId, cell_capacity: usize) -> SubClone {
        SubClone {
            cells: Vec::with_capacity(cell_capacity),
            id,
        }
    }

    pub fn get_cells(&self) -> &[StemCell] {
        &self.cells
    }

    pub fn is_empty(&self) -> bool {
        self.get_cells().is_empty()
    }

    pub fn assign_cell(&mut self, cell: StemCell) {
        self.cells.push(cell);
    }

    pub fn random_cell(&mut self, rng: &mut impl Rng) -> anyhow::Result<StemCell> {
        ensure!(!self.cells.is_empty());
        Ok(self.cells.swap_remove(rng.gen_range(0..self.cells.len())))
    }

    pub fn cell_count(&self) -> u64 {
        self.cells.len() as u64
    }
}

pub fn assign(
    subclones: &mut SubClones,
    old_subclone_id: CloneId,
    cell: StemCell,
    distr: &Distributions,
    rng: &mut impl Rng,
    verbosity: u8,
) {
    //! Assign `cell` to a subclone.
    //!
    //! There are two scenarios based on whether the proliferation has
    //! generated a new fit variant:
    //!
    //! 1. if it's the case, then `cell` will be assigned to a new **empty**
    //! random clone with an id different from `old_subclone_id` (panics if
    //! there aren't any empty subclones left).
    //! 2. else, reassign `cell` to the old subclone with id `old_subclone_id`
    // assign cell to a new random clone if `assign` returned some cell
    let id = if distr.acquire_p_mutation(rng) {
        let mut rnd_clone_id = rng.gen_range(0..MAX_SUBCLONES);
        let mut counter = 0;
        // the new random clone cannot have `subclone_id` id and must be
        // empty
        while (rnd_clone_id == old_subclone_id
            || !subclones.get_clone(rnd_clone_id).unwrap().is_empty())
            && counter <= MAX_SUBCLONES
        {
            rnd_clone_id = rng.gen_range(0..MAX_SUBCLONES);
            counter += 1;
        }
        assert!(counter <= MAX_SUBCLONES, "max number of clones reached");
        rnd_clone_id
    } else {
        if verbosity > 2 {
            println!("no new fit variants");
        }
        old_subclone_id
    };
    if verbosity > 2 {
        println!("assgin {:#?} to clone {:#?}", cell, subclones.get_clone(id));
    }

    subclones.get_mut_clone_unchecked(id).assign_cell(cell);
}

/// A collection of subclones each having their proliferative advantage.
/// This collection contains the neutral clone as well.
#[derive(Clone, Debug)]
pub struct SubClones([SubClone; MAX_SUBCLONES]);

impl SubClones {
    pub fn new(cells: Vec<StemCell>, capacity: usize) -> Self {
        //! Returns all the newly initiated subclones by assigning all `cells`
        //! to the neutral clone.
        //!
        //! All clones have their own `capacity`.
        // initial state
        let mut subclones: [SubClone; MAX_SUBCLONES] =
            std::array::from_fn(|i| SubClone::new(i, capacity));
        for cell in cells {
            subclones[0].assign_cell(cell);
        }

        Self(subclones)
    }

    pub fn compute_tot_cells(&self) -> u64 {
        self.0.iter().map(|subclone| subclone.cell_count()).sum()
    }

    pub fn the_only_one_subclone_present(&self) -> Option<CloneId> {
        //! returns `None` if more than one subclone is present else the clone
        //! id.
        for subclone in self.0.iter() {
            if subclone.cell_count() == self.compute_tot_cells() {
                return Some(subclone.id);
            }
        }
        None
    }

    pub fn get_clone(&self, id: usize) -> Option<&SubClone> {
        self.0.get(id)
    }

    pub fn get_mut_clone_unchecked(&mut self, id: usize) -> &mut SubClone {
        &mut self.0[id]
    }

    pub fn get_neutral_clone(&self) -> &SubClone {
        &self.0[0]
    }

    pub fn gillespie_set_of_reactions(&self) -> [CloneId; MAX_SUBCLONES] {
        core::array::from_fn(|i| i)
    }

    pub fn get_cells(&self) -> Vec<&StemCell> {
        self.0
            .iter()
            .flat_map(|subclone| subclone.get_cells())
            .collect()
    }

    pub fn get_cells_subsampled(&self, nb_cells: usize, rng: &mut impl Rng) -> Vec<&StemCell> {
        self.0
            .iter()
            .flat_map(|subclone| subclone.get_cells())
            .choose_multiple(rng, nb_cells)
    }

    pub fn gillespie_rates(
        &self,
        fitness: &Fitness,
        b0: f32,
        rng: &mut impl Rng,
    ) -> ReactionRates<MAX_SUBCLONES> {
        //! Create the Gillespie reaction rates according to the `fitness`
        //! model with `b0` being the proliferative rate of the wild-type clone.
        match fitness {
            Fitness::Fixed { s } => ReactionRates(core::array::from_fn(|i| {
                if i == 0 {
                    b0
                } else {
                    b0 * (1. + s)
                }
            })),
            &Fitness::GammaSampled { shape, scale } => {
                let gamma = Gamma::new(shape, scale).unwrap();
                ReactionRates(core::array::from_fn(|i| {
                    if i == 0 {
                        b0
                    } else {
                        b0 * (1. + gamma.sample(rng))
                    }
                }))
            }
            &Fitness::Neutral => ReactionRates(core::array::from_fn(|_| b0)),
        }
    }
}

impl Default for SubClones {
    fn default() -> Self {
        //! Create new subclones each having one cell (this creates also the cells).
        let subclones = std::array::from_fn(|i| {
            let mut subclone = SubClone::new(i, 2);
            subclone.assign_cell(StemCell::new());
            subclone
        });
        Self(subclones)
    }
}

/// Number of cells in subclones.
pub struct Variants {}

impl Variants {
    pub fn variant_counts(subclones: &SubClones) -> [u64; MAX_SUBCLONES] {
        //! The total variant count is the number of cells in all subclones.
        //! ```
        //! use hsc::MAX_SUBCLONES;
        //! # use hsc::{stemcell::StemCell, process::{Moran}, subclone::Variants};
        //! // create a process with one cell in each `MAX_SUBCLONES` subclones
        //! let mut hsc = Moran::default();
        //!
        //! assert_eq!(Variants::variant_counts(&hsc.subclones), [1; MAX_SUBCLONES]);
        //! ```
        std::array::from_fn(|i| subclones.0[i].cell_count())
    }

    pub fn variant_fractions(subclones: &SubClones, tot_cells: u64) -> Vec<f32> {
        //! The proportion of cells in all subclones.
        //!
        //! ```
        //! use hsc::MAX_SUBCLONES;
        //! # use hsc::{stemcell::StemCell, process::{Moran}, subclone::Variants};
        //! // create a process with one cell in each `MAX_SUBCLONES` subclones
        //! let mut hsc = Moran::default();
        //!
        //! assert_eq!(
        //!     Variants::variant_fractions(&hsc.subclones, hsc.subclones.compute_tot_cells()),
        //!     [1. / MAX_SUBCLONES as f32; MAX_SUBCLONES]
        //! );
        //! ```
        Variants::variant_counts(subclones)
            .into_iter()
            .map(|frac| (frac as f32 / tot_cells as f32))
            .collect()
    }
}

pub fn save_variant_fraction(
    subclones: &SubClones,
    path2file: &Path,
    verbosity: u8,
) -> anyhow::Result<()> {
    let path2file = path2file.with_extension("csv");
    let total_variant_frac = Variants::variant_fractions(subclones, subclones.compute_tot_cells());
    if verbosity > 1 {
        println!("total variant fraction in {:#?}", path2file)
    }
    write2file(&total_variant_frac, &path2file, None, false)?;
    Ok(())
}

pub fn proliferating_cell(
    subclones: &mut SubClones,
    subclone_id: CloneId,
    verbosity: u8,
    rng: &mut impl Rng,
) -> StemCell {
    //! Determine which cells will proliferate by randomly selecting a cell
    //! from the subclone with id `subclone_id`.
    //! This will also remove the random selected cell from the subclone,
    //! hence subclones are borrowed mutably.
    if verbosity > 2 {
        println!(
            "a cell from clone {:#?} will divide",
            subclones.0[subclone_id]
        );
    }
    subclones.0[subclone_id]
        .random_cell(rng)
        .with_context(|| "found empty subclone")
        .unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck_macros::quickcheck;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;
    use rand_distr::Uniform;

    #[quickcheck]
    fn assign_cell_test(id: usize) -> bool {
        let mut neutral_clone = SubClone { cells: vec![], id };
        let cell = StemCell::new();
        assert!(neutral_clone.cells.is_empty());

        neutral_clone.assign_cell(cell);
        !neutral_clone.cells.is_empty()
    }

    #[quickcheck]
    fn proliferating_cell_test(seed: u64, mut subclone_id: u8) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let mut subclones = SubClones::default();
        if subclone_id >= subclones.0.len() as u8 {
            subclone_id = rng.sample(Uniform::new(0, subclones.0.len())) as u8;
        }
        let number_of_cells = subclones.0[subclone_id as usize].get_cells().len();
        let tot_cells = subclones.compute_tot_cells();

        let _ = proliferating_cell(&mut subclones, subclone_id as usize, 1, &mut rng);

        let subclones_have_lost_cell =
            Variants::variant_counts(&subclones).iter().sum::<u64>() < tot_cells;
        let subclone_has_lost_cell =
            subclones.0[subclone_id as usize].get_cells().len() == number_of_cells - 1;
        subclones_have_lost_cell && subclone_has_lost_cell
    }

    #[test]
    fn variant_fraction() {
        let subclones = SubClones::default();
        assert_eq!(
            Variants::variant_fractions(&subclones, subclones.compute_tot_cells()),
            &[1. / MAX_SUBCLONES as f32; MAX_SUBCLONES]
        );
    }

    #[quickcheck]
    fn division_no_new_clone(seed: u64, cells_present: u8) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let distr = Distributions::new(1., 0f64, 0);

        let cells = vec![StemCell::new(); cells_present as usize];
        let before_assignment = cells.len();
        let cell2assign = StemCell::new();
        let mut subclones = SubClones::new(cells, cells_present as usize + 1);

        assign(&mut subclones, 0, cell2assign, &distr, &mut rng, 0);
        subclones.0[0].cell_count() as usize == before_assignment + 1
    }

    #[quickcheck]
    fn division_new_clone(seed: u64, cells_present: u8) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let distr = Distributions::new(1., 1f64, 0);

        let cells = vec![StemCell::new(); cells_present as usize];
        let before_assignment = cells.len();
        let cell2assign = StemCell::new();
        let mut subclones = SubClones::new(cells, cells_present as usize + 1);

        assign(&mut subclones, 0, cell2assign, &distr, &mut rng, 0);
        subclones.0[0].cell_count() as usize == before_assignment
    }

    #[test]
    #[should_panic]
    fn assign_all_clones_occupied() {
        let mut rng = ChaCha8Rng::seed_from_u64(26);
        let distr = Distributions::new(1., 1f64, 0);

        let mut subclones = SubClones::default();
        let cell = StemCell::new();

        assign(&mut subclones, 0, cell, &distr, &mut rng, 0);
    }
}
