use crate::genotype::NeutralMutationPoisson;
use crate::Probs;
use crate::{stemcell::StemCell, write2file, MAX_SUBCLONES};
use anyhow::{ensure, Context};
use log::{debug, trace};
use rand::seq::IteratorRandom;
use rand::Rng;
use rand_distr::{Bernoulli, Distribution, Gamma};
use sosa::ReactionRates;
use std::num::NonZeroUsize;
use std::path::Path;

/// Distribution probabilities for the simulations upon cell division.
#[derive(Debug, Default, Clone)]
pub struct Distributions {
    /// arrival rate of fit mutations per division per cell
    pub u: f32,
    /// neutral mutations
    pub neutral_poisson: NeutralMutationPoisson,
}

impl Distributions {
    pub fn new(probs: Probs) -> Self {
        let (u, background, division) = match probs {
            Probs::Symmetric { u, probs_per_year } => {
                (u, probs_per_year.mu_background, probs_per_year.mu_division)
            }
            Probs::Asymmetric {
                u, probs_per_year, ..
            } => (u, probs_per_year.mu_background, probs_per_year.mu_division),
        };
        debug!(
                "creating distributions with u: {u}, lambda_background: {background}, lambda_division: {division}"
            );

        Self {
            u,
            neutral_poisson: NeutralMutationPoisson::new(division, background).unwrap(),
        }
    }
}

pub fn from_mean_std_to_shape_scale(mean: f32, std: f32) -> (f32, f32) {
    (mean.powf(2.) / std.powf(2.), std.powf(2.) / mean)
}

/// Fitness models implemented so far.
#[derive(Clone, Copy, Debug)]
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

impl Fitness {
    pub fn get_mean_std(&self) -> (f32, f32) {
        match self {
            Fitness::Neutral => (0., 0.),
            Fitness::Fixed { s } => (*s, 0.),
            Fitness::GammaSampled { shape, scale } => {
                (shape * scale, f32::sqrt(shape * scale.powi(2)))
            }
        }
    }
}

/// Id of the [`SubClone`]s.
///
/// `u16` is wide enough for `MAX_SUBCLONES = 3000` (≪ `u16::MAX = 65535`) and
/// keeps `Option<CloneId>` at 4 bytes (vs 16 with `usize`), shrinking
/// `SubClone` accordingly. `usize` is reserved for array-indexing call sites,
/// which cast at the boundary.
pub type CloneId = u16;

fn default_parent_id(id: CloneId) -> Option<CloneId> {
    // wild-type has no parent; every other slot is, by default, a child of
    // the wild-type clone until a fit-mutation event reassigns it.
    if id == 0 {
        None
    } else {
        Some(0)
    }
}

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
    // `None` only for the wild-type clone (id 0), which has no parent.
    // Every other clone defaults to `Some(0)`: until a fit-mutation event
    // overwrites it, the slot is conceptually descended from the wild-type.
    parent_id: Option<CloneId>,
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
            parent_id: default_parent_id(id),
        }
    }

    pub fn with_capacity(id: CloneId, capacity: usize) -> SubClone {
        SubClone {
            cells: Vec::with_capacity(capacity),
            id,
            parent_id: default_parent_id(id),
        }
    }

    pub fn empty_with_id(id: CloneId) -> SubClone {
        SubClone {
            cells: vec![],
            id,
            parent_id: default_parent_id(id),
        }
    }

    pub fn get_mut_cells(&mut self) -> &mut [StemCell] {
        &mut self.cells
    }

    pub fn get_cells(&self) -> &[StemCell] {
        &self.cells
    }

    pub fn get_cells_subclones_idx(&self) -> Vec<(&StemCell, CloneId)> {
        self.cells.iter().map(|cell| (cell, self.id)).collect()
    }

    pub fn is_empty(&self) -> bool {
        self.get_cells().is_empty()
    }

    pub fn get_parent_id(&self) -> Option<CloneId> {
        self.parent_id
    }

    pub fn set_parent_id(&mut self, parent_id: CloneId) {
        self.parent_id = Some(parent_id);
    }

    pub fn assign_cell(&mut self, cell: StemCell) {
        self.cells.push(cell);
    }

    pub fn random_cell(&mut self, rng: &mut impl Rng) -> anyhow::Result<StemCell> {
        ensure!(!self.cells.is_empty());
        Ok(self
            .cells
            .swap_remove(rng.random_range(0..self.cells.len())))
    }

    pub fn cell_count(&self) -> u64 {
        self.cells.len() as u64
    }
}

pub fn next_clone(
    subclones: &SubClones,
    old_subclone_id: CloneId,
    p: f64,
    rng: &mut impl Rng,
) -> CloneId {
    //! Returns a [`CloneId`] to which the stem cell must be assigned by
    //! sampling a fit variant.
    //!
    //! This id can be either a new id if a new variant has been arisen, else
    //! it will be the old clone.
    // this will panic when u is not small and the cell hasn't divided much,
    // i.e. when p > 1
    if Bernoulli::new(p).unwrap().sample(rng) {
        let mut rnd_clone_id = rng.random_range(0..MAX_SUBCLONES) as CloneId;
        let mut counter = 0;
        // the new random clone cannot have `subclone_id` id and must be empty
        while (rnd_clone_id == old_subclone_id
            || !subclones.get_clone(rnd_clone_id).unwrap().is_empty())
            && counter <= MAX_SUBCLONES
        {
            rnd_clone_id = rng.random_range(0..MAX_SUBCLONES) as CloneId;
            counter += 1;
        }
        assert!(counter <= MAX_SUBCLONES, "max number of clones reached");
        debug!("new fit variant: assign cell to clone {rnd_clone_id}");
        return rnd_clone_id;
    }
    debug!("no new fit variants with p {p}",);
    old_subclone_id
}

/// A collection of subclones each having their proliferative advantage.
/// This collection contains the neutral clone as well.
#[derive(Clone, Debug)]
// Boxed slice rather than `[SubClone; MAX_SUBCLONES]` to avoid overflow in
// debug mode.
pub struct SubClones(Box<[SubClone]>);

impl SubClones {
    pub fn new(cells: Vec<StemCell>, capacity: usize) -> Self {
        //! Returns all the newly initiated subclones by assigning all `cells`
        //! to the neutral clone.
        //!
        //! All clones have their own `capacity`.
        // initial state
        let mut subclones: Box<[SubClone]> = (0..MAX_SUBCLONES)
            .map(|i| SubClone::new(i as CloneId, capacity))
            .collect();
        let tot_cells = cells.len();
        debug!("assigning {tot_cells} cells to the wild-type clone");
        for cell in cells {
            subclones[0].assign_cell(cell);
        }
        assert_eq!(subclones[0].cells.len(), tot_cells);

        Self(subclones)
    }

    pub fn new_empty() -> Self {
        SubClones(
            (0..MAX_SUBCLONES)
                .map(|i| SubClone::empty_with_id(i as CloneId))
                .collect(),
        )
    }

    pub fn with_capacity(capacity: usize) -> Self {
        SubClones(
            (0..MAX_SUBCLONES)
                .map(|id| SubClone::with_capacity(id as CloneId, capacity))
                .collect(),
        )
    }

    pub fn compute_tot_cells(&self) -> u64 {
        self.0.iter().map(|subclone| subclone.cell_count()).sum()
    }

    pub(crate) fn get_clone(&self, id: CloneId) -> Option<&SubClone> {
        self.0.get(id as usize)
    }

    pub(crate) fn get_mut_clone_unchecked(&mut self, id: CloneId) -> &mut SubClone {
        &mut self.0[id as usize]
    }

    pub fn get_neutral_clone(&self) -> &SubClone {
        &self.0[0]
    }

    pub fn array_of_gillespie_reactions(&self) -> [CloneId; MAX_SUBCLONES] {
        core::array::from_fn(|i| i as CloneId)
    }

    pub(crate) fn into_subsampled(self, nb_cells: NonZeroUsize, rng: &mut impl Rng) -> Self {
        //! Take a subsample of the population that is stored in these subclones.
        let parent_ids: [Option<CloneId>; MAX_SUBCLONES] =
            std::array::from_fn(|i| self.0[i].parent_id);

        let mut subclones = Self::new_empty();

        // `Box<[T]>::into_iter()` derefs to `[T]` and yields `&T`; round-trip
        // through `Vec` to consume the slice and get owned `SubClone`s.
        for (cell, clone_id) in Vec::from(self.0)
            .into_iter()
            .flat_map(|subclone| {
                // for all subclones get all the cells with the its subclone_id.
                // This will return all the population of cells, as a vec of
                // tuple with the first entry being a stem cell and the second
                // entry the id of the clone to which the cell belongs.
                // Then subsample this population, by getting a `nb_cells`
                // ammount of tuples (cell, id).
                subclone
                    .cells
                    .into_iter()
                    .map(|cell| (cell, subclone.id))
                    .collect::<Vec<(StemCell, CloneId)>>()
            })
            .choose_multiple(rng, nb_cells.get())
        {
            subclones.0[clone_id as usize].cells.push(cell);
        }

        // `parent_id` is preserved for clones that retain at least one cell
        // after subsampling: clones that lost all their cells revert to the
        // constructor default (`None` for the wild-type, `Some(0)` elsewhere).
        for (i, &pid) in parent_ids.iter().enumerate() {
            if !subclones.0[i].is_empty() {
                subclones.0[i].parent_id = pid;
            }
        }

        subclones
    }

    pub fn get_mut_cells(&mut self) -> Vec<&mut StemCell> {
        self.0
            .iter_mut()
            .flat_map(|subclone| subclone.get_mut_cells())
            .collect()
    }

    pub fn get_cells(&self) -> Vec<&StemCell> {
        self.0
            .iter()
            .flat_map(|subclone| subclone.get_cells())
            .collect()
    }

    pub(crate) fn get_cells_with_clones_idx(&self) -> Vec<(&StemCell, CloneId)> {
        self.0
            .iter()
            .flat_map(|subclone| subclone.get_cells_subclones_idx())
            .collect()
    }

    pub(crate) fn get_cells_subsampled_with_clones_idx(
        &self,
        nb_cells: usize,
        rng: &mut impl Rng,
    ) -> Vec<(&StemCell, CloneId)> {
        self.0
            .iter()
            .flat_map(|subclone| subclone.get_cells_subclones_idx())
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
        let subclones: Box<[SubClone]> = (0..MAX_SUBCLONES)
            .map(|i| {
                let mut subclone = SubClone::new(i as CloneId, 2);
                subclone.assign_cell(StemCell::new());
                subclone
            })
            .collect();
        Self(subclones)
    }
}

impl From<Vec<(StemCell, CloneId)>> for SubClones {
    fn from(cells: Vec<(StemCell, CloneId)>) -> Self {
        let mut subclones = SubClones::new_empty();
        for (cell, id) in cells.into_iter() {
            subclones.get_mut_clone_unchecked(id).assign_cell(cell);
        }
        subclones
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

    pub fn variant_fractions(subclones: &SubClones) -> Vec<f32> {
        //! The proportion of cells in all subclones.
        //!
        //! ```
        //! use hsc::MAX_SUBCLONES;
        //! use hsc::subclone::SubClones;
        //! # use hsc::{stemcell::StemCell, process::{Moran}, subclone::Variants};
        //! // create a process with one cell in each `MAX_SUBCLONES` subclones
        //! let mut hsc = Moran::default();
        //!
        //! assert_eq!(
        //!     Variants::variant_fractions(&hsc.subclones),
        //!     [1. / MAX_SUBCLONES as f32; MAX_SUBCLONES]
        //! );
        //!
        //! let subclones = SubClones::new(vec![StemCell::new()], MAX_SUBCLONES);
        //! let mut variant_fraction = [0.; MAX_SUBCLONES];
        //! variant_fraction[0] = 1.;
        //!
        //! assert_eq!(
        //!     Variants::variant_fractions(&subclones),
        //!     variant_fraction
        //! );
        //! ```
        let tot_cells = subclones.compute_tot_cells();
        Variants::variant_counts(subclones)
            .into_iter()
            .map(|frac| frac as f32 / tot_cells as f32)
            .collect()
    }
}

pub fn save_variant_fraction(subclones: &SubClones, path2file: &Path) -> anyhow::Result<()> {
    let path2file = path2file.with_extension("csv");
    let total_variant_frac = Variants::variant_fractions(subclones);
    debug!("total variant fraction in {path2file:#?}");
    write2file(&total_variant_frac, &path2file, None, false)?;
    Ok(())
}

/// Display wrapper around `Option<CloneId>` so the value can flow through
/// [`write2file`]. Serializes as the bare integer for `Some(id)` and the
/// literal `None` for the wild-type root.
struct ParentId(Option<CloneId>);

impl std::fmt::Display for ParentId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.0 {
            Some(id) => write!(f, "{id}"),
            None => write!(f, "None"),
        }
    }
}

pub fn save_variant_phylogeny(subclones: &SubClones, path2file: &Path) -> anyhow::Result<()> {
    let path2file = path2file.with_extension("csv");
    let parent_ids: Vec<ParentId> = (0..MAX_SUBCLONES)
        .map(|i| ParentId(subclones.get_clone(i as CloneId).unwrap().get_parent_id()))
        .collect();
    debug!("variant phylogeny in {path2file:#?}");
    write2file(&parent_ids, &path2file, None, false)?;
    Ok(())
}

pub fn proliferating_cell(
    subclones: &mut SubClones,
    subclone_id: CloneId,
    rng: &mut impl Rng,
) -> StemCell {
    //! Determine which cells will proliferate by randomly selecting a cell
    //! from the subclone with id `subclone_id`.
    //! This will also remove the random selected cell from the subclone,
    //! hence subclones are borrowed mutably.
    trace!(
        "a cell from clone {:#?} will divide",
        subclones.0[subclone_id as usize]
    );
    subclones.0[subclone_id as usize]
        .random_cell(rng)
        .with_context(|| "found empty subclone")
        .unwrap()
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroU8;

    use crate::tests::LambdaFromNonZeroU8;

    use super::*;
    use quickcheck_macros::quickcheck;
    use rand::{rngs::SmallRng, SeedableRng};

    #[quickcheck]
    fn assign_cell_test(id: CloneId) -> bool {
        let mut neutral_clone = SubClone {
            cells: vec![],
            id,
            parent_id: default_parent_id(id),
        };
        let cell = StemCell::new();
        assert!(neutral_clone.cells.is_empty());

        neutral_clone.assign_cell(cell);
        !neutral_clone.cells.is_empty()
    }

    #[quickcheck]
    fn proliferating_cell_test(seed: u64, mut subclone_id: u8) -> bool {
        let mut rng = SmallRng::seed_from_u64(seed);
        let mut subclones = SubClones::default();
        if subclone_id >= subclones.0.len() as u8 {
            subclone_id = 1;
        }
        let number_of_cells = subclones.0[subclone_id as usize].get_cells().len();
        let tot_cells = subclones.compute_tot_cells();

        let _ = proliferating_cell(&mut subclones, subclone_id as CloneId, &mut rng);

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
            Variants::variant_fractions(&subclones),
            &[1. / MAX_SUBCLONES as f32; MAX_SUBCLONES]
        );
    }

    #[test]
    fn into_subsampled_preserves_parent_id_for_surviving_clones() {
        // Set up: clone 0 (wild-type), clone 7 with parent_id=0, clone 42
        // with parent_id=7. Subsample to a size that keeps cells in clones 0
        // and 7 but is small enough to risk dropping clone 42 — we'll seed
        // the RNG and verify whichever clones survive carry their parent_id.
        let mut rng = SmallRng::seed_from_u64(123);
        let mut subclones = SubClones::new_empty();
        for _ in 0..20 {
            subclones
                .get_mut_clone_unchecked(0)
                .assign_cell(StemCell::new());
        }
        for _ in 0..5 {
            subclones
                .get_mut_clone_unchecked(7)
                .assign_cell(StemCell::new());
        }
        for _ in 0..2 {
            subclones
                .get_mut_clone_unchecked(42)
                .assign_cell(StemCell::new());
        }
        // Pick non-default `parent_id`s so the surviving-vs-drained
        // assertion is distinguishable from the constructor default
        // `Some(0)` that drained slots fall back to.
        subclones.get_mut_clone_unchecked(7).set_parent_id(3);
        subclones.get_mut_clone_unchecked(42).set_parent_id(7);

        let subsampled = subclones.into_subsampled(NonZeroUsize::new(10).unwrap(), &mut rng);

        // Invariant: clones that retain at least one cell keep their
        // `parent_id`; clones that lost every cell fall back to the
        // constructor default (`None` for slot 0, `Some(0)` elsewhere).
        for i in 0..MAX_SUBCLONES {
            let slot = subsampled.get_clone(i as CloneId).unwrap();
            let expected_parent: Option<CloneId> = match i {
                0 => None,
                7 if !slot.is_empty() => Some(3),
                42 if !slot.is_empty() => Some(7),
                _ => Some(0),
            };
            assert_eq!(
                slot.get_parent_id(),
                expected_parent,
                "slot {i} (empty={}) had wrong parent_id",
                slot.is_empty()
            );
        }
        // At least one of the non-trivial clones should still be alive,
        // otherwise the test is vacuous.
        assert!(
            !subsampled.get_clone(7).unwrap().is_empty()
                || !subsampled.get_clone(42).unwrap().is_empty(),
            "subsample dropped both non-neutral clones; test is vacuous"
        );
    }

    #[test]
    fn save_variant_phylogeny_writes_parent_ids() {
        let mut subclones = SubClones::new_empty();
        subclones.get_mut_clone_unchecked(7).set_parent_id(3);
        subclones.get_mut_clone_unchecked(42).set_parent_id(7);

        let path =
            std::env::temp_dir().join(format!("hsc-phylogeny-test-{}.csv", std::process::id()));
        let _ = std::fs::remove_file(&path);

        save_variant_phylogeny(&subclones, &path).unwrap();

        let body = std::fs::read_to_string(&path).unwrap();
        let values: Vec<&str> = body.trim_end_matches(',').split(',').collect();
        assert_eq!(values.len(), MAX_SUBCLONES);
        for (i, v) in values.iter().enumerate() {
            let expected: &str = match i {
                0 => "None",
                7 => "3",
                42 => "7",
                _ => "0",
            };
            assert_eq!(*v, expected, "slot {i} mismatched");
        }

        std::fs::remove_file(&path).unwrap();
    }

    #[quickcheck]
    fn division_no_new_clone(seed: u64, cells_present: NonZeroU8) -> bool {
        let mut rng = SmallRng::seed_from_u64(seed);
        let mut cell = StemCell::new();
        cell.set_last_division_time(1.1).unwrap();

        let cells = vec![cell; cells_present.get() as usize];
        let subclones = SubClones::new(cells, cells_present.get() as usize + 1);

        let old_id = 0;
        let id = next_clone(&subclones, old_id, 0., &mut rng);
        id == old_id
    }

    #[quickcheck]
    fn division_new_clone(seed: u64, cells_present: NonZeroU8) -> bool {
        let mut rng = SmallRng::seed_from_u64(seed);
        let mut cell = StemCell::new();
        cell.set_last_division_time(1.1).unwrap();

        let cells = vec![cell; cells_present.get() as usize];
        let before_assignment = cells.len();
        let subclones = SubClones::new(cells, cells_present.get() as usize + 1);

        next_clone(&subclones, 0, 1., &mut rng);
        subclones.0[0].cell_count() as usize == before_assignment
    }

    #[test]
    #[should_panic]
    fn assign_all_clones_occupied() {
        let mut rng = SmallRng::seed_from_u64(26);

        let subclones = SubClones::default();
        let mut cell = StemCell::new();
        cell.set_last_division_time(1.).unwrap();

        next_clone(&subclones, 0, 1., &mut rng);
    }

    #[quickcheck]
    fn new_distribution_test(
        lambda_division: LambdaFromNonZeroU8,
        lambda_background: LambdaFromNonZeroU8,
    ) -> bool {
        let probs = Probs::new(lambda_background.0, lambda_division.0, 0.01, 0., 10);
        let distrs = Distributions::new(probs);
        distrs.neutral_poisson.eq(&NeutralMutationPoisson::new(
            lambda_division.0,
            lambda_background.0,
        )
        .unwrap())
    }
}
