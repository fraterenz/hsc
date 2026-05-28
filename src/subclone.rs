use crate::Probs;
use crate::genotype::NeutralMutationPoisson;
use crate::{MAX_SUBCLONES, stemcell::StemCell, write2file};
use anyhow::{Context, ensure};
use log::{debug, trace};
use rand::Rng;
use rand::RngExt;
use rand::seq::{IndexedRandom, IteratorRandom};
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

    /// Draw a Gillespie rate for a non-wild-type subclone under this fitness
    /// model. `b0` is the wild-type birth rate; the returned value is
    /// `b0 * (1 + s)` where `s` is fixed, sampled, or zero (for `Neutral`)
    /// depending on the variant.
    pub fn sample_rate(&self, b0: f32, rng: &mut impl Rng) -> f32 {
        match self {
            Fitness::Neutral => b0,
            Fitness::Fixed { s } => b0 * (1. + s),
            &Fitness::GammaSampled { shape, scale } => {
                let gamma = Gamma::new(shape, scale).unwrap();
                b0 * (1. + gamma.sample(rng))
            }
        }
    }
}

/// Number of fit-mutation hits accumulated by a subclone lineage.
///
/// Saturates at [`HitCount::Fourth`]: any further hit beyond the 4th maps back
/// onto the 4th tier. Used by the `--multihits` feature to pick a fitness tier
/// when a new subclone is born.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum HitCount {
    #[default]
    WildType,
    First,
    Second,
    Third,
    Fourth,
}

impl HitCount {
    /// The following `HitCount` with a saturating successor when the hit reaches
    /// `Fourth`: `WildType -> First -> Second -> Third -> Fourth -> Fourth`.
    pub fn child(self) -> HitCount {
        match self {
            HitCount::WildType => HitCount::First,
            HitCount::First => HitCount::Second,
            HitCount::Second => HitCount::Third,
            HitCount::Third => HitCount::Fourth,
            HitCount::Fourth => HitCount::Fourth,
        }
    }
}

/// Id of the [`SubClone`]s.
/// The id of zero indicates the wild-type clone
pub type CloneId = u16;

/// Display wrapper for a subclone's optional parent.
///
/// `Some(p)` formats as the integer `p`, `None` formats as the empty string
/// (read as `NaN` by pandas and other CSV readers). The `:.6` precision
/// specifier used by [`crate::write2file`] is silently ignored.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ParentId(pub(crate) Option<CloneId>);

impl std::fmt::Display for ParentId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.0 {
            Some(p) => write!(f, "{p}"),
            None => Ok(()),
        }
    }
}

#[derive(Debug, Clone, Default)]
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
    parent_id: Option<CloneId>,
    hit_count: HitCount,
}

impl Iterator for SubClone {
    type Item = StemCell;

    fn next(&mut self) -> Option<Self::Item> {
        self.cells.pop()
    }
}

impl SubClone {
    pub fn new(id: CloneId, cell_capacity: usize, parent_id: Option<CloneId>) -> SubClone {
        SubClone {
            cells: Vec::with_capacity(cell_capacity),
            id,
            parent_id,
            hit_count: HitCount::WildType,
        }
    }

    pub fn with_capacity(id: CloneId, capacity: usize) -> SubClone {
        SubClone {
            cells: Vec::with_capacity(capacity),
            id,
            parent_id: None,
            hit_count: HitCount::WildType,
        }
    }

    pub fn empty_with_id(id: CloneId) -> SubClone {
        SubClone {
            cells: vec![],
            id,
            parent_id: None,
            hit_count: HitCount::WildType,
        }
    }

    pub fn parent_id(&self) -> Option<CloneId> {
        self.parent_id
    }

    pub(crate) fn set_parent_id(&mut self, parent: CloneId) {
        self.parent_id = Some(parent);
    }

    pub(crate) fn hit_count(&self) -> HitCount {
        self.hit_count
    }

    pub(crate) fn set_hit_count(&mut self, hit_count: HitCount) {
        self.hit_count = hit_count;
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
        //! A non-wild type subclone is empty is it's parent id is None.
        //! If the subclone is the wild-type, it's empty if it has no cells.
        if self.id == 0 {
            self.get_cells().is_empty()
        } else {
            self.parent_id.is_none()
        }
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
        // unrwap is fine bc subclones has fixed size at comp time `MAX_SUBCLONES`
        let mut rnd_clone_id = subclones.0.as_slice().choose(rng).unwrap().id;
        let mut counter = 0;
        // the new random clone cannot have `subclone_id` id and must be empty
        while (rnd_clone_id == old_subclone_id
            || !subclones.get_clone(rnd_clone_id).unwrap().is_empty())
            && counter <= MAX_SUBCLONES
        {
            rnd_clone_id = subclones.0.as_slice().choose(rng).unwrap().id;
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
pub struct SubClones([SubClone; MAX_SUBCLONES]);

impl SubClones {
    pub fn new(cells: Vec<StemCell>, capacity: usize) -> Self {
        //! Returns all the newly initiated subclones by assigning all `cells`
        //! to the neutral clone.
        //!
        //! All clones have the same `capacity`.
        // initial state
        let mut subclones: [SubClone; MAX_SUBCLONES] =
            std::array::from_fn(|i| SubClone::new(i as CloneId, capacity, None));
        let tot_cells = cells.len();
        debug!("assigning {tot_cells} cells to the wild-type clone");
        for cell in cells {
            subclones[0].assign_cell(cell);
        }
        assert_eq!(subclones[0].cells.len(), tot_cells);

        Self(subclones)
    }

    pub fn new_empty() -> Self {
        SubClones(std::array::from_fn(|i| {
            SubClone::empty_with_id(i as CloneId)
        }))
    }

    pub fn with_capacity(capacity: usize) -> Self {
        SubClones(std::array::from_fn(|id| {
            SubClone::with_capacity(id as CloneId, capacity)
        }))
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
        let mut subclones = Self::new_empty();

        for (cell, clone_id) in self
            .0
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
            .sample(rng, nb_cells.get())
        {
            subclones.0[clone_id as usize].cells.push(cell);
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
            .sample(rng, nb_cells)
    }

    pub fn gillespie_rates(
        &self,
        fitness: &Fitness,
        b0: f32,
        rng: &mut impl Rng,
    ) -> ReactionRates<MAX_SUBCLONES> {
        //! Create the Gillespie reaction rates according to the `fitness`
        //! model with `b0` being the proliferative rate of the wild-type clone.
        ReactionRates(core::array::from_fn(|i| {
            if i == 0 {
                b0
            } else {
                fitness.sample_rate(b0, rng)
            }
        }))
    }
}

impl Default for SubClones {
    fn default() -> Self {
        //! Create new subclones each having one cell (this creates also the cells).
        let subclones = std::array::from_fn(|i| {
            let parent_id = if i == 0 { None } else { Some(0) };
            let mut subclone = SubClone::new(i as CloneId, 2, parent_id);
            subclone.assign_cell(StemCell::new());
            subclone
        });
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

pub(crate) fn save_variant_fraction(subclones: &SubClones, path2file: &Path) -> anyhow::Result<()> {
    let path2file = path2file.with_extension("csv");
    let total_variant_frac = Variants::variant_fractions(subclones);
    debug!("total variant fraction in {path2file:#?}");
    write2file(&total_variant_frac, &path2file, None, false)?;
    Ok(())
}

pub(crate) fn save_phylogeny(subclones: &SubClones, path2file: &Path) -> anyhow::Result<()> {
    //! Save the per-clone parent-id row for the current `subclones` state.
    //!
    //! Emits a dense row of `MAX_SUBCLONES` comma-separated entries; clones
    //! without a parent (the wild type and any never-born slot) appear as
    //! empty fields.
    let path2file = path2file.with_extension("csv");
    let row: Vec<ParentId> = (0..MAX_SUBCLONES as CloneId)
        .map(|id| ParentId(subclones.get_clone(id).unwrap().parent_id()))
        .collect();
    debug!("phylogeny in {path2file:#?}");
    write2file(&row, &path2file, None, false)?;
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
    use rand::{SeedableRng, rngs::SmallRng};

    #[test]
    fn assign_cell_test() {
        let mut neutral_clone = SubClone::default();
        let cell = StemCell::new();
        assert!(neutral_clone.cells.is_empty());

        neutral_clone.assign_cell(cell);
        assert!(!neutral_clone.cells.is_empty());
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
        next_clone(&SubClones::default(), 0, 1., &mut rand::rng());
    }

    #[test]
    fn parent_id_defaults_to_none() {
        assert!(SubClone::default().parent_id().is_none());
        assert!(SubClone::empty_with_id(13).parent_id().is_none());

        for clone in SubClones::new(vec![StemCell::new()], 4).0.iter() {
            assert!(clone.parent_id().is_none());
        }
        for clone in SubClones::new_empty().0.iter() {
            assert!(clone.parent_id().is_none());
        }
        for clone in SubClones::with_capacity(4).0.iter() {
            assert!(clone.parent_id().is_none());
        }
        for clone in SubClones::default().0.iter() {
            if clone.id == 0 {
                assert!(clone.parent_id().is_none());
            } else {
                assert!(clone.parent_id().unwrap() == 0);
            }
        }
    }

    #[test]
    fn sample_rate_neutral_returns_b0() {
        let mut rng = SmallRng::seed_from_u64(7);
        let b0 = 1.7_f32;
        for _ in 0..32 {
            assert_eq!(Fitness::Neutral.sample_rate(b0, &mut rng), b0);
        }
    }

    #[test]
    fn sample_rate_fixed_returns_b0_times_one_plus_s() {
        let mut rng = SmallRng::seed_from_u64(7);
        let b0 = 1.7_f32;
        let s = 0.3_f32;
        for _ in 0..32 {
            assert_eq!(
                Fitness::Fixed { s }.sample_rate(b0, &mut rng),
                b0 * (1.0 + s),
            );
        }
    }

    #[test]
    fn sample_rate_gamma_returns_b0_times_one_plus_positive_draw() {
        let mut rng = SmallRng::seed_from_u64(7);
        let b0 = 1.7_f32;
        let fitness = Fitness::GammaSampled {
            shape: 2.0,
            scale: 0.5,
        };
        for _ in 0..32 {
            let rate = fitness.sample_rate(b0, &mut rng);
            // Gamma(shape, scale) is supported on (0, +inf), so the rate is
            // strictly greater than b0.
            assert!(
                rate > b0,
                "expected rate > b0={b0}, got {rate} (b0*(1+gamma_sample))",
            );
            assert!(rate.is_finite());
        }
    }

    #[test]
    fn hit_count_child_saturates() {
        assert_eq!(HitCount::WildType.child(), HitCount::First);
        assert_eq!(HitCount::First.child(), HitCount::Second);
        assert_eq!(HitCount::Second.child(), HitCount::Third);
        assert_eq!(HitCount::Third.child(), HitCount::Fourth);
        assert_eq!(HitCount::Fourth.child(), HitCount::Fourth);
    }

    #[test]
    fn hit_count_defaults_to_wild_type() {
        assert_eq!(HitCount::default(), HitCount::WildType);
        assert_eq!(SubClone::default().hit_count(), HitCount::WildType);
        assert_eq!(SubClone::empty_with_id(13).hit_count(), HitCount::WildType);
        assert_eq!(
            SubClone::with_capacity(4, 2).hit_count(),
            HitCount::WildType
        );
        assert_eq!(SubClone::new(5, 2, Some(0)).hit_count(), HitCount::WildType);

        for clone in SubClones::new(vec![StemCell::new()], 4).0.iter() {
            assert_eq!(clone.hit_count(), HitCount::WildType);
        }
        for clone in SubClones::new_empty().0.iter() {
            assert_eq!(clone.hit_count(), HitCount::WildType);
        }
        for clone in SubClones::with_capacity(4).0.iter() {
            assert_eq!(clone.hit_count(), HitCount::WildType);
        }
        for clone in SubClones::default().0.iter() {
            assert_eq!(clone.hit_count(), HitCount::WildType);
        }
    }

    #[test]
    fn hit_count_resets_on_subsample_and_from() {
        // Build a SubClones, set a non-default hit_count, then run the
        // reconstruction helpers and confirm hit_count resets to WildType.
        let mut subclones = SubClones::new(vec![StemCell::new(); 4], 4);
        subclones
            .get_mut_clone_unchecked(2)
            .set_hit_count(HitCount::Third);
        assert_eq!(subclones.get_clone(2).unwrap().hit_count(), HitCount::Third);

        let mut rng = SmallRng::seed_from_u64(42);
        let resampled = subclones
            .clone()
            .into_subsampled(NonZeroUsize::new(2).unwrap(), &mut rng);
        for clone in resampled.0.iter() {
            assert_eq!(clone.hit_count(), HitCount::WildType);
        }

        let cells: Vec<(StemCell, CloneId)> = subclones
            .0
            .iter()
            .flat_map(|s| s.get_cells().iter().map(move |c| (c.clone(), s.id)))
            .collect();
        let rebuilt = SubClones::from(cells);
        for clone in rebuilt.0.iter() {
            assert_eq!(clone.hit_count(), HitCount::WildType);
        }
    }

    #[test]
    fn parent_id_resets_on_subsample_and_from() {
        // Build a SubClones, set a parent_id, then run the reconstruction
        // helpers and confirm the parent state is reset.
        let mut subclones = SubClones::new(vec![StemCell::new(); 4], 4);
        subclones.get_mut_clone_unchecked(2).set_parent_id(0);
        assert_eq!(subclones.get_clone(2).unwrap().parent_id(), Some(0));

        let mut rng = SmallRng::seed_from_u64(42);
        let resampled = subclones
            .clone()
            .into_subsampled(NonZeroUsize::new(2).unwrap(), &mut rng);
        for clone in resampled.0.iter() {
            assert!(clone.parent_id().is_none());
        }

        let cells: Vec<(StemCell, CloneId)> = subclones
            .0
            .iter()
            .flat_map(|s| s.get_cells().iter().map(move |c| (c.clone(), s.id)))
            .collect();
        let rebuilt = SubClones::from(cells);
        for clone in rebuilt.0.iter() {
            assert!(clone.parent_id().is_none());
        }
    }

    #[test]
    fn parent_id_display_format() {
        // None -> empty string; Some -> integer; the `:.6` precision specifier
        // used by write2file is silently ignored for our Display impl.
        assert_eq!(format!("{}", ParentId(None)), "");
        assert_eq!(format!("{}", ParentId(Some(42))), "42");
        assert_eq!(format!("{:.6}", ParentId(None)), "");
        assert_eq!(format!("{:.6}", ParentId(Some(42))), "42");
    }

    fn unique_tmpdir(name: &str) -> std::path::PathBuf {
        // Per-test unique directory under the system temp dir, cleaned at the
        // start of the test (best-effort) so reruns start fresh.
        let path = std::env::temp_dir().join(format!(
            "hsc-test-{name}-{}-{:?}",
            std::process::id(),
            std::thread::current().id(),
        ));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).unwrap();
        path
    }

    #[test]
    fn save_phylogeny_writes_dense_row() {
        let dir = unique_tmpdir("save_phylogeny");
        let path = dir.join("row.csv");
        let mut subclones = SubClones::new(vec![StemCell::new()], 4);
        subclones.get_mut_clone_unchecked(2).set_parent_id(0);
        subclones.get_mut_clone_unchecked(7).set_parent_id(2);

        save_phylogeny(&subclones, &path).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        let fields: Vec<&str> = content
            .strip_suffix(',')
            .unwrap_or(&content)
            .split(',')
            .collect();
        assert_eq!(fields.len(), MAX_SUBCLONES);
        assert_eq!(fields[0], ""); // wild-type root
        assert_eq!(fields[1], "");
        assert_eq!(fields[2], "0");
        assert_eq!(fields[7], "2");
        // every other slot stayed at None -> empty
        for (i, f) in fields.iter().enumerate() {
            if i != 2 && i != 7 {
                assert!(f.is_empty(), "slot {i} expected empty, got {f:?}");
            }
        }

        let _ = std::fs::remove_dir_all(&dir);
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
