use crate::stemcell::{NbPoissonMutations, Sfs, StemCell};
use crate::subclone::{CloneId, SubClone};
use crate::{write2file, MAX_SUBCLONES};
use anyhow::Context;
use rand::Rng;
use rand_distr::{Bernoulli, Distribution, Poisson, WeightedIndex};
use sosa::{AdvanceStep, CurrentState, NextReaction};
use std::collections::{HashMap, VecDeque};
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Hash, PartialEq, Eq)]
enum Stats2Save {
    VariantFraction,
    SfsNeutral,
    Sfs,
}

/// Number of cells in subclones.
pub struct Variants {}

impl Variants {
    pub fn variant_counts(subclones: &[SubClone; MAX_SUBCLONES]) -> Vec<u64> {
        //! The total variant count is the number of cells in all subclones.
        //! ```
        //! use hsc::MAX_SUBCLONES;
        //! # use hsc::{stemcell::StemCell, process::{HSCProcess, Variants}};
        //! // create a process with one cell in each `MAX_SUBCLONES` subclones
        //! let mut hsc = HSCProcess::default();
        //!
        //! assert_eq!(Variants::variant_counts(&hsc.subclones), [1; MAX_SUBCLONES]);
        //! ```
        let mut variant_counts = Vec::with_capacity(MAX_SUBCLONES);
        // skip the neutral one
        for subclone in subclones {
            variant_counts.push(subclone.cell_count());
        }
        variant_counts
    }

    pub fn variant_fractions(subclones: &[SubClone; MAX_SUBCLONES], tot_cells: u64) -> Vec<f32> {
        //! The proportion of cells in all subclones.
        //!
        //! ```
        //! use hsc::MAX_SUBCLONES;
        //! # use hsc::{stemcell::StemCell, process::{HSCProcess, Variants}};
        //! // create a process with one cell in each `MAX_SUBCLONES` subclones
        //! let mut hsc = HSCProcess::default();
        //!
        //! assert_eq!(
        //!     Variants::variant_fractions(&hsc.subclones, hsc.compute_tot_cells()),
        //!     [1. / MAX_SUBCLONES as f32; MAX_SUBCLONES]
        //! );
        //! ```
        Variants::variant_counts(subclones)
            .into_iter()
            .map(|frac| (frac as f32 / tot_cells as f32))
            .collect()
    }
}

/// The Moran process saves the state of the agents and simulates new
/// proliferative events at each timestep according to the Gillespie algorithm,
/// see [`HSCProcess::advance_step`].
#[derive(Debug, Clone)]
pub struct HSCProcess {
    /// A collection of clones having a proliferative advantage.
    pub subclones: [SubClone; MAX_SUBCLONES],
    /// The counter for the number of proliferative events.
    counter_divisions: usize,
    id: usize,
    time: f32,
    pub snapshot: VecDeque<f32>,
    pub path2dir: PathBuf,
    pub verbosity: u8,
    distributions: Distributions,
}

impl Default for HSCProcess {
    fn default() -> Self {
        let subclones = std::array::from_fn(|i| {
            let mut subclone = SubClone::new(i, 2);
            subclone.assign_cell(StemCell::new());
            subclone
        });
        HSCProcess::new(
            CellDivisionProbabilities {
                p_asymmetric: 0.5,
                lambda_poisson: 1.,
                p: 0.1,
            },
            subclones,
            vec![0.01, 0.1],
            PathBuf::from("./output"),
            0,
            0.,
            1,
        )
    }
}

#[derive(Debug, Clone)]
pub struct CellDivisionProbabilities {
    pub p_asymmetric: f64,
    /// Rate of neutral mutations per cell-division
    pub lambda_poisson: f32,
    /// Probability of getting a fit mutant upon cell division
    pub p: f64,
}

impl HSCProcess {
    /// A Moran process with wild-type subclone with neutral fitness, to which
    /// all cells will be assigned.
    pub fn new(
        probabilities: CellDivisionProbabilities,
        initial_subclones: [SubClone; MAX_SUBCLONES],
        mut snapshot: Vec<f32>,
        path2dir: PathBuf,
        id: usize,
        time: f32,
        verbosity: u8,
    ) -> HSCProcess {
        let distributions = Distributions::new(
            probabilities.p_asymmetric,
            probabilities.lambda_poisson,
            probabilities.p,
            verbosity,
        );
        snapshot.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let snapshot = VecDeque::from(snapshot);
        let hsc = HSCProcess {
            subclones: initial_subclones,
            distributions,
            counter_divisions: 0,
            id,
            path2dir,
            time,
            snapshot,
            verbosity,
        };
        if verbosity > 1 {
            println!("process created: {:#?}", hsc);
        }
        hsc
    }

    pub fn compute_tot_cells(&self) -> u64 {
        self.subclones
            .iter()
            .map(|subclone| subclone.cell_count())
            .sum()
    }

    pub fn the_only_one_subclone_present(&self) -> Option<usize> {
        //! returns `None` if more than one subclone is present else the clone
        //! id.
        for subclone in self.subclones.iter() {
            if subclone.cell_count() == self.compute_tot_cells() {
                return Some(subclone.id);
            }
        }
        None
    }

    fn assign(&mut self, subclone_id: usize, stem_cell: StemCell, rng: &mut impl Rng) {
        //! The `stem_cell` will be assigned to a new clone
        //! `subclone_id` or not. If that's the case, perform assignment.
        let cell = assign(
            &mut self.subclones[subclone_id],
            stem_cell,
            &self.distributions,
            rng,
        );

        // assign cell to a new random clone if `assign` returned some cell
        if let Some(cell) = cell {
            let mut rnd_clone_id = rng.gen_range(0..self.subclones.len());
            let mut counter = 0;
            // the new random clone cannot have `subclone_id` id and must be
            // empty
            while (rnd_clone_id == subclone_id || !self.subclones[rnd_clone_id].is_empty())
                && counter <= self.subclones.len()
            {
                rnd_clone_id = rng.gen_range(0..self.subclones.len());
                counter += 1;
            }
            assert!(
                counter <= self.subclones.len(),
                "max number of clones reached"
            );
            if self.verbosity > 2 {
                println!(
                    "assgin {:#?} to clone {:#?}",
                    cell, self.subclones[rnd_clone_id]
                );
            }
            self.subclones[rnd_clone_id].assign_cell(cell);
        } else if cell.is_none() && self.verbosity > 2 {
            println!("no new fit variants");
        }
    }

    fn proliferating_cells(&mut self, subclone_id: usize, rng: &mut impl Rng) -> Vec<StemCell> {
        //! Determine which cells will proliferate by randomly selecting a cell
        //! from the subclone with id `subclone_id`.
        if self.verbosity > 2 {
            println!(
                "a cell from clone {:#?} will divide",
                self.subclones[subclone_id]
            );
        }
        let mut proliferating_cells = Vec::with_capacity(2);
        if self.distributions.bern_asymmetric.sample(rng) {
            proliferating_cells.push(
                self.subclones[subclone_id]
                    .random_cell(rng)
                    .with_context(|| "found empty subclone")
                    .unwrap(),
            )
        } else {
            // remove a cell from a random subclone based on the frequencies of
            // the clones at the current state
            let id2remove = if let Some(id) = self.the_only_one_subclone_present() {
                id
            } else {
                let mut variants = Variants::variant_counts(&self.subclones);
                // do not sample the cell that will proliferate
                variants[subclone_id] -= 1;
                WeightedIndex::new(variants).unwrap().sample(rng)
            };
            // remove a cell from the random subclone
            let _ = self.subclones[id2remove]
                .random_cell(rng)
                .with_context(|| "found empty subclone")
                .unwrap();
            if self.verbosity > 2 {
                println!("removing one cell from clone {}", id2remove);
            }
            let cell1 = self.subclones[subclone_id]
                .random_cell(rng)
                .with_context(|| "found empty subclone")
                .unwrap();
            let cell2 = cell1.clone();
            proliferating_cells.push(cell1);
            proliferating_cells.push(cell2);
        }
        if self.verbosity > 2 {
            println!("proliferating cells {:#?}", proliferating_cells)
        }
        proliferating_cells
    }

    fn save_variant_fraction(&self, path2file: &Path) -> anyhow::Result<()> {
        let path2file = path2file.with_extension("csv");
        let total_variant_frac =
            Variants::variant_fractions(&self.subclones, self.compute_tot_cells());
        if self.verbosity > 1 {
            println!("total variant fraction in {:#?}", path2file)
        }
        write2file(&total_variant_frac, &path2file, None, false)?;
        Ok(())
    }

    fn save_sfs(&self, path2file: &Path, neutral: bool, rng: &mut impl Rng) -> anyhow::Result<()> {
        let path2file = path2file.with_extension("json");
        if neutral {
            // subclone 0 is the neutral one
            let cells = self.subclones[0].get_cells();
            // some simulations might no cells in the wild-type clone
            if let Ok(sfs) =
                &Sfs::from_cells(cells, &self.distributions.poisson, self.verbosity, rng)
                    .with_context(|| "cannot construct the sfs neutral")
            {
                let sfs = serde_json::to_string(sfs)
                    .with_context(|| "cannot serialize the neutral sfs")?;
                fs::write(path2file, sfs)
                    .with_context(|| "cannot save the neutral SFS".to_string())?;
            }
        } else {
            let cells: Vec<StemCell> = self
                .subclones
                .iter()
                .flat_map(|subclone| subclone.get_cells().to_vec())
                .collect();
            // some simulations might have only cells in the wild-type clone
            if let Ok(sfs) =
                &Sfs::from_cells(&cells, &self.distributions.poisson, self.verbosity, rng)
            {
                let sfs = serde_json::to_string(sfs).with_context(|| "cannot serialize the sfs")?;
                fs::write(path2file, sfs)
                    .with_context(|| "Cannot save the total SFS".to_string())?;
            }
        }

        Ok(())
    }

    pub fn save(&self, timepoint: usize, rng: &mut impl Rng) -> anyhow::Result<()> {
        //! Save the SFS and SFS neutral.
        let make_path = |tosave: Stats2Save| -> anyhow::Result<PathBuf> {
            let path2file = match tosave {
                Stats2Save::VariantFraction => self
                    .path2dir
                    .join("variant_fraction")
                    .join(timepoint.to_string()),
                Stats2Save::Sfs => self.path2dir.join("sfs").join(timepoint.to_string()),
                Stats2Save::SfsNeutral => self
                    .path2dir
                    .join("sfs_neutral")
                    .join(timepoint.to_string()),
            };
            fs::create_dir_all(&path2file).with_context(|| "Cannot create dir")?;
            Ok(path2file.join(self.id.to_string()))
        };
        let mut paths = HashMap::with_capacity(3);
        paths.insert(
            Stats2Save::VariantFraction,
            make_path(Stats2Save::VariantFraction)?,
        );
        paths.insert(Stats2Save::Sfs, make_path(Stats2Save::Sfs)?);
        paths.insert(Stats2Save::SfsNeutral, make_path(Stats2Save::SfsNeutral)?);

        self.save_variant_fraction(&paths[&Stats2Save::VariantFraction])?;
        self.save_sfs(&paths[&Stats2Save::SfsNeutral], true, rng)?;
        self.save_sfs(&paths[&Stats2Save::Sfs], false, rng)?;

        Ok(())
    }
}

impl AdvanceStep<MAX_SUBCLONES> for HSCProcess {
    /// The id of a subclone [`CloneId`].
    type Reaction = CloneId;

    fn advance_step(&mut self, reaction: NextReaction<Self::Reaction>, rng: &mut impl Rng) {
        //! Update the process by simulating the next proliferative event
        //! according to the next `reaction` determined by the Gillespie
        //! algorithm.
        //!
        //! The proliferation step is implemented as following:
        //!
        //! 1. select the cell that will proliferate next from the clone
        //! with id `reaction` determined by the Gillespie algorithm
        //!
        //! 2. if the next proliferation is performed symmetricaly, clone the
        //! proliferating cell, else continue to step 3
        //!
        //! 3. for all proliferating cells (1 cell in case of a asymmetric
        //! division, 2 cells in the case of a symmetric division):
        //!     * mutate genome by storing the division id, see
        //!     [`crate::stemcell::Sfs`]
        //!
        //!     * assign to new subclone with a probability determined by the
        //!     rate of mutations conferring a proliferative advantage
        // The Gillespie sampler samples the clone that will proliferate next,
        // that is the clone with id `reaction.event`.
        // Pick random proliferating cells from this clone. **Note that this
        // removes the cells from the clone with id `reaction.event`.**
        let stem_cells = self.proliferating_cells(reaction.event, rng);

        if self.verbosity > 2 {
            println!("{:#?} cells are dividing", stem_cells);
        }

        for mut stem_cell in stem_cells.into_iter() {
            // perform division
            self.counter_divisions += 1;
            // mutate cell
            stem_cell.record_proliferation_event(self.counter_divisions);
            // assign cell to clone. This can have two outcomes based on a
            // Bernouilli trial see `assign` and `self.assign`:
            // 1. the stem cell is re-assigned to the old clone with id
            // `reaction.event` (remember that `self.proliferating_cells` has
            // removed the cells from the clone)
            // 2. the stem cell is assigned to a new clone with id different
            // from `reaction.event`.
            self.assign(reaction.event, stem_cell, rng);
        }

        // take snapshot
        self.time += reaction.time;
        if self.verbosity > 1 {
            println!("time: {}", self.time);
        }
        if self.verbosity > 1 {
            println!(
                "{:#?} timepoints to save, time now is {}",
                self.snapshot, self.time
            );
        }
        if let Some(&time) = self.snapshot.front() {
            if self.time >= time {
                if self.verbosity > 0 {
                    println!(
                        "saving variant fraction for {:#?} at time {}",
                        time, self.time
                    );
                }
                self.save(self.snapshot.len(), rng)
                    .expect("cannot save snapshot");
                self.snapshot.pop_front();
            }
        }
    }

    fn update_state(&self, state: &mut CurrentState<MAX_SUBCLONES>) {
        for i in 0..state.population.len() {
            state.population[i] = self.subclones[i].cell_count();
        }
    }
}

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
        mutations as u8
    }
}

/// Distribution probabilities used to model acquisition of neutral and
/// conferring-fitness mutations upon cell division.
#[derive(Debug, Clone)]
pub struct Distributions {
    poisson: NeutralMutationPoisson,
    bern: Bernoulli,
    bern_asymmetric: Bernoulli,
}

impl Distributions {
    pub fn new(p_asymmetric: f64, lambda_poisson: f32, p_fitness: f64, verbosity: u8) -> Self {
        if verbosity > 1 {
            println!("creating distributions with parameters asymmetric: {}, neutral_lambda: {}, p_fitness: {}", p_asymmetric, lambda_poisson, p_fitness);
        }
        Self {
            poisson: NeutralMutationPoisson(
                Poisson::new(lambda_poisson).expect("Invalid lambda found: lambda <= 0 or nan"),
            ),
            bern: Bernoulli::new(p_fitness).expect("Invalid p: p<0 or p>1"),
            bern_asymmetric: Bernoulli::new(p_asymmetric).expect("Invalid p: p<0 or p>1"),
        }
    }

    pub fn acquire_p_mutation(&self, rng: &mut impl Rng) -> bool {
        self.bern.sample(rng)
    }
}

fn assign(
    subclone: &mut SubClone,
    cell: StemCell,
    distr: &Distributions,
    rng: &mut impl Rng,
) -> Option<StemCell> {
    //! Check if `cell` will be assigned to `subclone`, according to a
    //! Bernouilli trial with probability `p` (see [`Distributions::new`]).
    //! Assign cell to `subclone` if no fit variant has been generated.
    //!
    //! ## Returns
    //! If the cell gets one fitness advantage mutation, then the function
    //! returns the cell, otherwise it returns None.
    if distr.acquire_p_mutation(rng) {
        return Some(cell);
    }
    subclone.assign_cell(cell);
    None
}

#[cfg(test)]
mod tests {
    use crate::tests::LambdaFromNonZeroU8;

    use super::*;
    use quickcheck_macros::quickcheck;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;
    use rand_distr::Uniform;

    #[should_panic]
    #[test]
    fn new_distribution_wrong_lambda_test() {
        Distributions::new(1., f32::NAN, 0.9, 0);
    }

    #[should_panic]
    #[test]
    fn new_distribution_wrong_lambda_neg_test() {
        Distributions::new(1., -1.0, 0.9, 0);
    }

    #[should_panic]
    #[test]
    fn new_distribution_wrong_p_test() {
        Distributions::new(1., 1., f64::NAN, 0);
    }

    #[should_panic]
    #[test]
    fn new_distribution_wrong_p_inf_test() {
        Distributions::new(1., 1., f64::INFINITY, 0);
    }

    #[should_panic]
    #[test]
    fn new_distribution_wrong_p_neg_test() {
        Distributions::new(1., 1.0, -0.9, 0);
    }

    #[quickcheck]
    fn test_nb_neutral_mutations(lambda_poisson: LambdaFromNonZeroU8) {
        let mut rng = ChaCha8Rng::seed_from_u64(26);
        NeutralMutationPoisson(Poisson::new(lambda_poisson.0).unwrap())
            .nb_neutral_mutations(&mut rng);
    }

    #[quickcheck]
    fn division_no_new_clone(lambda: LambdaFromNonZeroU8, seed: u64) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let distr = Distributions::new(1., lambda.0, 0f64, 0);

        let mut neutral_clone = SubClone::new(1, 2);
        let cell = StemCell::new();

        assign(&mut neutral_clone, cell, &distr, &mut rng).is_none()
    }

    #[quickcheck]
    fn division_new_clone(lambda: LambdaFromNonZeroU8, seed: u64) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let distr = Distributions::new(1., lambda.0, 1f64, 0);

        let mut neutral_clone = SubClone::new(1, 2);
        let cell = StemCell::new();

        assign(&mut neutral_clone, cell, &distr, &mut rng).is_some()
    }

    #[test]
    fn variant_fraction() {
        let hsc = HSCProcess::default();
        assert_eq!(
            Variants::variant_fractions(&hsc.subclones, hsc.compute_tot_cells()),
            [1. / MAX_SUBCLONES as f32; MAX_SUBCLONES]
        );
    }

    #[quickcheck]
    fn proliferating_cells_test(seed: u64, mut subclone_id: u8) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let mut hsc = HSCProcess::default();
        if subclone_id >= hsc.subclones.len() as u8 {
            subclone_id = rng.sample(Uniform::new(0, hsc.subclones.len())) as u8;
        }
        let tot_cells = hsc.compute_tot_cells();

        let proliferating_cells = hsc.proliferating_cells(subclone_id as usize, &mut rng);

        let mut subclone_id_has_lost_cell = true;
        if proliferating_cells.len() == 2 {
            subclone_id_has_lost_cell = hsc.subclones[subclone_id as usize].cell_count() == 0;
        }
        let subclone_has_lost_cell =
            Variants::variant_counts(&hsc.subclones).iter().sum::<u64>() < tot_cells;
        subclone_has_lost_cell
            && proliferating_cells.iter().all(|cell| !cell.has_mutations())
            && subclone_id_has_lost_cell
    }

    #[quickcheck]
    fn poisson_neutral_mutations(seed: u64) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let poisson = Poisson::new(0.0001).unwrap();
        NeutralMutationPoisson(poisson).nb_neutral_mutations(&mut rng) < 2
    }

    #[test]
    #[should_panic]
    fn assign_all_clones_occupied() {
        let mut process = HSCProcess::default();
        process.distributions.bern = Bernoulli::new(1.).unwrap();
        let cell = StemCell::new();
        let mut rng = ChaCha8Rng::seed_from_u64(26);
        process.assign(0, cell, &mut rng);
    }
}
