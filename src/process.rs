use crate::genotype::NbPoissonMutations;
use crate::stemcell::{save_cells, StemCell};
use crate::subclone::{assign, CloneId, SubClone};
use crate::{write2file, MAX_SUBCLONES};
use anyhow::Context;
use rand::Rng;
use rand_distr::{Bernoulli, Distribution, Poisson, WeightedIndex};
use sosa::{AdvanceStep, CurrentState, NextReaction};
use std::collections::{HashMap, VecDeque};
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Hash, PartialEq, Eq)]
pub enum Stats2Save {
    VariantFraction,
    Sfs,
    Genotypes,
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
    pub distributions: Distributions,
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
                p_asymmetric: 0.,
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

    fn proliferating_cell(&mut self, subclone_id: usize, rng: &mut impl Rng) -> StemCell {
        //! Determine which cells will proliferate by randomly selecting a cell
        //! from the subclone with id `subclone_id`.
        self.counter_divisions += 1;
        if self.verbosity > 2 {
            println!(
                "a cell from clone {:#?} will divide",
                self.subclones[subclone_id]
            );
        }
        self.subclones[subclone_id]
            .random_cell(rng)
            .with_context(|| "found empty subclone")
            .unwrap()
    }

    fn keep_const_population_upon_symmetric_division(
        &mut self,
        subclone_id: usize,
        rng: &mut impl Rng,
    ) {
        //! If an symmetric division is performed, need to remove a random cell
        //! from the population.
        //!
        //! We proceed as following:
        //!     1. check if there is only one clone in the population, then
        //!     remove a cell from this clone
        //!     2. else, compute the variant counts (removing one cell to avoid
        //!     sampling the cell that will proliferate)
        //!     3. and sample from any clone based on the weights defined by
        //!     the variant counts
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

    pub fn make_path(&self, tosave: Stats2Save, timepoint: usize) -> anyhow::Result<PathBuf> {
        let path2file = match tosave {
            Stats2Save::VariantFraction => self.path2dir.join("variant_fraction"),
            Stats2Save::Genotypes => self.path2dir.join("genotypes"),
            Stats2Save::Sfs => self.path2dir.join("burden"),
        };
        let path2file = path2file.join(timepoint.to_string());
        fs::create_dir_all(&path2file).with_context(|| "Cannot create dir")?;
        Ok(path2file.join(self.id.to_string()))
    }

    pub fn save(&mut self, timepoint: usize) -> anyhow::Result<()> {
        //! Save the cells with their genotypes. Save also the variant fraction
        //! that is the fraction of subclones present in the population.
        let mut paths = HashMap::with_capacity(3);
        paths.insert(
            Stats2Save::VariantFraction,
            self.make_path(Stats2Save::VariantFraction, timepoint)?,
        );
        self.save_variant_fraction(&paths[&Stats2Save::VariantFraction])?;

        paths.insert(
            Stats2Save::Genotypes,
            self.make_path(Stats2Save::Genotypes, timepoint)?,
        );
        let cells: Vec<StemCell> = self
            .subclones
            .iter()
            .flat_map(|subclone| subclone.get_cells().to_vec())
            .collect();
        save_cells(
            &cells,
            &paths[&Stats2Save::Genotypes].with_extension("json"),
        )?;

        if self.verbosity > 0 {
            println!(
                "saving measurements after {} mutational events",
                self.counter_divisions
            );
        }

        Ok(())
    }
}

impl AdvanceStep<MAX_SUBCLONES> for HSCProcess {
    /// The id of a subclone [`CloneId`].
    type Reaction = CloneId;

    fn advance_step(&mut self, reaction: NextReaction<Self::Reaction>, rng: &mut impl Rng) {
        //! Update the process by simulating the next proliferative event
        //! according to the next `reaction` determined by the [Gillespie
        //! algorithm](`sosa`).
        //!
        //! The proliferation step is implemented as following:
        //!
        //! 1. select the cell that will proliferate next from the clone
        //! with id `reaction` determined by the Gillespie algorithm
        //!
        //! 2. if the next proliferation is performed symmetricaly, clone the
        //! proliferating cell and assign it to the subclone, else continue to
        //! step 3
        //!
        //! 3. for the proliferating cell only (not the clone from 2):
        //!     * mutate genome by storing the division id, see
        //!     [`crate::genotype::StatisticsMutations`]
        //!
        //!     * assign to new subclone with a probability determined by the
        //!     rate of mutations conferring a proliferative advantage
        // The Gillespie sampler samples the clone that will proliferate next,
        // that is the clone with id `reaction.event`.
        // Pick random proliferating cells from this clone. **Note that this
        // removes the cells from the clone with id `reaction.event`.**
        let mut stem_cell = self.proliferating_cell(reaction.event, rng);

        if self.distributions.can_only_be_symmetric()
            || !self.distributions.bern_asymmetric.sample(rng)
        {
            if self.verbosity > 2 {
                println!("keeping the cell population constant");
            }
            // remove a cell from the population
            self.keep_const_population_upon_symmetric_division(reaction.event, rng);
            self.subclones[reaction.event].assign_cell(stem_cell.clone());
        }

        if self.verbosity > 2 {
            println!("assigning mutations to cell {:#?}", stem_cell)
        }
        stem_cell.record_division(self.counter_divisions);

        // assign cell to clone. This can have two outcomes based on a
        // Bernouilli trial see `assign` and `self.assign`:
        // 1. the stem cell is re-assigned to the old clone with id
        // `reaction.event` (remember that `self.proliferating_cells` has
        // removed the cells from the clone)
        // 2. the stem cell is assigned to a new clone with id different
        // from `reaction.event`.
        self.assign(reaction.event, stem_cell, rng);

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
                        "saving variant fraction for {:#?} at time {} for timepoint {}",
                        time,
                        self.time,
                        self.snapshot.len()
                    );
                }
                self.save(self.snapshot.len())
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
        mutations as u16
    }
}

/// Distribution probabilities used to model acquisition of neutral and
/// conferring-fitness mutations upon cell division.
#[derive(Debug, Clone)]
pub struct Distributions {
    pub poisson: NeutralMutationPoisson,
    bern: Bernoulli,
    bern_asymmetric: Bernoulli,
    no_asymmetric_division: bool,
}

impl Distributions {
    pub fn new(p_asymmetric: f64, lambda_poisson: f32, p_fitness: f64, verbosity: u8) -> Self {
        if verbosity > 1 {
            println!("creating distributions with parameters asymmetric: {}, neutral_lambda: {}, p_fitness: {}", p_asymmetric, lambda_poisson, p_fitness);
        }
        let no_asymmetric_division = (p_asymmetric - 0.).abs() < f64::EPSILON;
        Self {
            poisson: NeutralMutationPoisson(
                Poisson::new(lambda_poisson).expect("Invalid lambda found: lambda <= 0 or nan"),
            ),
            bern: Bernoulli::new(p_fitness).expect("Invalid p: p<0 or p>1"),
            bern_asymmetric: Bernoulli::new(p_asymmetric).expect("Invalid p: p<0 or p>1"),
            no_asymmetric_division,
        }
    }

    pub fn acquire_p_mutation(&self, rng: &mut impl Rng) -> bool {
        self.bern.sample(rng)
    }

    pub fn can_only_be_symmetric(&self) -> bool {
        self.no_asymmetric_division
    }
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
    fn proliferating_cell_test(seed: u64, mut subclone_id: u8) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let mut hsc = HSCProcess::default();
        if subclone_id >= hsc.subclones.len() as u8 {
            subclone_id = rng.sample(Uniform::new(0, hsc.subclones.len())) as u8;
        }
        let number_of_cells = hsc.subclones[subclone_id as usize].get_cells().len();
        let tot_cells = hsc.compute_tot_cells();

        let _ = hsc.proliferating_cell(subclone_id as usize, &mut rng);

        let subclones_have_lost_cell =
            Variants::variant_counts(&hsc.subclones).iter().sum::<u64>() < tot_cells;
        let subclone_has_lost_cell =
            hsc.subclones[subclone_id as usize].get_cells().len() == number_of_cells - 1;
        subclones_have_lost_cell && subclone_has_lost_cell
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
