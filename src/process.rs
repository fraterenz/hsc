use crate::neutral::{Genotype, NbPoissonNeutralMutations, StemCell};
use crate::sfs::{sfs, CloneId, SubClone};
use crate::{write2file, MAX_SUBCLONES};
use anyhow::{ensure, Context};
use rand::Rng;
use rand_distr::{Bernoulli, Distribution, Poisson, WeightedIndex};
use sosa::{AdvanceStep, CurrentState, NextReaction};
use std::fs;
use std::path::PathBuf;

pub struct Variants {}

impl Variants {
    pub fn variant_counts(subclones: &[SubClone; MAX_SUBCLONES]) -> Vec<u64> {
        //! The total variant count is the number of cells in all subclones.
        //! ```
        //! use hsc::MAX_SUBCLONES;
        //! # use hsc::{neutral::StemCell, process::{HSCProcess, Variants}};
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
        //! # use hsc::{neutral::StemCell, process::{HSCProcess, Variants}};
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

#[derive(Debug, Clone)]
pub struct HSCProcess {
    /// A collection of clones having a proliferative advantage.
    pub subclones: [SubClone; MAX_SUBCLONES],
    /// The counter for the number of proliferative events.
    pub counter_divisions: usize,
    pub id: usize,
    pub time: Vec<f32>,
    pub snapshot: Option<Vec<usize>>,
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
            None,
            PathBuf::from("./output"),
            0,
            vec![0.],
            1,
        )
    }
}

pub struct CellDivisionProbabilities {
    pub p_asymmetric: f64,
    /// Rate of
    pub lambda_poisson: f32,
    /// Probability of getting a fit mutant upon cell division
    pub p: f64,
}

impl HSCProcess {
    /// A Moran process with wild-type subclone with neutral fitness, to which
    /// all cells will be assigned.
    ///
    /// # Arguments
    /// * `p_asymmetric` - The parameter of a Bernouilli trial modelling the
    /// probability of a cell to undergo asymmetric division (1 - p symmetric).
    /// * `lambda_poisson` - The parameter of the Poisson distribution
    /// modelling the acquisition of neutral mutations (the average number of
    /// neutral mutations acquired by cells upon cell division).
    /// * `p` - The parameter of a Bernouilli trial modelling the probability
    /// of acquiring a mutation conferring a positive advantage to the cell.
    /// * `initial_subclones` - The collection of subclones used to start the simulation
    ///
    /// # Panics
    /// When p<0 or p>1 and when lambda < 0 or nan.
    /// When `p_asymmetric` is < 0 or > 1
    pub fn new(
        probabilities: CellDivisionProbabilities,
        initial_subclones: [SubClone; MAX_SUBCLONES],
        snapshot: Option<Vec<usize>>,
        path2dir: PathBuf,
        id: usize,
        time: Vec<f32>,
        verbosity: u8,
    ) -> HSCProcess {
        let distributions = Distributions::new(
            probabilities.p_asymmetric,
            probabilities.lambda_poisson,
            probabilities.p,
            verbosity,
        );
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

    fn mutate_and_assign(
        &mut self,
        subclone_id: usize,
        stem_cells: Vec<StemCell>,
        rng: &mut impl Rng,
    ) {
        //! Mutate proliferating cells and assign them to new subclones if needed
        for stem_cell in stem_cells.into_iter() {
            self.counter_divisions += 1;
            let cell = division(
                &mut self.subclones[subclone_id],
                stem_cell,
                &self.distributions,
                self.counter_divisions,
                rng,
            );
            let mut rnd_clone_id = rng.gen_range(0..self.subclones.len());
            let mut counter = 0;
            while rnd_clone_id == subclone_id && counter <= self.subclones.len() {
                rnd_clone_id = rng.gen_range(0..self.subclones.len());
                counter += 1;
            }
            assert!(
                counter <= self.subclones.len(),
                "max number of clones reached"
            );
            if let Some(cell) = cell {
                if self.verbosity > 1 {
                    println!(
                        "division generated {:#?} cells to assign to clone {:#?}",
                        cell, self.subclones[rnd_clone_id]
                    );
                }
                self.subclones[rnd_clone_id].assign_cell(cell);
            } else if cell.is_none() && self.verbosity > 1 {
                println!("division resulted in no new fit variants");
            }
        }
    }

    fn proliferating_cells(&mut self, subclone_id: usize, rng: &mut impl Rng) -> Vec<StemCell> {
        if self.distributions.bern_asymmetric.sample(rng) {
            vec![asymmetric_division(&mut self.subclones[subclone_id], rng).unwrap()]
        } else {
            // remove a cell from a subclone first and then divide
            let id2remove = if let Some(id) = self.the_only_one_subclone_present() {
                id
            } else {
                let mut variants = Variants::variant_counts(&self.subclones);
                // do not sample the cell that will proliferate
                variants[subclone_id] -= 1;
                WeightedIndex::new(variants).unwrap().sample(rng)
            };
            let _ = self.subclones[id2remove]
                .random_cell(rng)
                .with_context(|| "found empty subclone")
                .unwrap();
            if self.verbosity > 1 {
                println!("removing one cell from clone {}", id2remove);
            }

            symmetric_division(&mut self.subclones[subclone_id], rng)
                .unwrap()
                .to_vec()
        }
    }

    pub fn save_variant_fraction(&self, timepoint: usize) -> anyhow::Result<()> {
        let path2file = self
            .path2dir
            .join("variant_fraction")
            .join(timepoint.to_string());
        fs::create_dir_all(&path2file).expect("Cannot create dir");
        let path2file = path2file.join(self.id.to_string()).with_extension("csv");
        let total_variant_frac =
            Variants::variant_fractions(&self.subclones, self.compute_tot_cells());
        if self.verbosity > 1 {
            println!("total variant fraction in {:#?}", path2file)
        }

        write2file(&total_variant_frac, &path2file, None, false)?;
        Ok(())
    }

    pub fn save(&self, rng: &mut impl Rng) -> anyhow::Result<()> {
        //! Save the process.
        //! Save the SFS, SFS neutral and time.
        let path2sfs = self.path2dir.join("time");
        fs::create_dir_all(&path2sfs).expect("Cannot create dir");
        let path2file = path2sfs.join(self.id.to_string()).with_extension("csv");
        write2file(&self.time, &path2file, None, false)?;

        let path2sfs = self.path2dir.join("sfs");
        fs::create_dir_all(&path2sfs).expect("Cannot create dir");
        let path2file = path2sfs.join(self.id.to_string()).with_extension("json");
        // save subclones first
        let sfs = serde_json::to_string(&sfs(&self.subclones)).expect("cannot serialize the sfs");
        fs::write(path2file, sfs).with_context(|| "Cannot save the neutral SFS".to_string())?;

        // save SFS neutral
        let path2sfs = self.path2dir.join("sfs_neutral");
        fs::create_dir_all(&path2sfs).expect("Cannot create dir");
        let path2file = path2sfs.join(self.id.to_string()).with_extension("json");
        let cells: Vec<StemCell> = self
            .subclones
            .iter()
            .flat_map(|subclone| subclone.get_cells().to_vec())
            .collect();
        let sfs = serde_json::to_string(&Genotype::sfs_neutral(
            &cells,
            &self.distributions,
            self.verbosity,
            rng,
        ))
        .expect("Cannot serialize the neutral sfs");
        fs::write(path2file, sfs).with_context(|| "Cannot save the neutral SFS".to_string())?;

        Ok(())
    }
}

impl AdvanceStep<MAX_SUBCLONES> for HSCProcess {
    /// The id of a subclone [`CloneId`].
    type Reaction = CloneId;

    fn advance_step(&mut self, reaction: NextReaction<Self::Reaction>, rng: &mut impl Rng) {
        //! Update the process according to the next `reaction`.
        //!
        //! Select the next cell that will proliferate which belongs to the
        //! `CloneId` found by the next `reaction`.
        //! This cell will generate either a symmetric or asymmetric division.
        self.time.push(self.time.last().unwrap() + reaction.time);
        let stem_cells = self.proliferating_cells(reaction.event, rng);

        self.mutate_and_assign(reaction.event, stem_cells, rng);
    }

    fn update_state(&self, state: &mut CurrentState<MAX_SUBCLONES>) {
        for i in 0..state.population.len() {
            state.population[i] = self.subclones[i].cell_count();
        }
        // take snapshot
        let iteration = self
            .snapshot
            .iter()
            .flatten()
            .find(|&iteration| iteration == &self.counter_divisions);

        if let Some(&iteration) = iteration {
            self.save_variant_fraction(iteration)
                .expect("cannot save snapshot");
        }
    }
}

/// Distribution probabilities used to model acquisition of neutral and
/// conferring-fitness mutations upon cell division.
#[derive(Debug, Clone)]
pub struct Distributions {
    poisson: Poisson<f32>,
    bern: Bernoulli,
    bern_asymmetric: Bernoulli,
}

impl Distributions {
    pub fn new(p_asymmetric: f64, lambda_poisson: f32, p_fitness: f64, verbosity: u8) -> Self {
        if verbosity > 1 {
            println!("creating distributions with parameters asymmetric: {}, neutral_lambda: {}, p_fitness: {}", p_asymmetric, lambda_poisson, p_fitness);
        }
        Self {
            poisson: Poisson::new(lambda_poisson)
                .expect("Invalid lambda found: lambda <= 0 or nan"),
            bern: Bernoulli::new(p_fitness).expect("Invalid p: p<0 or p>1"),
            bern_asymmetric: Bernoulli::new(p_asymmetric).expect("Invalid p: p<0 or p>1"),
        }
    }

    pub fn nb_neutral_mutations(&self, rng: &mut impl Rng) -> NbPoissonNeutralMutations {
        //! The number of neutral mutations acquired upon cell division.
        let mut mutations = self.poisson.sample(rng);
        while mutations >= u8::MAX as f32 || mutations.is_sign_negative() || !mutations.is_normal()
        {
            mutations = self.poisson.sample(rng);
        }
        NbPoissonNeutralMutations::new(mutations as u8).unwrap()
    }

    pub fn acquire_p_mutation(&self, rng: &mut impl Rng) -> bool {
        self.bern.sample(rng)
    }
}

fn division(
    subclone: &mut SubClone,
    mut cell: StemCell,
    distr: &Distributions,
    proliferative_events: usize,
    rng: &mut impl Rng,
) -> Option<StemCell> {
    //! Mutate a cell: add a poisson number of neutral mutations (with rate
    //! `lambda`, see [`Distributions::new`]) and a variant conferring a
    //! fitness advantage (with probability `p`, see [`Distributions::new`]).
    //! Assign cell to `subclone` if no fit variant has been generated.
    //!
    //! ## Returns
    //! If the cell gets one fitness advantage mutation, then the function
    //! returns the cell, otherwise it returns None.
    Genotype::mutate(&mut cell, proliferative_events);
    if distr.acquire_p_mutation(rng) {
        return Some(cell);
    }
    subclone.assign_cell(cell);
    None
}

fn asymmetric_division(subclone: &mut SubClone, rng: &mut impl Rng) -> anyhow::Result<StemCell> {
    //! An asymmetric division generates a differentiated cell and a stem cell.
    //! It removes a cell from the subclone.
    subclone.random_cell(rng)
}

fn symmetric_division(
    subclone: &mut SubClone,
    rng: &mut impl Rng,
) -> anyhow::Result<[StemCell; 2]> {
    //! A symmetric division generates two stem cells belonging to `subclone`.
    //! It removes a cell from the subclone.
    //!
    //! Since we are simulating a Moran process, another stem cell must be
    //! removed from the system each time a symmetric division occurs.
    ensure!(!subclone.is_empty());
    let cell1 = subclone
        .random_cell(rng)
        .with_context(|| "found empty subclone")?;
    let cell2 = cell1.clone();

    Ok([cell1, cell2])
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
        Distributions::new(1., lambda_poisson.0, 0., 0).nb_neutral_mutations(&mut rng);
    }

    #[quickcheck]
    fn division_no_new_clone(
        lambda: LambdaFromNonZeroU8,
        proliferative_events: usize,
        seed: u64,
    ) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let distr = Distributions::new(1., lambda.0, 0f64, 0);

        let mut neutral_clone = SubClone::new(1, 2);
        let cell = StemCell::new();

        division(
            &mut neutral_clone,
            cell,
            &distr,
            proliferative_events,
            &mut rng,
        )
        .is_none()
    }

    #[quickcheck]
    fn division_new_clone(
        lambda: LambdaFromNonZeroU8,
        proliferative_events: usize,
        seed: u64,
    ) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let distr = Distributions::new(1., lambda.0, 1f64, 0);

        let mut neutral_clone = SubClone::new(1, 2);
        let cell = StemCell::new();

        division(
            &mut neutral_clone,
            cell,
            &distr,
            proliferative_events,
            &mut rng,
        )
        .is_some()
    }

    #[test]
    #[should_panic]
    fn asymmetric_should_panic_empty_subclone_test() {
        let mut rng = ChaCha8Rng::seed_from_u64(26);
        let mut neutral_clone = SubClone::new(1, 2);
        asymmetric_division(&mut neutral_clone, &mut rng).unwrap();
    }

    #[quickcheck]
    fn asymmetric_division_test(seed: u64, id: usize) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let mut neutral_clone = SubClone::new(id, 2);
        let cell = StemCell::new();
        neutral_clone.assign_cell(cell);

        let stem_cell = asymmetric_division(&mut neutral_clone, &mut rng).unwrap();
        neutral_clone.cell_count() == 0 && !stem_cell.has_mutations()
    }

    #[test]
    #[should_panic]
    fn symmetric_should_panic_empty_subclone_test() {
        let mut rng = ChaCha8Rng::seed_from_u64(26);
        let mut neutral_clone = SubClone::new(1, 2);
        symmetric_division(&mut neutral_clone, &mut rng).unwrap();
    }

    #[quickcheck]
    fn symmetric_division_test(seed: u64, id: usize) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let mut neutral_clone = SubClone::new(id, 2);
        let cell = StemCell::new();
        neutral_clone.assign_cell(cell);
        assert_eq!(neutral_clone.cell_count(), 1);

        let stem_cells = symmetric_division(&mut neutral_clone, &mut rng).unwrap();
        neutral_clone.cell_count() == 0
            && stem_cells.iter().all(|cell| !cell.has_mutations())
            && stem_cells.len() == 2
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

        dbg!(&hsc.subclones[subclone_id as usize]);
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
}
