use anyhow::{ensure, Context};
use rand::Rng;

use crate::genotype::{NeutralMutationPoisson, Variant};

/// Hematopoietic stem and progenitor cells (HSPCs) are a rare population of
/// precursor cells that possess the capacity for self-renewal and multilineage
/// differentiation.
///
/// They carry a set of neutral mutations and are assigned to [`crate::subclone::SubClone`].
#[derive(Debug, Clone)]
pub struct StemCell {
    pub variants: Vec<Variant>,
    /// the last time at which the cell has divided
    last_division_t: f32,
}

impl Default for StemCell {
    fn default() -> Self {
        //! Creates a stem cell without any neutral mutations.
        Self::new()
    }
}

impl StemCell {
    pub fn new() -> StemCell {
        //! Construct a new cell without any neutral mutations.
        StemCell {
            variants: Vec::new(),
            last_division_t: 0.,
        }
    }

    pub fn with_mutations(mutations: Vec<Variant>) -> StemCell {
        assert!(!mutations.is_empty());
        let mut cell = StemCell::new();
        cell.variants = mutations;
        cell
    }

    pub fn has_mutations(&self) -> bool {
        !self.variants.is_empty()
    }

    pub fn burden(&self) -> usize {
        self.variants.len()
    }

    pub fn get_last_division_time(&mut self) -> &f32 {
        &self.last_division_t
    }

    pub fn set_last_division_time(&mut self, last_division_t: f32) -> anyhow::Result<()> {
        ensure!(last_division_t.is_sign_positive());
        self.last_division_t = last_division_t;
        Ok(())
    }

    pub fn interdivision_time(&self, time: f32) -> anyhow::Result<f32> {
        ensure!(
            time >= self.last_division_t,
            "found a cell that has divided in the future! {} last_division vs {} time",
            self.last_division_t,
            time,
        );
        Ok(time - self.last_division_t)
    }
}

fn mutate(cell: &mut StemCell, mut mutations: Vec<Variant>) {
    cell.variants.append(&mut mutations);
}

pub fn assign_divisional_mutations(
    stem_cell: &mut StemCell,
    neutral_poisson: &NeutralMutationPoisson,
    rng: &mut impl Rng,
    verbosity: u8,
) {
    let mutations = neutral_poisson.new_muts_upon_division(rng);
    if let Some(mutations) = mutations {
        if verbosity > 2 {
            println!("assigning {mutations:#?} to cell {stem_cell:#?}");
        }
        mutate(stem_cell, mutations);
    } else if verbosity > 2 {
        println!("no mutations to assign to cell {stem_cell:#?}");
    }
}

pub fn assign_background_mutations(
    stem_cell: &mut StemCell,
    time: f32,
    neutral_poisson: &NeutralMutationPoisson,
    rng: &mut impl Rng,
    verbosity: u8,
) {
    //! Assign background mutations to all cells in the system based on the
    //! current simulation time.
    //!
    //! This updates also the time of the last division for the cell.
    let interdivison_time = stem_cell
        .interdivision_time(time)
        .with_context(|| "wrong interdivision time")
        .unwrap();
    if verbosity > 1 {
        println!("assigning background mutations with interdivision time {interdivison_time}");
    }
    // 2. draw background mutations and assign them to `c`
    if let Some(background) = neutral_poisson.new_muts_background(interdivison_time, rng, verbosity)
    {
        if verbosity > 2 {
            println!(
                "assigning {} background mutations to cell {:#?}",
                background.len(),
                stem_cell
            )
        }
        mutate(stem_cell, background);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck_macros::quickcheck;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;
    use std::num::NonZeroU8;
    use uuid::Uuid;

    #[should_panic]
    #[test]
    fn new_cell_with_empty_mutations_test() {
        StemCell::with_mutations(vec![]);
    }

    #[quickcheck]
    fn new_cell_with_mutations_test(nb_mutations: NonZeroU8) -> bool {
        let mutations = (0..nb_mutations.get()).map(|_| Uuid::new_v4()).collect();
        let cell = StemCell::with_mutations(mutations);
        cell.has_mutations() && nb_mutations.get() as usize == cell.burden()
    }

    #[quickcheck]
    fn mutate_test(nb_mutations: NonZeroU8) -> bool {
        let mut cell = StemCell::new();
        let mutations = (0..nb_mutations.get()).map(|_| Uuid::new_v4()).collect();
        mutate(&mut cell, mutations);
        nb_mutations.get() as usize == cell.burden()
    }

    #[quickcheck]
    fn assign_background_mutations_test(seed: u64) -> bool {
        let rng = &mut ChaCha8Rng::seed_from_u64(seed);
        let mutations = vec![Variant::new_v4()];
        let time = 9.1;
        let mut stem_cell = StemCell::with_mutations(mutations.clone());
        let poissons = NeutralMutationPoisson::new(1.1, 12f32).unwrap();

        assign_background_mutations(&mut stem_cell, time, &poissons, rng, 0);
        mutations != stem_cell.variants && mutations.len() < stem_cell.variants.len()
    }

    #[quickcheck]
    fn interdivision_time_test(time: NonZeroU8) -> bool {
        let mut stem_cell = StemCell::new();
        stem_cell.last_division_t = time.get() as f32;

        (stem_cell.interdivision_time(time.get() as f32).unwrap()).abs() < f32::EPSILON
    }

    #[test]
    #[should_panic]
    fn interdivision_time_panic_test() {
        let mut stem_cell = StemCell::new();
        stem_cell.last_division_t = 1.;

        stem_cell.interdivision_time(0.1).unwrap();
    }
}
