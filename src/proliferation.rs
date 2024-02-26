use enum_dispatch::enum_dispatch;
use rand::Rng;

use crate::{
    genotype::Variant,
    stemcell::{assign_background_mutations, mutate, StemCell},
    subclone::{assign, proliferating_cell, CloneId, Distributions, SubClones},
};

fn assign_mutations(stem_cell: &mut StemCell, mutations: Option<Vec<Variant>>, verbosity: u8) {
    if let Some(mutations) = mutations {
        if verbosity > 2 {
            println!("assigning {:#?} to cell {:#?}", mutations, stem_cell);
        }
        mutate(stem_cell, mutations);
    } else if verbosity > 2 {
        println!("no mutations to assign to cell {:#?}", stem_cell);
    }
}

fn proliferation(
    stem_cell: StemCell,
    subclones: &mut SubClones,
    proliferating_subclone: CloneId,
    distributions: &Distributions,
    rng: &mut impl Rng,
    verbosity: u8,
) {
    if verbosity > 1 {
        println!("cell {:#?} is dividing", stem_cell);
    }

    let mut new_cell = stem_cell.clone();
    // the new cell hasn't divided yet
    new_cell.last_division_t = 0f32;
    for mut cell in [new_cell, stem_cell] {
        let division = distributions.neutral_poisson.new_muts_upon_division(rng);
        if verbosity > 2 {
            println!(
                "assigning {} mutations upon cell division to cell {:#?}",
                division.as_ref().unwrap_or(&vec![]).len(),
                cell
            );
        }
        assign_mutations(&mut cell, division, verbosity);
        assign(
            subclones,
            proliferating_subclone,
            cell,
            distributions,
            rng,
            verbosity,
        );
    }
}

#[derive(Clone, Debug)]
/// Simulate fit mutations, neutral background mutations and neutral mutations
/// upon division.
pub struct DivisionAndBackgroundMutationsProliferation;

#[derive(Clone, Debug)]
/// Do not simulate neutral background mutations but only neutral mutations
/// and fit mutations upon division
pub struct DivisionMutationsProliferation;

#[enum_dispatch]
/// All the updates and changes in the system upon proliferation.
///
/// The proliferation step is implemented as following:
///
/// 1. select the cell `c` that will proliferate next from the clone
/// with id `reaction` determined by the Gillespie algorithm
///
/// 2. draw and assign mb neutral background mutations to `c` from
/// Poisson(DeltaT * mub), where mub is the backgound mutation rate
///
/// 3. draw and assign md neutral division mutations to `c` from
/// Poisson(mud), where mud is the division mutation rate
///
/// 4. draw from a Bernoulli trial a fit mutation and in case assign
/// `c` to a new clone. This clone must be empty and different from the
/// old clone which `c` belonged to.
///
/// 5. clone the proliferating cell `c` into `c1` and repeat step 3, 4
/// with `c1`.
pub trait Proliferation {
    fn duplicate_cell_assign_mutations(
        &self,
        subclones: &mut SubClones,
        time: f32,
        proliferating_subclone: CloneId,
        distributions: &Distributions,
        rng: &mut impl Rng,
        verbosity: u8,
    );
}

#[derive(Clone, Debug)]
#[enum_dispatch(Proliferation)]
/// Duplicate cells and deal with subclones.
/// Also specifies whether to simulate background neutral mutations or only
/// mutations upon divisions, see [`Proliferation`].
pub enum MutateUponDivision {
    DivisionAndBackgroundMutations(DivisionAndBackgroundMutationsProliferation),
    DivisionMutations(DivisionMutationsProliferation),
}

impl Default for MutateUponDivision {
    fn default() -> Self {
        MutateUponDivision::DivisionAndBackgroundMutations(
            DivisionAndBackgroundMutationsProliferation,
        )
    }
}

impl Proliferation for DivisionMutationsProliferation {
    fn duplicate_cell_assign_mutations(
        &self,
        subclones: &mut SubClones,
        time: f32,
        proliferating_subclone: CloneId,
        distributions: &Distributions,
        rng: &mut impl Rng,
        verbosity: u8,
    ) {
        let stem_cell = proliferating_cell(subclones, proliferating_subclone, verbosity, rng);
        if verbosity > 1 {
            println!("proliferation at time {}", time);
        }
        proliferation(
            stem_cell,
            subclones,
            proliferating_subclone,
            distributions,
            rng,
            verbosity,
        )
    }
}

impl Proliferation for DivisionAndBackgroundMutationsProliferation {
    fn duplicate_cell_assign_mutations(
        &self,
        subclones: &mut SubClones,
        time: f32,
        proliferating_subclone: CloneId,
        distributions: &Distributions,
        rng: &mut impl Rng,
        verbosity: u8,
    ) {
        let mut stem_cell = proliferating_cell(subclones, proliferating_subclone, verbosity, rng);
        if verbosity > 1 {
            println!("proliferation at time {}", time);
        }
        assign_background_mutations(
            &mut stem_cell,
            time,
            &distributions.neutral_poisson,
            rng,
            verbosity,
        );
        proliferation(
            stem_cell,
            subclones,
            proliferating_subclone,
            distributions,
            rng,
            verbosity,
        )
    }
}
