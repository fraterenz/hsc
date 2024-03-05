use rand::Rng;
use rand_distr::{Bernoulli, Distribution};

use crate::{
    stemcell::{assign_background_mutations, assign_divisional_mutations},
    subclone::{assign_fit_mutations, proliferating_cell, CloneId, Distributions, SubClones},
};

#[derive(Debug, Clone, Default)]
/// Stem cells can either self-renew by symmetric division or differentiate by
/// asymmetric division.
pub enum Division {
    /// A stem cell proliferates and gives rise to a differentiate cell (which
    /// we discard in our simulations)
    Asymmetric(Bernoulli),
    /// A stem cell proliferates and gives to a new stem cell with the same
    /// genetic material as the mother cell plus some additional mutations
    /// acquired upon division (fit divisions and neutral ones)
    #[default]
    Symmetric,
}

#[derive(Clone, Debug, Default)]
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
/// with `c1` only if we simulated a symmetric division, see [`Division`].
pub struct Proliferation {
    pub neutral_mutation: NeutralMutations,
    pub division: Division,
}

impl Proliferation {
    pub fn new(neutral_mutation: NeutralMutations, division: Division) -> Self {
        Proliferation {
            neutral_mutation,
            division,
        }
    }

    pub fn proliferate(
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
            println!(
                "cell with last division time {} is dividing",
                stem_cell.last_division_t
            );
            if verbosity > 2 {
                println!("cell {:#?} is dividing", stem_cell);
            }
        }
        if let NeutralMutations::UponDivisionAndBackground = self.neutral_mutation {
            assign_background_mutations(
                &mut stem_cell,
                time,
                &distributions.neutral_poisson,
                rng,
                verbosity,
            );
        }
        let cells = match self.division {
            Division::Asymmetric(p) => {
                // true means asymmetric
                if p.sample(rng) {
                    if verbosity > 1 {
                        println!("asymmetric division");
                    }
                    vec![stem_cell]
                } else {
                    vec![stem_cell.clone(), stem_cell]
                }
            }
            Division::Symmetric => vec![stem_cell.clone(), stem_cell],
        };
        for mut stem_cell in cells {
            assign_divisional_mutations(
                &mut stem_cell,
                &distributions.neutral_poisson,
                rng,
                verbosity,
            );
            assign_fit_mutations(
                subclones,
                proliferating_subclone,
                stem_cell,
                distributions,
                rng,
                verbosity,
            );
        }
    }
}

#[derive(Clone, Debug, Default)]
/// Specifies the kind of neutral mutations to simulate.
pub enum NeutralMutations {
    /// Neutral mutations upon division due to errors in the DNA duplication
    UponDivision,
    /// Neutral mutations upon division and background mutations (mutations
    /// appearing during the lifetime of the cell)
    #[default]
    UponDivisionAndBackground,
}
