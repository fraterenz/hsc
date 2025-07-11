use anyhow::Context;
use rand::Rng;
use rand_distr::{Bernoulli, Distribution};

use crate::{
    stemcell::{assign_background_mutations, assign_divisional_mutations},
    subclone::{next_clone, proliferating_cell, CloneId, Distributions, SubClones},
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
///    with id `reaction` determined by the Gillespie algorithm
///
/// 2. draw and assign mb neutral background mutations to `c` from
///    Poisson(DeltaT * mub), where mub is the backgound mutation rate
///
/// 3. draw and assign md neutral division mutations to `c` from
///    Poisson(mud), where mud is the division mutation rate
///
/// 4. draw from a Bernoulli trial a fit mutation and in case assign
///    `c` to a new clone. This clone must be empty and different from the
///    old clone which `c` belonged to.
///
/// 5. clone the proliferating cell `c` into `c1` and repeat step 3, 4
///    with `c1` only if we simulated a symmetric division, see [`Division`].
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
        //! Sample a random cell from the subclones with id
        //! `proliferating_subclone`, assign neutral mutations (background and
        //! those upon proliferation) and finally assign `cell` to a subclone.
        //!
        //! For the last step, there are two possible scenarios:
        //!
        //! 1. if there is a new fit variant, then `cell` will be assigned to a
        //!    new **empty** random clone with an id different from `old_subclone_id`
        //!    and panics if there aren't any empty subclones left
        //! 2. else, reassign `cell` to the old subclone with id `old_subclone_id`
        let mut stem_cell = proliferating_cell(subclones, proliferating_subclone, verbosity, rng);
        let p = (distributions.u
            * stem_cell
                .interdivision_time(time)
                .with_context(|| "wrong interdivision time")
                .unwrap()) as f64;
        // the fit variant is sampled here
        let clone_id = next_clone(subclones, proliferating_subclone, p, rng, verbosity);
        if verbosity > 1 {
            println!("proliferation at time {time}");
            println!(
                "cell with last division time {} is dividing",
                stem_cell.get_last_division_time()
            );
            if verbosity > 2 {
                println!("cell {stem_cell:#?} is dividing");
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
            Division::Asymmetric(prob) => {
                if prob.sample(rng) {
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
            stem_cell
                .set_last_division_time(time)
                .with_context(|| "wrong time")
                .unwrap();
            subclones
                .get_mut_clone_unchecked(clone_id)
                .assign_cell(stem_cell);
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
