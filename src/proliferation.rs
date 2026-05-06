use anyhow::Context;
use log::{debug, trace};
use rand::Rng;
use rand_distr::{Bernoulli, Distribution};

use crate::{
    stemcell::{assign_background_mutations, assign_divisional_mutations, StemCell},
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
        let mut stem_cell = proliferating_cell(subclones, proliferating_subclone, rng);
        let p = (distributions.u
            * stem_cell
                .interdivision_time(time)
                .with_context(|| "wrong interdivision time")
                .unwrap()) as f64;
        // the fit variant is sampled here
        let clone_id = next_clone(subclones, proliferating_subclone, p, rng);
        debug!("proliferation at time {time}");
        debug!(
            "cell with last division time {} is dividing",
            stem_cell.get_last_division_time()
        );
        trace!("cell {stem_cell:#?} is dividing");
        self.realise_background_for_cell(&mut stem_cell, time, distributions, rng);
        let cells = match self.division {
            Division::Asymmetric(prob) => {
                if prob.sample(rng) {
                    debug!("asymmetric division");
                    vec![stem_cell]
                } else {
                    vec![stem_cell.clone(), stem_cell]
                }
            }
            Division::Symmetric => vec![stem_cell.clone(), stem_cell],
        };
        for mut stem_cell in cells {
            assign_divisional_mutations(&mut stem_cell, &distributions.neutral_poisson, rng);
            stem_cell
                .set_last_division_time(time)
                .with_context(|| "wrong time")
                .unwrap();
            subclones
                .get_mut_clone_unchecked(clone_id)
                .assign_cell(stem_cell);
        }
    }

    pub fn realise_background_mutations(
        &self,
        subclones: &mut SubClones,
        time: f32,
        distributions: &Distributions,
        rng: &mut impl Rng,
    ) {
        //! Force the lazily-deferred background-mutation draw for every cell
        //! whose `last_division_t` lies strictly before `time`.
        //!
        //! Background mutations are normally drawn at the next division (see
        //! [`Proliferation::proliferate`]); snapshots and phase transitions
        //! need to read mutation state without waiting for that division, so
        //! they call into this method to bring every cell up to the current
        //! `time`. Cells whose last division is at or after `time` are left
        //! untouched.
        //!
        //! No-op when `self.neutral_mutation` is [`NeutralMutations::UponDivision`].
        debug!("realise and assign the neutral background mutations to all cells");
        for stem_cell in subclones.get_mut_cells() {
            self.realise_background_for_cell(stem_cell, time, distributions, rng);
        }
    }

    fn realise_background_for_cell(
        &self,
        cell: &mut StemCell,
        time: f32,
        distributions: &Distributions,
        rng: &mut impl Rng,
    ) {
        //! Per-cell variant of [`Self::realise_background_mutations`], shared
        //! by the bulk sweep and by [`Self::proliferate`] (which operates on a
        //! single cell detached from [`SubClones`]).
        if let NeutralMutations::UponDivisionAndBackground = self.neutral_mutation {
            if cell.get_last_division_time() < &time {
                assign_background_mutations(cell, time, &distributions.neutral_poisson, rng);
            }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stemcell::StemCell;
    use crate::Probs;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    fn make_distributions() -> Distributions {
        Distributions::new(Probs::new(50., 50., 0.01, 0., 100))
    }

    fn make_subclones(n: usize, last_division_t: f32) -> SubClones {
        let mut cells = vec![StemCell::new(); n];
        for cell in &mut cells {
            cell.set_last_division_time(last_division_t).unwrap();
        }
        SubClones::new(cells, n)
    }

    #[test]
    fn realise_background_mutations_is_noop_for_upon_division() {
        let rng = &mut ChaCha8Rng::seed_from_u64(42);
        let mut subclones = make_subclones(5, 0.0);
        let proliferation = Proliferation::new(NeutralMutations::UponDivision, Division::Symmetric);

        proliferation.realise_background_mutations(
            &mut subclones,
            10.0,
            &make_distributions(),
            rng,
        );

        for cell in subclones.get_neutral_clone().get_cells() {
            assert!(cell.mutations.is_empty());
        }
    }

    #[test]
    fn realise_background_mutations_assigns_with_background_enabled() {
        let rng = &mut ChaCha8Rng::seed_from_u64(42);
        let mut subclones = make_subclones(5, 0.0);
        let proliferation = Proliferation::new(
            NeutralMutations::UponDivisionAndBackground,
            Division::Symmetric,
        );

        proliferation.realise_background_mutations(
            &mut subclones,
            10.0,
            &make_distributions(),
            rng,
        );

        let total: usize = subclones
            .get_neutral_clone()
            .get_cells()
            .iter()
            .map(|c| c.burden())
            .sum();
        assert!(total > 0);
    }

    #[test]
    fn realise_background_mutations_skips_cells_past_time() {
        // last_division_t is in the future relative to the realise time:
        // the method must skip those cells without panicking.
        let rng = &mut ChaCha8Rng::seed_from_u64(42);
        let mut subclones = make_subclones(3, 5.0);
        let proliferation = Proliferation::new(
            NeutralMutations::UponDivisionAndBackground,
            Division::Symmetric,
        );

        proliferation.realise_background_mutations(&mut subclones, 1.0, &make_distributions(), rng);

        for cell in subclones.get_neutral_clone().get_cells() {
            assert!(cell.mutations.is_empty());
        }
    }
}
