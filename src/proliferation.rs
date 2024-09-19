use anyhow::{ensure, Context};
use rand::Rng;
use rand_distr::{Bernoulli, Distribution};

use crate::{
    stemcell::{assign_background_mutations, assign_divisional_mutations, StemCell, StemCellId},
    subclone::{next_clone, proliferating_cell, CloneId, Distributions, SubClones},
    tree::PhyloTree,
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
    // counter of cells used to generate new cells
    last_id: StemCellId,
}

impl Proliferation {
    pub fn new(
        neutral_mutation: NeutralMutations,
        division: Division,
        cells: &[&StemCell],
    ) -> Self {
        //! Create a new Proliferation scheme.
        //!
        //! ## Panics
        //! When the slice `cells` is empty.
        assert!(!cells.is_empty());
        let last_id = cells.iter().map(|&cell| cell.id).max().unwrap();
        Proliferation {
            neutral_mutation,
            division,
            last_id,
        }
    }

    pub fn proliferate(
        &mut self,
        subclones: &mut SubClones,
        tree: &mut PhyloTree,
        time: f32,
        proliferating_subclone: CloneId,
        distributions: &Distributions,
        rng: &mut impl Rng,
        verbosity: u8,
    ) -> anyhow::Result<()> {
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
        // the fit variant is sampled assuming a background mutation
        let clone_id = next_clone(subclones, proliferating_subclone, p, rng, verbosity);
        if verbosity > 1 {
            println!("proliferation at time {}", time);
            println!(
                "cell with id {} and last division time {} is dividing",
                stem_cell.id,
                stem_cell.get_last_division_time()
            );
            if verbosity > 2 {
                println!("cell {:#?} is dividing", stem_cell);
            }
        }
        // the mutations of the ancestral cell TODO THINK
        let nb_muts_old = stem_cell.variants.len();
        if let NeutralMutations::UponDivisionAndBackground = self.neutral_mutation {
            assign_background_mutations(
                &mut stem_cell,
                time,
                &distributions.neutral_poisson,
                rng,
                verbosity,
            );
            // update node edge with the background mutations
            tree.assign_background_muts_to_tree(
                stem_cell.id,
                stem_cell.variants.len() - nb_muts_old,
                verbosity,
            )?;
        }

        // perform the division
        let cells = match self.division {
            Division::Asymmetric(prob) => {
                if prob.sample(rng) {
                    if verbosity > 1 {
                        println!("asymmetric division");
                    }
                    assign_divisional_mutations(
                        &mut stem_cell,
                        &distributions.neutral_poisson,
                        rng,
                        verbosity,
                    );
                    let mutational_distance = (stem_cell.variants.len() - nb_muts_old) as f64;
                    ensure!(
                        tree.put_cell_on_new_node(stem_cell.id, mutational_distance, verbosity)?,
                        "the proliferating cell is not registered as a leaf"
                    );
                    vec![stem_cell]
                } else {
                    let new_cell = self.symmtetric_division_helper(
                        &mut stem_cell,
                        tree,
                        distributions,
                        rng,
                        verbosity,
                    )?;
                    vec![stem_cell, new_cell]
                }
            }
            Division::Symmetric => {
                let new_cell = self.symmtetric_division_helper(
                    &mut stem_cell,
                    tree,
                    distributions,
                    rng,
                    verbosity,
                )?;
                vec![stem_cell, new_cell]
            }
        };
        // assign new cells to subclones: if new fit variant present, move cell
        // to a new subclone, otherwise reassign both cells to the subclone of
        // the ancestral cell (i.e. assume fit variants are only background
        // mutations)
        for mut stem_cell in cells {
            stem_cell
                .set_last_division_time(time)
                .with_context(|| "wrong time")
                .unwrap();
            subclones
                .get_mut_clone_unchecked(clone_id)
                .assign_cell(stem_cell);
        }
        //dbg!(tree.tree.to_newick().unwrap());
        Ok(())
    }

    fn symmtetric_division_helper(
        &mut self,
        stem_cell: &mut StemCell,
        tree: &mut PhyloTree,
        distributions: &Distributions,
        rng: &mut impl Rng,
        verbosity: u8,
    ) -> anyhow::Result<StemCell> {
        //! Creates a new leaf node to which `stem_cell` will be assigned with
        //! a distance equal to the new mutations acquired.
        //! The background mutations are not simulated but assumed to be
        //! acquired somewhere else, whereas the divisional ones are simulated
        //! here.
        //!
        //! Creates also a new cell by cloning the ancestral cell (i.e.
        //! `stem_cell` before assigning divisional mutations) and put this
        //! cell on a leaf node whose parent is the node of the ancestral cell.
        //! As a results, the new cell will be a sibling of the new `stem_cell`,
        //! i.e. the cell with divisional mutations assigned to a new node.
        //!
        //! ## Returns
        //! A new cell which is a sibling of `stem_cell` and daughter of the
        //! ancestral cell shared with `stem_cell`.
        // store mutations before division which are the same for both cells
        let nb_muts_old = stem_cell.variants.len();
        // clone before assigning divisional mutations
        let mut new_cell = stem_cell.clone();
        self.last_id += 1;
        new_cell.id = self.last_id;

        assign_divisional_mutations(stem_cell, &distributions.neutral_poisson, rng, verbosity);
        if verbosity > 1 {
            println!(
                "putting dividing cell with id {} into new node",
                stem_cell.id
            );
        }
        ensure!(
            tree.put_cell_on_new_node(
                stem_cell.id,
                (stem_cell.variants.len() - nb_muts_old) as f64,
                verbosity
            )?,
            "the proliferating cell is not registered as a leaf"
        );

        assign_divisional_mutations(
            &mut new_cell,
            &distributions.neutral_poisson,
            rng,
            verbosity,
        );
        if verbosity > 1 {
            println!(
                "adding the sibling cell with id {} on a new node",
                self.last_id
            );
        }
        // a new cell is not registered as a leaf node yet
        ensure!(
            !tree.assign_sibling(
                new_cell.id,
                stem_cell.id,
                (new_cell.variants.len() - nb_muts_old) as f64,
                verbosity
            )?,
            "the new cell is already registered as a leaf"
        );
        Ok(new_cell)
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
