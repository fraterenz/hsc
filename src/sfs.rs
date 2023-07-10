use crate::{neutral::StemCell, MAX_SUBCLONES};
use anyhow::ensure;
use rand::Rng;
use sosa::NbIndividuals;
use std::collections::HashMap;

/// Id of the [`SubClone`]s used for the computation of the [`sfs`].
pub type CloneId = usize;

#[derive(Debug, Clone)]
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
}

impl SubClone {
    pub fn new(id: CloneId, cell_capacity: usize) -> SubClone {
        SubClone {
            cells: Vec::with_capacity(cell_capacity),
            id,
        }
    }

    pub fn get_cells(&self) -> &[StemCell] {
        &self.cells
    }

    pub fn is_empty(&self) -> bool {
        self.get_cells().is_empty()
    }

    pub fn assign_cell(&mut self, cell: StemCell) {
        self.cells.push(cell);
    }

    pub fn random_cell(&mut self, rng: &mut impl Rng) -> anyhow::Result<StemCell> {
        ensure!(!self.cells.is_empty());
        Ok(self.cells.swap_remove(rng.gen_range(0..self.cells.len())))
    }

    pub fn cell_count(&self) -> u64 {
        self.cells.len() as u64
    }
}

pub fn sfs(subclones: &[SubClone]) -> HashMap<CloneId, NbIndividuals> {
    //! Compute the site frequency spectrum for clones that have a
    //! proliferative advantage.
    //! Each key of the sfs is the `CloneId` while the value is the number of
    //! cells corresponding the to clone identified by the key.
    let mut clones = HashMap::<CloneId, NbIndividuals>::with_capacity(MAX_SUBCLONES);
    for clone in subclones {
        clones.insert(clone.id, clone.cell_count());
    }
    clones
}

//
// keep only (filter) cells that have the good clone
// random select one cell from those cels
// add mutations
// check if assign clone

// fn main() {
//
//
//     asymmetric_division(&mut cell_1);
//     let population = vec![cell_0, cell_1];
//     println!("{:#?}", compute_reactions(&population));
//     println!("{:#?}", &population);
// }
//
#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck_macros::quickcheck;

    #[quickcheck]
    fn assign_cell_test(id: usize) -> bool {
        let mut neutral_clone = SubClone { cells: vec![], id };
        let cell = StemCell::new();
        assert!(neutral_clone.cells.is_empty());

        neutral_clone.assign_cell(cell);
        !neutral_clone.cells.is_empty()
    }

    #[quickcheck]
    fn sfs_test(id0: usize) -> bool {
        let mut clone0 = SubClone {
            cells: Vec::with_capacity(2),
            id: id0,
        };
        let id1 = id0 / 2 + 2;
        let clone1 = SubClone {
            cells: Vec::with_capacity(2),
            id: id1,
        };

        let cell = StemCell::new();

        clone0.assign_cell(cell);

        let sfs = sfs(&[clone0, clone1]);
        sfs[&id0] == 1 && sfs[&id1] == 0
    }
}
