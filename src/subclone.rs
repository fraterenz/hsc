use crate::stemcell::StemCellId;
use anyhow::{bail, ensure};
use rand::{seq::SliceRandom, Rng};

/// Id of the [`SubClone`]s.
pub type CloneId = usize;

/// A group of cells sharing the same genetic background with a specific
/// proliferation rate.
///
/// The main loop of the simulation delegates the proliferation of cells to
/// this structure, meaning that the `SubClone` will randomly pick one of its
/// cells and make it proliferate.
/// Upon proliferation, the cell can be assigned to a new clone with
/// probability `p` (see [`crate::process::CellDivisionProbabilities`]).
#[derive(Debug, Clone)]
pub struct SubClone {
    // TODO: transform this into reference
    cells: Vec<StemCellId>,
    pub id: CloneId,
}

impl SubClone {
    pub fn new(id: CloneId, cell_capacity: usize) -> SubClone {
        SubClone {
            cells: Vec::with_capacity(cell_capacity),
            id,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }

    pub fn assign_cell(&mut self, cell: StemCellId) {
        self.cells.push(cell);
    }

    fn gen_rand_id(&self, rng: &mut impl Rng) -> StemCellId {
        self.cells
            .choose(rng)
            .expect("found empty subclone")
            .to_owned()
    }

    fn remove_cell_with_id(&mut self, id: &StemCellId) -> anyhow::Result<StemCellId> {
        if let Some(index) = self.cells.iter().position(|cellid| cellid == id) {
            return Ok(self.cells.swap_remove(index));
        }
        bail!(format!("cannot find cell with id {} in subclone", id))
    }

    pub fn random_cell(&mut self, rng: &mut impl Rng) -> anyhow::Result<StemCellId> {
        //! Pick a random cell from the clone
        ensure!(!self.cells.is_empty());
        let id = self.gen_rand_id(rng);
        self.remove_cell_with_id(&id)
    }

    pub fn cell_count(&self) -> u64 {
        self.cells.len() as u64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck_macros::quickcheck;

    #[quickcheck]
    fn assign_cell_test(id: usize) -> bool {
        let mut neutral_clone = SubClone { cells: vec![], id };
        let cell = StemCellId::new_v4();
        assert!(neutral_clone.cells.is_empty());

        neutral_clone.assign_cell(cell);
        !neutral_clone.cells.is_empty()
    }
}
