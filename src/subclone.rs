use crate::stemcell::StemCell;
use anyhow::ensure;
use rand::Rng;

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
}
