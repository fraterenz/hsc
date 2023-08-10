use anyhow::Context;
use rustc_hash::FxHashSet;
use serde::{Deserialize, Serialize};
use std::{collections::HashSet, fs, path::Path};

use crate::genotype::GenotypeId;

/// Hematopoietic stem and progenitor cells (HSPCs) are a rare population of
/// precursor cells that possess the capacity for self-renewal and multilineage
/// differentiation.
///
/// They carry a set of neutral mutations and are assigned to [`crate::subclone::SubClone`].
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(transparent)]
pub struct StemCell {
    /// A collection of ids which identify the iterations upon which the cell
    /// prolfierates during the simulation.
    pub proliferation_events_id: FxHashSet<GenotypeId>,
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
            proliferation_events_id: FxHashSet::default(),
        }
    }

    pub fn with_set_of_mutations(mutation_set: Vec<usize>) -> StemCell {
        //! Create a stem cell with a set of neutral mutations.
        //!
        //! The mutation set is not a collection of mutations but a collection
        //! of genotypes that will be converted later on into mutations, see
        //! [`crate::genotype::StatisticsMutations`].
        let mut cell = StemCell::new();
        let mutation_set = HashSet::from_iter(mutation_set.into_iter());
        cell.proliferation_events_id = mutation_set;
        cell
    }

    pub fn has_mutations(&self) -> bool {
        !self.proliferation_events_id.is_empty()
    }

    pub fn record_division(&mut self, event_id: usize) {
        //! Record the iteration when the cell has undergone cell-division
        //! (proliferation).
        //!
        //! Since we store the divisions performed by each cell, instead of the
        //! mutations (to avoid generating a random number at every division?),
        //! we mutate a cell by storing the iterations at which the cell
        //! proliferates.
        //! The mapping between genotypes (i.e. proliferation id) and mutations
        //! is performed by [`crate::genotype::StatisticsMutations`].
        self.proliferation_events_id.insert(event_id);
    }
}

pub fn load_cells(path2file: &Path) -> anyhow::Result<Vec<StemCell>> {
    Ok(serde_json::from_slice(&fs::read(path2file)?)?)
}

pub fn save_cells(cells: &[StemCell], path2file: &Path) -> anyhow::Result<()> {
    let cells_sr = serde_json::to_string(cells).with_context(|| "cannot serialize cells")?;
    fs::write(path2file, cells_sr)
        .with_context(|| format!("Cannot save cells to {:#?}", path2file))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck_macros::quickcheck;

    #[quickcheck]
    fn genotype_mutate_test(divisions: usize) -> bool {
        let mut cell = StemCell::new();
        cell.record_division(divisions);
        cell.has_mutations() && cell.proliferation_events_id.len() == 1
    }
}
