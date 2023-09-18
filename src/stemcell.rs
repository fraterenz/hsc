use crate::genotype::Variant;

/// Hematopoietic stem and progenitor cells (HSPCs) are a rare population of
/// precursor cells that possess the capacity for self-renewal and multilineage
/// differentiation.
///
/// They carry a set of neutral mutations and are assigned to [`crate::subclone::SubClone`].
#[derive(Debug, Clone)]
pub struct StemCell {
    pub variants: Vec<Variant>,
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
}

pub fn mutate(cell: &mut StemCell, mut mutations: Vec<Variant>) {
    cell.variants.append(&mut mutations);
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck_macros::quickcheck;
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
}
