use anyhow::Context;
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use uuid::Uuid;

pub type Variant = Uuid;
pub type StemCellId = Uuid;

#[derive(Clone, Debug)]
pub struct StemCellPopulation(FxHashMap<StemCellId, StemCell>);

impl StemCellPopulation {
    pub fn new(cells: Vec<StemCell>) -> Self {
        //! Create a population of cells, with all cells being unrelated, and
        //! assign a random idx.
        let mut population = FxHashMap::default();
        population.shrink_to(cells.len());
        for cell in cells.into_iter() {
            population.insert(Uuid::new_v4(), cell);
        }
        Self(population)
    }

    pub fn get_mut_cell(&mut self, id: &StemCellId) -> &mut StemCell {
        self.0.get_mut(id).unwrap()
    }

    pub fn retrieve_all_variants(mut self) -> Vec<Variant> {
        let mut variants = Vec::new();
        for cell in self.0.values_mut() {
            variants.append(&mut cell.private_variants);
            variants.append(&mut cell.shared_variants);
        }
        variants
    }

    pub fn duplicate_cell(&mut self, id: &StemCellId, id2remove: &StemCellId) -> Vec<StemCellId> {
        //! Duplicate a stem cell and remove another stem cellfrom the
        //! population (Moran process).
        //!
        //! ## Panics
        //! Panics when `id` is equal to `id2remove`.
        assert_ne!(id, id2remove, "cannot duplicate and remove the same cell");
        self.remove_cell(id2remove).unwrap();
        let new_cell = self.0[id].clone();
        todo!();
        // self.0[id].children
    }

    pub fn remove_cell(&mut self, id: &StemCellId) -> anyhow::Result<()> {
        let cell = self
            .0
            .remove(id)
            .with_context(|| format!("cannot find cell with id {} in the population", id))?;
        for parent_id in cell.parents.iter() {
            let parent = self.0.get_mut(parent_id).expect(&format!(
                "cannot find parent with id {} in the population",
                id,
            ));
            for variant in cell.shared_variants.iter() {
                parent.shared_variants.push(*variant);
            }
            if let Some(index) = parent
                .children
                .iter_mut()
                .position(|parent_id| parent_id == id)
            {
                parent.children.swap_remove(index);
            }
        }

        for child_id in cell.children.iter() {
            let child = self.0.get_mut(child_id).unwrap();
            if let Some(index) = child
                .parents
                .iter_mut()
                .position(|parent_id| parent_id == id)
            {
                child.parents.swap_remove(index);
            }
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Default)]
pub struct StemCell {
    private_variants: Vec<Variant>,
    shared_variants: Vec<Variant>,
    children: Vec<StemCellId>,
    parents: Vec<StemCellId>,
}

impl StemCell {
    pub fn new(variants: Vec<Variant>) -> Self {
        //! Construct a new stem cell with some neutral mutations
        StemCell {
            private_variants: variants,
            shared_variants: Vec::with_capacity(5000),
            children: Vec::with_capacity(200),
            parents: Vec::with_capacity(200),
        }
    }

    pub fn has_mutations(&self) -> bool {
        //! Returns true if a stem cell has **private** mutations
        !self.private_variants.is_empty()
    }
}

pub struct VariantCount;

impl VariantCount {
    pub fn from_cells(cells: &StemCellPopulation) -> HashMap<Variant, u16> {
        let mut variant_count = HashMap::new();
        for cell in cells.0.values() {
            // private mutations owned by the stem cell
            for variant in cell.private_variants.iter() {
                variant_count
                    .entry(*variant)
                    .and_modify(|counter| *counter += 1)
                    .or_insert(1);
            }

            for variant in cell.shared_variants.iter() {
                variant_count
                    .entry(*variant)
                    .and_modify(|counter| *counter += 1)
                    .or_insert(1);
            }
        }
        variant_count
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arbitrary::{Arbitrary, Result, Unstructured};
    use std::collections::HashSet;

    #[derive(arbitrary::Arbitrary, Clone)]
    pub struct FourCellIdx {
        pub cell0: Uuid,
        pub cell1: Uuid,
        pub cell2: Uuid,
        pub cell3: Uuid,
    }

    #[derive(arbitrary::Arbitrary, Clone)]
    pub struct FiveMutationIdx {
        // shared mut between cell 0 and cell 1
        pub shared0: Uuid,
        // shared mut between cell 0, cell 1 and cell2
        pub shared1: Uuid,
        // private mutation of cell 0
        pub private0: Uuid,
        // private mutation of cell 2
        pub private2: Uuid,
        //Fiveate mutation of cell 3
        pub private3: Uuid,
    }

    #[derive(Clone)]
    struct PopulationTest {
        mutations_idx: FiveMutationIdx,
        cells_idx: FourCellIdx,
        cells: StemCellPopulation,
    }

    impl PopulationTest {
        fn get_unique_mut_idx(&self) -> HashSet<Uuid> {
            HashSet::from([
                self.mutations_idx.private0,
                self.mutations_idx.private2,
                self.mutations_idx.private3,
                self.mutations_idx.shared0,
                self.mutations_idx.shared1,
            ])
        }

        fn get_cell_counts() -> Vec<u16> {
            vec![1, 1, 1, 2, 3]
        }
    }

    impl<'a> arbitrary::Arbitrary<'a> for PopulationTest {
        fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
            let mutations_idx = FiveMutationIdx::arbitrary(u).unwrap();
            let cells_idx = FourCellIdx::arbitrary(u).unwrap();

            let cell3 = StemCell::new(vec![mutations_idx.private3]);

            let mut cell2 = StemCell::new(vec![mutations_idx.private2]);
            cell2.shared_variants = vec![mutations_idx.shared1];
            cell2.parents = vec![cells_idx.cell0, cells_idx.cell1];

            let mut cell1 = StemCell::new(vec![]);
            cell1.shared_variants = vec![mutations_idx.shared0, mutations_idx.shared1];
            cell1.parents = vec![cells_idx.cell0];
            cell1.children = vec![cells_idx.cell2];

            let mut cell0 = StemCell::new(vec![mutations_idx.private0]);
            cell0.shared_variants = vec![mutations_idx.shared0, mutations_idx.shared1];
            cell0.children = vec![cells_idx.cell1, cells_idx.cell2];

            let mut population = FxHashMap::default();
            population.insert(cells_idx.cell0, cell0);
            population.insert(cells_idx.cell1, cell1);
            population.insert(cells_idx.cell2, cell2);
            population.insert(cells_idx.cell3, cell3);

            Ok(PopulationTest {
                mutations_idx,
                cells_idx,
                cells: StemCellPopulation(population),
            })
        }
    }

    #[test]
    fn remove_cell0_test() {
        let raw_data: [u8; 200] = core::array::from_fn(|i| (i + 1) as u8);
        let mut u = Unstructured::new(&raw_data);
        let mut population = PopulationTest::arbitrary(&mut u).unwrap();

        let id2remove = population.cells_idx.cell0;
        population.cells.remove_cell(&id2remove).unwrap();

        // check that cell is not in the population anymore
        assert!(population
            .cells
            .0
            .keys()
            .all(|cell_id| cell_id != &id2remove));
        // check that cell is not in any child
        assert!(population
            .cells
            .0
            .values()
            .all(|cell| cell.children.iter().all(|c| c != &id2remove)));
        // check that cell is not in any parent
        assert!(population
            .cells
            .0
            .values()
            .all(|cell| cell.parents.iter().all(|p| p != &id2remove)));
        // check that the mutations of cell 0 are not there anymore
        let mutations_cell0 = population.mutations_idx.private0;
        dbg!(&population.cells);
        assert!(population
            .cells
            .retrieve_all_variants()
            .into_iter()
            .all(|variant| variant != mutations_cell0));
    }

    #[test]
    fn remove_cell1_test() {
        let raw_data: [u8; 200] = core::array::from_fn(|i| (i + 1) as u8);
        let mut u = Unstructured::new(&raw_data);
        let mut population = PopulationTest::arbitrary(&mut u).unwrap();
        let mut variants = HashSet::with_capacity(5);
        for variant in population.clone().cells.retrieve_all_variants().into_iter() {
            variants.insert(variant);
        }

        let id2remove = population.cells_idx.cell1;
        population.cells.remove_cell(&id2remove).unwrap();

        // check that cell is not in the population anymore
        assert!(population
            .cells
            .0
            .keys()
            .all(|cell_id| cell_id != &id2remove));
        // check that cell is not in any child
        assert!(population
            .cells
            .0
            .values()
            .all(|cell| cell.children.iter().all(|c| c != &id2remove)));
        // check that cell is not in any parent
        assert!(population
            .cells
            .0
            .values()
            .all(|cell| cell.parents.iter().all(|p| p != &id2remove)));

        // check that the parent aka cell 0 has the mutations from cell 1 (
        // which has just been removed)
        let mut expected_muts = HashSet::with_capacity(2);
        expected_muts.insert(&population.mutations_idx.shared0);
        expected_muts.insert(&population.mutations_idx.shared1);

        let mut got_muts = HashSet::with_capacity(3);
        for var in population.cells.0[&population.cells_idx.cell0]
            .shared_variants
            .iter()
        {
            got_muts.insert(var);
        }
        assert_eq!(expected_muts, got_muts);

        // check that the mutations have not changed
        let mut got_variants = HashSet::with_capacity(5);
        for variant in population.cells.retrieve_all_variants().into_iter() {
            got_variants.insert(variant);
        }
        assert_eq!(got_variants, variants);
    }

    #[test]
    fn remove_cell2_test() {
        let raw_data: [u8; 200] = core::array::from_fn(|i| (i + 1) as u8);
        let mut u = Unstructured::new(&raw_data);
        let mut population = PopulationTest::arbitrary(&mut u).unwrap();

        let id2remove = population.cells_idx.cell2;
        population.cells.remove_cell(&id2remove).unwrap();

        // check that cell is not in the population anymore
        assert!(population
            .cells
            .0
            .keys()
            .all(|cell_id| cell_id != &id2remove));
        // check that cell is not in any child
        assert!(population
            .cells
            .0
            .values()
            .all(|cell| cell.children.iter().all(|c| c != &id2remove)));
        // check that cell is not in any parent
        assert!(population
            .cells
            .0
            .values()
            .all(|cell| cell.parents.iter().all(|p| p != &id2remove)));
        // check that the mutations of cell 2 are not there anymore
        let mutations_cell2 = population.mutations_idx.private2;
        dbg!(&population.cells);
        assert!(population
            .cells
            .retrieve_all_variants()
            .into_iter()
            .all(|variant| variant != mutations_cell2));
    }

    #[test]
    fn remove_cell3_test() {
        let raw_data: [u8; 200] = core::array::from_fn(|i| (i + 1) as u8);
        let mut u = Unstructured::new(&raw_data);
        let mut population = PopulationTest::arbitrary(&mut u).unwrap();

        let id2remove = population.cells_idx.cell3;
        population.cells.remove_cell(&id2remove).unwrap();

        // check that cell is not in the population anymore
        assert!(population
            .cells
            .0
            .keys()
            .all(|cell_id| cell_id != &id2remove));
        // check that cell is not in any child
        assert!(population
            .cells
            .0
            .values()
            .all(|cell| cell.children.iter().all(|c| c != &id2remove)));
        // check that cell is not in any parent
        assert!(population
            .cells
            .0
            .values()
            .all(|cell| cell.parents.iter().all(|p| p != &id2remove)));
        // check that the mutations of cell 3 are not there anymore
        let mutations_cell3 = population.mutations_idx.private3;
        assert!(population
            .cells
            .retrieve_all_variants()
            .into_iter()
            .all(|variant| variant != mutations_cell3));
    }

    #[test]
    fn variant_count_from_cells_test() {
        let raw_data: [u8; 200] = core::array::from_fn(|i| (i + 1) as u8);
        let population = PopulationTest::arbitrary(&mut Unstructured::new(&raw_data)).unwrap();
        let expected_mut_idx = population.get_unique_mut_idx();
        let expected_cell_counts = PopulationTest::get_cell_counts();

        let cells = population.clone().cells;
        let variant_count = VariantCount::from_cells(&cells);
        let unique_mut_idx = variant_count.clone().into_keys().collect::<HashSet<Uuid>>();
        let mut cell_counts = variant_count.into_values().collect::<Vec<u16>>();
        cell_counts.sort_unstable();

        assert_eq!(unique_mut_idx, expected_mut_idx);
        assert_eq!(cell_counts, expected_cell_counts);
    }
}
