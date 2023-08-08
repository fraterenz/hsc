use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

#[derive(Clone, Debug, Default)]
pub struct StemCell {
    id: usize,
    variants: Vec<u16>,
    children: RefCell<Vec<Rc<StemCell>>>,
}

impl StemCell {
    pub fn new() -> Self {
        StemCell::default()
    }

    pub fn has_mutations(&self) -> bool {
        //! Returns true if a stem cell has **private** mutations
        !self.variants.is_empty()
    }
}

pub struct Sfs;

impl Sfs {
    pub fn from_cells(cells: &[Rc<StemCell>]) -> HashMap<u16, u16> {
        let mut sfs = HashMap::new();
        for cell in cells.into_iter() {
            // private mutations owned by the stem cell
            for variant in cell.variants.iter() {
                sfs.entry(*variant)
                    .and_modify(|counter| *counter += 1)
                    .or_insert(1);
            }

            // shared mutations appearing in all children
            for child in cell.children.borrow().iter() {
                for variant in child.variants.iter() {
                    sfs.entry(*variant)
                        .and_modify(|counter| *counter += 1)
                        .or_insert(1);
                }
            }
        }
        sfs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct PopulationTest {}

    impl PopulationTest {
        fn construct_population_three_cells() -> [Rc<StemCell>; 3] {
            // cell1: [3, *cell2], cell2: [4], cell3: [1]
            let cell3 = Rc::new(StemCell {
                id: 3,
                variants: vec![1],
                children: RefCell::new(vec![]),
            });

            let cell2 = Rc::new(StemCell {
                id: 2,
                variants: vec![4],
                children: RefCell::new(vec![]),
            });

            let cell1 = Rc::new(StemCell {
                id: 1,
                variants: vec![3],
                children: RefCell::new(vec![Rc::clone(&cell2)]),
            });

            [cell1, cell2, cell3]
        }

        fn get_sorted_unique_mut_idx() -> [u16; 3] {
            [1, 3, 4]
        }

        fn get_sorted_cell_counts() -> [u16; 3] {
            [1, 1, 2]
        }
    }

    #[test]
    fn sfs_from_cells_test() {
        let cells = PopulationTest::construct_population_three_cells().to_vec();
        let sfs = Sfs::from_cells(&cells);
        let expected_mut_idx = PopulationTest::get_sorted_unique_mut_idx();
        let expected_cell_counts = PopulationTest::get_sorted_cell_counts();

        let mut unique_mut_idx = sfs.clone().into_keys().collect::<Vec<u16>>();
        unique_mut_idx.sort_unstable();
        let mut cell_counts = sfs.into_values().collect::<Vec<u16>>();
        cell_counts.sort_unstable();

        assert_eq!(unique_mut_idx, expected_mut_idx);
        assert_eq!(cell_counts, expected_cell_counts);
    }

    #[test]
    fn sfs_from_cells_without_cell1_test() {
        // removing variant 3 and 4 from cell1
        let mut cells = PopulationTest::construct_population_three_cells().to_vec();
        // rm cell1
        cells.swap_remove(0);
        let sfs = Sfs::from_cells(&cells);
        let mut expected_mut_idx = PopulationTest::get_sorted_unique_mut_idx().to_vec();
        // remove id 3
        expected_mut_idx.swap_remove(1);
        let mut expected_cell_counts = PopulationTest::get_sorted_cell_counts().to_vec();
        // remove count 2
        expected_cell_counts.swap_remove(2);

        let mut unique_mut_idx = sfs.clone().into_keys().collect::<Vec<u16>>();
        unique_mut_idx.sort_unstable();
        let mut cell_counts = sfs.into_values().collect::<Vec<u16>>();
        cell_counts.sort_unstable();

        assert_eq!(unique_mut_idx, expected_mut_idx);
        assert_eq!(cell_counts, expected_cell_counts);
    }

    #[test]
    fn sfs_from_cells_without_cell2_test() {
        // removing variant 4 from cell2
        let mut cells = PopulationTest::construct_population_three_cells().to_vec();
        // rm cell2
        cells.swap_remove(1);
        cells[0].children.replace(vec![]);

        let sfs = Sfs::from_cells(&cells);
        let mut expected_mut_idx = PopulationTest::get_sorted_unique_mut_idx().to_vec();
        // remove id 4
        expected_mut_idx.swap_remove(2);
        let mut expected_cell_counts = PopulationTest::get_sorted_cell_counts().to_vec();
        // remove count 2
        expected_cell_counts.swap_remove(2);

        let mut unique_mut_idx = sfs.clone().into_keys().collect::<Vec<u16>>();
        unique_mut_idx.sort_unstable();
        let mut cell_counts = sfs.into_values().collect::<Vec<u16>>();
        cell_counts.sort_unstable();

        assert_eq!(unique_mut_idx, expected_mut_idx);
        assert_eq!(cell_counts, expected_cell_counts);
    }

    #[test]
    fn sfs_from_cells_without_cell3_test() {
        // removing variant 1 from cell3
        let mut cells = PopulationTest::construct_population_three_cells().to_vec();
        // rm cell3
        cells.swap_remove(2);

        let sfs = Sfs::from_cells(&cells);
        let mut expected_mut_idx = PopulationTest::get_sorted_unique_mut_idx().to_vec();
        // remove id 1
        expected_mut_idx.remove(0);
        let mut expected_cell_counts = PopulationTest::get_sorted_cell_counts().to_vec();
        // remove count 1
        expected_cell_counts.swap_remove(1);

        let mut unique_mut_idx = sfs.clone().into_keys().collect::<Vec<u16>>();
        unique_mut_idx.sort_unstable();
        let mut cell_counts = sfs.into_values().collect::<Vec<u16>>();
        cell_counts.sort_unstable();

        assert_eq!(unique_mut_idx, expected_mut_idx);
        assert_eq!(cell_counts, expected_cell_counts);
    }
}
