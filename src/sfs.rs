use ssa::ReactionRates;
use std::collections::HashMap;
use std::rc::Rc;
use uuid::Uuid;

pub type CloneId = Uuid;
/// We assume that each division generates a different set of neutral
/// mutations (inifinte-site assumption?).
pub type NeutralMutation = Uuid;

#[derive(Debug, Clone, Copy)]
pub struct SubClone {
    pub id: CloneId,
    /// The fitness coefficient
    pub s: f32,
}

#[derive(Debug, Clone)]
pub struct StemCell {
    pub id: u64,
    pub subclone: Rc<SubClone>,
    pub mutations: Vec<NeutralMutation>,
}

pub fn assign_cell(cell: &mut StemCell, clone: &Rc<SubClone>) {
    cell.subclone = Rc::clone(clone);
}

pub fn compute_reactions<const MAX_SUBCLONES: usize>(
    cells: &[StemCell],
) -> ReactionRates<MAX_SUBCLONES> {
    // sfs(cells)
    todo!();
}

pub fn sfs(cells: &[StemCell]) -> HashMap<CloneId, u16> {
    //! Compute the site frequency spectrum for clones that have a
    //! proliferative advantage.
    let mut clones = HashMap::<CloneId, u16>::new();
    for cell in cells {
        clones
            .entry(cell.subclone.id)
            .and_modify(|count| *count += 1)
            .or_insert(1);
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
    use quickcheck::{Arbitrary, Gen};
    use quickcheck_macros::quickcheck;
    use rand::SeedableRng;
    use std::num::NonZeroU8;

    #[test]
    fn assign_cell_test() {
        let clone = SubClone {
            id: CloneId::new_v4(),
            s: 1.9,
        };
        let neutral_clone = SubClone {
            id: CloneId::new_v4(),
            s: 0.,
        };
        let mut cell = StemCell {
            id: 0,
            subclone: Rc::new(neutral_clone),
            mutations: Vec::new(),
        };
        assert_eq!(neutral_clone.id, cell.subclone.id);
        assert_ne!(clone.id, cell.subclone.id);

        assign_cell(&mut cell, &Rc::new(clone));
        assert_eq!(clone.id, cell.subclone.id);
        assert_ne!(neutral_clone.id, cell.subclone.id);
    }

    #[quickcheck]
    fn sfs_test(s0: f32, s1: f32) -> bool {
        let clone0 = SubClone {
            id: CloneId::new_v4(),
            s: s0,
        };
        let clone1 = SubClone {
            id: CloneId::new_v4(),
            s: s1,
        };

        let cell = StemCell {
            id: 0,
            subclone: Rc::new(clone0),
            mutations: Vec::new(),
        };
        let cell_1 = StemCell {
            id: 1,
            subclone: Rc::new(clone1),
            mutations: Vec::new(),
        };

        let sfs = sfs(&[cell, cell_1]);
        sfs[&clone0.id] == 1u16 && sfs[&clone1.id] == 1
    }

    #[quickcheck]
    fn sfs_after_assignment_test(s0: f32, s1: f32) -> bool {
        let clone0 = SubClone {
            id: CloneId::new_v4(),
            s: s0,
        };
        let clone1 = SubClone {
            id: CloneId::new_v4(),
            s: s1,
        };

        let cell = StemCell {
            id: 0,
            subclone: Rc::new(clone0),
            mutations: Vec::new(),
        };
        let mut cell_1 = StemCell {
            id: 1,
            subclone: Rc::new(clone1),
            mutations: Vec::new(),
        };

        assign_cell(&mut cell_1, &Rc::new(clone0));

        let sfs = sfs(&[cell, cell_1]);
        sfs[&clone0.id] == 2 && sfs.get(&clone1.id) == None
    }
}
