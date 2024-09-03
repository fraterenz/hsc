use anyhow::{ensure, Context};
use phylotree::tree::{Node, NodeId, Tree};
use std::collections::HashMap;

use crate::stemcell::{StemCell, StemCellId};

#[derive(Debug, Clone)]
pub struct PhyloTree {
    pub tree: Tree,
    last_id: NodeId,
    /// this places stem cells on the tree
    leaves: HashMap<StemCellId, NodeId>,
}

impl Default for PhyloTree {
    fn default() -> Self {
        Self::new()
    }
}

impl PhyloTree {
    pub fn new() -> Self {
        let mut tree = Tree::new();
        let last_id = tree.add(Node::new());

        Self {
            tree,
            last_id,
            leaves: HashMap::new(),
        }
    }

    pub fn with_cell(cell_id: StemCellId) -> Self {
        let mut tree = Tree::new();
        let last_id = tree.add(Node::new());
        // assert!(tree.is_rooted().unwrap());
        let leaves = HashMap::from([(cell_id, last_id)]);

        Self {
            tree,
            last_id,
            leaves,
        }
    }

    pub fn assign_sibling(
        &mut self,
        cell_id: StemCellId,
        sibling_id: StemCellId,
    ) -> anyhow::Result<()> {
        let parent_id = self
            .tree
            .get(self.leaves.get(&sibling_id).unwrap())?
            .parent
            .unwrap();
        self.last_id += 1;

        let mut node = Node::new();
        node.id = self.last_id;
        self.tree.add_child(node, parent_id, None)?;
        ensure!(self.leaves.insert(cell_id, self.last_id).is_none());
        Ok(())
    }

    pub fn put_cell_on_new_node(&mut self, cell_id: StemCellId) -> anyhow::Result<()> {
        //! Update the tree with a new node and affect the cell with `cell_id`
        //! to it.
        // create new node
        let parent_id = self.tree.get(self.leaves.get(&cell_id).unwrap())?.id;
        self.last_id += 1;
        let mut node = Node::new();
        node.id = self.last_id;
        self.tree.add_child(node, parent_id, None)?;

        // update leaves
        self.leaves
            .insert(cell_id, self.last_id)
            .with_context(|| "the proliferating cell is not registered as a leaf")?;
        Ok(())
    }

    pub fn remove_cell_from_tree(&mut self, cell_id: StemCellId) -> anyhow::Result<()> {
        //! Remove a cell from the tree, do not change the tree's topology as
        //! we still the node where the cell was for the other alive lineages.
        self.leaves
            .remove(&cell_id)
            .with_context(|| "the cell to remove is not registered as a leaf")?;
        Ok(())
    }

    pub fn create_tree_without_dead_cells(
        &self,
        cells: &[&StemCell],
        verbosity: u8,
    ) -> anyhow::Result<Tree> {
        //! Use `cells` to select only the leaves of interest, dropping all the
        //! other leaves.
        let mut cell2prune = self.leaves.clone();
        cell2prune.retain(|&k, _| !cells.iter().any(|x| x.id == k));
        let mut pruned = self.tree.clone();
        for node_id in cell2prune.values() {
            if verbosity > 1 {
                println!("removing node with id {}", node_id);
            }
            pruned.prune(node_id)?;
        }
        for (cell_id, node_id) in self.leaves.iter() {
            pruned.get_mut(node_id)?.name = Some(cell_id.to_string());
        }
        pruned.compress()?;
        Ok(pruned)
    }
}

pub fn proliferate_cell_with_id(
    parent: &StemCell,
    daughter: &StemCell,
    tree: &mut PhyloTree,
) -> anyhow::Result<()> {
    tree.put_cell_on_new_node(parent.id)?;
    tree.assign_sibling(daughter.id, parent.id)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck_macros::quickcheck;
    use uuid::Uuid;

    // Both `Clone` and `Debug` are required by `quickcheck`
    #[derive(Debug, Clone)]
    struct TwoIDs((usize, usize));

    impl quickcheck::Arbitrary for TwoIDs {
        fn arbitrary(g: &mut quickcheck::Gen) -> Self {
            let id1 = *g.choose(&[usize::MIN, usize::MAX]).unwrap();
            let mut id2 = *g.choose(&[usize::MIN, usize::MAX]).unwrap();
            while id2 == id1 {
                id2 = *g.choose(&[usize::MIN, usize::MAX]).unwrap();
            }
            TwoIDs((id1, id2))
        }
    }

    #[quickcheck]
    fn with_cell_test(id: usize) -> bool {
        let parent = StemCell::with_mutations(vec![Uuid::new_v4()], id);
        let tree = PhyloTree::with_cell(parent.id);
        tree.leaves.contains_key(&id) && tree.leaves.len() == 1
    }

    #[quickcheck]
    fn put_cell_on_new_node_test(id: usize) -> bool {
        let parent = StemCell::with_mutations(vec![Uuid::new_v4(), Uuid::new_v4()], id);
        let mut tree = PhyloTree::with_cell(parent.id);
        tree.put_cell_on_new_node(id).is_ok()
    }

    #[test]
    #[should_panic]
    fn put_cell_on_new_node_wrong_id_test() {
        let parent = StemCell::with_mutations(vec![Uuid::new_v4(), Uuid::new_v4()], 0);
        let mut tree = PhyloTree::with_cell(parent.id);
        tree.put_cell_on_new_node(4).unwrap();
    }

    #[test]
    #[should_panic]
    fn assign_sibling_already_there_test() {
        let parent = StemCell::with_mutations(vec![Uuid::new_v4(), Uuid::new_v4()], 0);
        let mut tree = PhyloTree::with_cell(parent.id);
        tree.assign_sibling(parent.id, parent.id).unwrap();
    }

    #[quickcheck]
    fn assign_sibling_test(idx: TwoIDs) -> bool {
        let (id1, id2) = (idx.0 .0, idx.0 .1);
        let stem_cell = StemCell::with_mutations(vec![Uuid::new_v4(), Uuid::new_v4()], id1);
        let sibling = StemCell::with_mutations(vec![Uuid::new_v4()], id2);
        let mut tree = PhyloTree::with_cell(sibling.id);
        tree.put_cell_on_new_node(sibling.id).unwrap();

        tree.assign_sibling(stem_cell.id, sibling.id).unwrap();
        let node1 = tree.tree.get(tree.leaves.get(&id1).unwrap()).unwrap();
        let parent1 = node1.parent.unwrap();
        let node2 = tree.tree.get(tree.leaves.get(&id2).unwrap()).unwrap();
        let parent2 = node2.parent.unwrap();
        parent1 == parent2
            && node1.id != node2.id
            && tree.tree.get_leaves().len() == 2
            && tree.tree.is_rooted().is_ok()
            && tree.tree.is_binary().is_ok()
    }

    #[quickcheck]
    fn remove_cell_from_tree_test(id: usize) -> bool {
        let cell = StemCell::with_mutations(vec![Uuid::new_v4(), Uuid::new_v4()], id);
        let mut tree = PhyloTree::with_cell(cell.id);
        tree.remove_cell_from_tree(id).is_ok()
    }

    #[test]
    #[should_panic]
    fn remove_cell_from_tree_wrong_id_test() {
        let cell = StemCell::with_mutations(vec![Uuid::new_v4(), Uuid::new_v4()], 0);
        let mut tree = PhyloTree::with_cell(cell.id);
        tree.remove_cell_from_tree(4).unwrap();
    }

    #[quickcheck]
    fn proliferate_cell_with_id_test(idx: TwoIDs) -> bool {
        let (id1, id2) = (idx.0 .0, idx.0 .1);
        let mut parent = StemCell::new(id1);
        let mut tree = PhyloTree::with_cell(parent.id);

        // proliferation
        let mut daughter = parent.clone();
        daughter.id = id2;
        daughter.variants.push(Uuid::new_v4());
        parent.variants.push(Uuid::new_v4());

        proliferate_cell_with_id(&parent, &daughter, &mut tree).unwrap();
        let mut leaves = tree.tree.get_leaves();
        leaves.sort_unstable();
        let exp_leaves = [1, 2];
        leaves
            .iter()
            .enumerate()
            .all(|(i, &ele)| ele == exp_leaves[i])
            && tree.leaves.len() == 2
            && tree.leaves.contains_key(&id1)
            && tree.leaves.contains_key(&(id2))
            && tree.tree.is_rooted().is_ok()
            && tree.tree.is_binary().is_ok()
    }

    #[test]
    fn create_tree_without_dead_cells_keep_all() {
        let original = "((A,(C,E)D)B,(H,Z)G)F;";
        let phylotree = Tree::from_newick(original).unwrap();
        let cell_idx = vec![2, 40, 10, 1, 2];
        let node_idx = vec![2, 4, 5, 7, 8];
        assert_eq!(phylotree.get_leaves(), node_idx);
        let leaves: HashMap<StemCellId, NodeId> =
            cell_idx.clone().into_iter().zip(node_idx).collect();

        let tree = PhyloTree {
            tree: phylotree,
            last_id: 0,
            leaves,
        };

        let cells2keep: Vec<StemCell> = cell_idx.into_iter().map(StemCell::new).collect();
        let mut cells2keep_ref: Vec<&StemCell> = Vec::new();
        for cell in cells2keep.iter() {
            cells2keep_ref.push(cell);
        }
        let pruned = tree
            .create_tree_without_dead_cells(&cells2keep_ref, 0)
            .unwrap();

        assert_eq!(pruned.to_newick().unwrap(), tree.tree.to_newick().unwrap());
    }

    #[test]
    fn create_tree_without_dead_cells_keep_all_collapsed() {
        let original = "((A,(C,E)D)B,((H),(Z)I)G)F;";
        let phylotree = Tree::from_newick(original).unwrap();
        let cell_idx = vec![2, 40, 10, 1, 2];
        let node_idx = vec![2, 4, 5, 8, 10];
        assert_eq!(phylotree.get_leaves(), node_idx);
        let leaves: HashMap<StemCellId, NodeId> =
            cell_idx.clone().into_iter().zip(node_idx).collect();

        let tree = PhyloTree {
            tree: phylotree,
            last_id: 0,
            leaves,
        };

        let cells2keep: Vec<StemCell> = cell_idx.into_iter().map(StemCell::new).collect();
        let mut cells2keep_ref: Vec<&StemCell> = Vec::new();
        for cell in cells2keep.iter() {
            cells2keep_ref.push(cell);
        }
        let pruned = tree
            .create_tree_without_dead_cells(&cells2keep_ref, 0)
            .unwrap();

        assert_eq!(pruned.to_newick().unwrap(), "((A,(C,E)D)B,(H,Z)G)F;",);
    }

    #[test]
    fn create_tree_without_dead_cells_remove_h_all() {
        let original = "((A,(C,E)D)B,(H,Z)G)F;";
        let phylotree = Tree::from_newick(original).unwrap();
        let cell_idx = vec![2, 40, 10, 1];
        let node_idx = vec![2, 4, 5, 7, 8];
        assert_eq!(phylotree.get_leaves(), node_idx);
        let leaves: HashMap<StemCellId, NodeId> =
            cell_idx.clone().into_iter().zip(node_idx).collect();

        let tree = PhyloTree {
            tree: phylotree,
            last_id: 0,
            leaves,
        };

        let cells2keep: Vec<StemCell> = cell_idx.into_iter().map(StemCell::new).collect();
        let mut cells2keep_ref: Vec<&StemCell> = Vec::new();
        for cell in cells2keep.iter() {
            // remove H
            if cell.id != 1 {
                cells2keep_ref.push(cell);
            }
        }
        let pruned = tree
            .create_tree_without_dead_cells(&cells2keep_ref, 0)
            .unwrap();

        assert_eq!(pruned.to_newick().unwrap(), "((A,(C,E)D)B,Z)F;");
    }

    #[test]
    fn create_tree_without_dead_cells_remove_c_all() {
        let original = "((A,(C,E)D)B,(H,Z)G)F;";
        let phylotree = Tree::from_newick(original).unwrap();
        let cell_idx = vec![2, 40, 10, 1];
        let node_idx = vec![2, 4, 5, 7, 8];
        assert_eq!(phylotree.get_leaves(), node_idx);
        let leaves: HashMap<StemCellId, NodeId> =
            cell_idx.clone().into_iter().zip(node_idx).collect();

        let tree = PhyloTree {
            tree: phylotree,
            last_id: 0,
            leaves,
        };

        let cells2keep: Vec<StemCell> = cell_idx.into_iter().map(StemCell::new).collect();
        let mut cells2keep_ref: Vec<&StemCell> = Vec::new();
        for cell in cells2keep.iter() {
            // remove c
            if cell.id != 40 {
                cells2keep_ref.push(cell);
            }
        }
        let pruned = tree
            .create_tree_without_dead_cells(&cells2keep_ref, 0)
            .unwrap();

        assert_eq!(pruned.to_newick().unwrap(), "((A,E)B,(H,Z)G)F;");
    }

    #[test]
    fn create_tree_without_dead_cells_remove_a_all() {
        let original = "((A,(C,E)D)B,(H,Z)G)F;";
        let phylotree = Tree::from_newick(original).unwrap();
        let cell_idx = vec![2, 40, 10, 1];
        let node_idx = vec![2, 4, 5, 7, 8];
        assert_eq!(phylotree.get_leaves(), node_idx);
        let leaves: HashMap<StemCellId, NodeId> =
            cell_idx.clone().into_iter().zip(node_idx).collect();

        let tree = PhyloTree {
            tree: phylotree,
            last_id: 0,
            leaves,
        };

        let cells2keep: Vec<StemCell> = cell_idx.into_iter().map(StemCell::new).collect();
        let mut cells2keep_ref: Vec<&StemCell> = Vec::new();
        for cell in cells2keep.iter() {
            // remove a
            if cell.id != 2 {
                cells2keep_ref.push(cell);
            }
        }
        let pruned = tree
            .create_tree_without_dead_cells(&cells2keep_ref, 0)
            .unwrap();

        assert_eq!(pruned.to_newick().unwrap(), "((H,Z)G,(C,E)D)F;");
    }
}
