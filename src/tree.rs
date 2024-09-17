use anyhow::{bail, ensure, Context};
use phylotree::tree::{Node, NodeId, Tree};
use std::{
    collections::{HashMap, HashSet},
    path::Path,
};

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

    fn add_cell_from_parent_node(
        &mut self,
        cell_id: StemCellId,
        parent_id: NodeId,
        edge: f64,
        verbosity: u8,
    ) -> anyhow::Result<bool> {
        //! Returns whether the cell was already present in the system as a
        //! leaf node.
        if verbosity > 1 {
            println!("create new node");
        }
        self.last_id += 1;

        let mut node = Node::new();
        node.id = self.last_id;
        if verbosity > 1 {
            println!(
                "create new node with id {}, with parent {} located at {} mutational distance",
                node.id, parent_id, edge
            );
        }
        self.tree.add_child(node, parent_id, Some(edge))?;
        Ok(self.leaves.insert(cell_id, self.last_id).is_some())
        // .with_context(|| "the proliferating cell is not registered as a leaf")?;
    }

    pub fn assign_sibling(
        &mut self,
        cell_id: StemCellId,
        sibling_id: StemCellId,
        edge: f64,
        verbosity: u8,
    ) -> anyhow::Result<bool> {
        //! Returns whether the cell was already present in the system as a
        //! leaf node.
        let sibling = self.tree.get(self.leaves.get(&sibling_id).unwrap())?;
        ensure!(sibling.is_tip());
        let parent_id = sibling.parent.unwrap();
        if verbosity > 1 {
            println!(
                "assign cell with id {} to new node which is sibling of cell with id {}",
                cell_id, sibling_id,
            );
        }
        self.add_cell_from_parent_node(cell_id, parent_id, edge, verbosity)
    }

    pub fn put_cell_on_new_node(
        &mut self,
        cell_id: StemCellId,
        edge: f64,
        verbosity: u8,
    ) -> anyhow::Result<bool> {
        //! Update the tree with a new node and affect the cell with `cell_id`
        //! to it.
        //!
        //! ## Returns
        //! Whether the cell to add on the new node was already present in the
        //! system as a leaf node.
        let parent = self
            .tree
            .get(
                self.leaves
                    .get(&cell_id)
                    .with_context(|| "cannot get cell from leaves")?,
            )
            .with_context(|| "cannot get node from the tree")?;
        ensure!(
            parent.is_tip(),
            format!("parent node {} is not a leaf node", parent)
        );
        self.add_cell_from_parent_node(cell_id, parent.id, edge, verbosity)
    }

    pub fn remove_cell_from_tree(
        &mut self,
        cell_id: StemCellId,
        verbosity: u8,
    ) -> anyhow::Result<()> {
        //! Remove a cell from the tree, do not change the tree's topology as
        //! we still need the node which acts as an ancestral cell.
        if verbosity > 1 {
            println!("removing cell with id {cell_id} from the leaves");
        }
        let removed = self
            .leaves
            .remove(&cell_id)
            .with_context(|| "the cell to remove is not registered as a leaf")?;
        ensure!(self.tree.get(&removed)?.is_tip());
        Ok(())
    }

    pub fn name_leaves_with_cell_idx(&mut self) -> anyhow::Result<()> {
        for (cell_id, node_id) in self.leaves.iter() {
            self.tree.get_mut(node_id)?.name = Some(cell_id.to_string());
        }
        Ok(())
    }

    pub fn create_tree_without_dead_cells(
        self,
        cells: &[&StemCell],
        verbosity: u8,
    ) -> anyhow::Result<Self> {
        //! Use `cells` to select only the leaves of interest, dropping all the
        //! other leaves and consuming `self`.
        // select leaves in self.leaves that are in `cells`
        let mut leaves = self.leaves;
        leaves.retain(|&k, _| cells.iter().any(|x| x.id == k));
        if verbosity > 0 {
            println!("creating new tree with {} cells", cells.len());
        }
        // get the node idx of those subset of leaves
        let leaves_idx: HashSet<&NodeId> = leaves.values().collect();

        let mut pruned = self.tree.clone();
        while pruned
            .get_leaves()
            .iter()
            .any(|leaf| !leaves_idx.contains(leaf))
        {
            for leaf in pruned.get_leaves() {
                if !leaves_idx.contains(&leaf) {
                    pruned
                        .prune(&leaf)
                        .with_context(|| format!("cannot remove node {}", leaf))?;
                }
            }
        }

        if verbosity > 0 {
            println!("compressing the tree")
        }
        pruned.compress().with_context(|| "cannot compress")?;
        if verbosity > 0 {
            println!("end tree compression")
        }
        // pruned.rescale(
        //     1. / pruned
        //         .length()
        //         .with_context(|| "cannot get lenght of the tree")?,
        // );
        Ok(PhyloTree {
            tree: pruned,
            last_id: self.last_id,
            leaves,
        })
    }

    pub fn get_mut_leaf(&mut self, cell_id: StemCellId) -> anyhow::Result<&mut Node> {
        let node = self
            .tree
            .get_mut(
                self.leaves
                    .get_mut(&cell_id)
                    .with_context(|| "cell not registered as leaf node")?,
            )
            .with_context(|| "cannot find leaf tree node")?;
        if !node.is_tip() {
            bail!("node found is not a leaf node")
        }
        Ok(node)
    }
}

pub fn save_tree(
    tree: PhyloTree,
    cells: &[&StemCell],
    path: &Path,
    verbosity: u8,
) -> anyhow::Result<()> {
    //! Consumes the tree as we need to create a new tree with only some `cells`
    if verbosity > 0 {
        println!("saving the tree");
    }
    let mut pruned = tree.create_tree_without_dead_cells(cells, verbosity)?;
    pruned.name_leaves_with_cell_idx()?;
    // let pruned = tree.clone();
    if verbosity > 1 {
        println!("Name leaf nodes with idx of cells");
    }
    pruned.tree.to_file(path)?;
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
    fn put_cell_on_new_node_test(id: usize, edge: f64) -> bool {
        let parent = StemCell::with_mutations(vec![Uuid::new_v4(), Uuid::new_v4()], id);
        let mut tree = PhyloTree::with_cell(parent.id);
        tree.put_cell_on_new_node(id, edge, 0).is_ok()
    }

    #[test]
    #[should_panic]
    fn put_cell_on_new_node_wrong_id_test() {
        let parent = StemCell::with_mutations(vec![Uuid::new_v4(), Uuid::new_v4()], 0);
        let mut tree = PhyloTree::with_cell(parent.id);
        tree.put_cell_on_new_node(4, 3., 0).unwrap();
    }

    #[test]
    #[should_panic]
    fn assign_sibling_already_there_test() {
        let parent = StemCell::with_mutations(vec![Uuid::new_v4(), Uuid::new_v4()], 0);
        let mut tree = PhyloTree::with_cell(parent.id);
        tree.assign_sibling(parent.id, parent.id, 3., 0).unwrap();
    }

    #[quickcheck]
    fn assign_sibling_test(idx: TwoIDs, edge: f64) -> bool {
        let (id1, id2) = (idx.0 .0, idx.0 .1);
        let stem_cell = StemCell::with_mutations(vec![Uuid::new_v4(), Uuid::new_v4()], id1);
        let sibling = StemCell::with_mutations(vec![Uuid::new_v4()], id2);
        let mut tree = PhyloTree::with_cell(sibling.id);
        tree.put_cell_on_new_node(sibling.id, edge, 0).unwrap();

        tree.assign_sibling(stem_cell.id, sibling.id, edge, 0)
            .unwrap();
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
        tree.remove_cell_from_tree(id, 0).is_ok()
    }

    #[test]
    #[should_panic]
    fn remove_cell_from_tree_wrong_id_test() {
        let cell = StemCell::with_mutations(vec![Uuid::new_v4(), Uuid::new_v4()], 0);
        let mut tree = PhyloTree::with_cell(cell.id);
        tree.remove_cell_from_tree(4, 0).unwrap();
    }

    #[quickcheck]
    fn proliferate_cell_with_id_test(idx: TwoIDs, edge: f64) -> bool {
        let (id1, id2) = (idx.0 .0, idx.0 .1);
        let mut parent = StemCell::new(id1);
        let mut tree = PhyloTree::with_cell(parent.id);

        // proliferation
        let mut daughter = parent.clone();
        daughter.id = id2;
        daughter.variants.push(Uuid::new_v4());
        parent.variants.push(Uuid::new_v4());

        tree.put_cell_on_new_node(parent.id, edge, 0).unwrap();
        tree.assign_sibling(daughter.id, parent.id, edge, 0)
            .unwrap();
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

    #[derive(Clone, Debug)]
    struct TreeTest {
        phylo: PhyloTree,
        cell_idx: Vec<usize>,
    }
    impl TreeTest {
        fn new() -> Self {
            let original = "((A:1.2,(C:1.2,E:0.1)D:0.05)B:1.2,(H:1.1,(Z:2)I:3)G:2)F:0;";
            let phylotree = Tree::from_newick(original).unwrap();
            let cell_idx = vec![2, 40, 10, 1, 3];
            let node_idx = vec![2, 4, 5, 7, 9];
            assert_eq!(phylotree.get_leaves(), node_idx);
            let leaves: HashMap<StemCellId, NodeId> =
                cell_idx.clone().into_iter().zip(node_idx).collect();

            TreeTest {
                phylo: PhyloTree {
                    tree: phylotree,
                    last_id: 0,
                    leaves,
                },
                cell_idx,
            }
        }
    }

    fn get_leaves_4test(
        tree1: TreeTest,
        tree2: &PhyloTree,
        cells: Vec<&StemCell>,
    ) -> (Vec<NodeId>, Vec<NodeId>) {
        let mut leaves_exp = cells
            .into_iter()
            .map(|cell| *tree1.phylo.leaves.get(&cell.id).unwrap())
            .collect::<Vec<usize>>();
        leaves_exp.sort_unstable();
        let mut leaves = tree2.tree.get_leaves();
        leaves.sort_unstable();
        (leaves, leaves_exp)
    }

    #[test]
    fn create_tree_without_dead_cells_keep_all() {
        let phylo = TreeTest::new();
        let cells2keep: Vec<StemCell> = phylo
            .clone()
            .cell_idx
            .into_iter()
            .map(StemCell::new)
            .collect();
        let mut cells2keep_ref: Vec<&StemCell> = Vec::new();
        for cell in cells2keep.iter() {
            cells2keep_ref.push(cell);
        }
        let pruned = phylo
            .clone()
            .phylo
            .create_tree_without_dead_cells(&cells2keep_ref, 0)
            .unwrap();

        let mut expected = phylo.phylo.tree.clone();
        expected.compress().unwrap();
        assert_eq!(
            pruned.tree.to_newick().unwrap(),
            expected.to_newick().unwrap()
        );
        let (leaves, leaves_exp) = get_leaves_4test(phylo, &pruned, cells2keep_ref);
        assert_eq!(leaves, leaves_exp);
    }

    #[test]
    fn create_tree_without_dead_cells_remove_h_all() {
        let phylo = TreeTest::new();
        let cells2keep: Vec<StemCell> = phylo
            .clone()
            .cell_idx
            .into_iter()
            .map(StemCell::new)
            .collect();
        let mut cells2keep_ref: Vec<&StemCell> = Vec::new();
        for cell in cells2keep.iter() {
            // remove H
            if cell.id != 1 {
                cells2keep_ref.push(cell);
            }
        }
        let pruned = phylo
            .clone()
            .phylo
            .create_tree_without_dead_cells(&cells2keep_ref, 0)
            .unwrap();

        assert_eq!(
            pruned
                .tree
                .to_formatted_newick(phylotree::tree::NewickFormat::OnlyNames)
                .unwrap(),
            "((A,(C,E)D)B,Z)F;"
        );
        let (leaves, leaves_exp) = get_leaves_4test(phylo, &pruned, cells2keep_ref);
        assert_eq!(leaves, leaves_exp);
    }

    #[test]
    fn create_tree_without_dead_cells_remove_c_all() {
        let phylo = TreeTest::new();
        let cells2keep: Vec<StemCell> = phylo
            .clone()
            .cell_idx
            .into_iter()
            .map(StemCell::new)
            .collect();
        let mut cells2keep_ref: Vec<&StemCell> = Vec::new();
        for cell in cells2keep.iter() {
            // remove c
            if cell.id != 40 {
                cells2keep_ref.push(cell);
            }
        }
        let pruned = phylo
            .clone()
            .phylo
            .create_tree_without_dead_cells(&cells2keep_ref, 0)
            .unwrap();

        assert_eq!(
            pruned
                .tree
                .to_formatted_newick(phylotree::tree::NewickFormat::OnlyNames)
                .unwrap(),
            "((A,E)B,(H,Z)G)F;"
        );
        let (leaves, leaves_exp) = get_leaves_4test(phylo, &pruned, cells2keep_ref);
        assert_eq!(leaves, leaves_exp);
    }

    #[test]
    fn create_tree_without_dead_cells_remove_a_all() {
        let phylo = TreeTest::new();
        let cells2keep: Vec<StemCell> = phylo
            .clone()
            .cell_idx
            .into_iter()
            .map(StemCell::new)
            .collect();
        let mut cells2keep_ref: Vec<&StemCell> = Vec::new();
        for cell in cells2keep.iter() {
            // remove a
            if cell.id != 2 {
                cells2keep_ref.push(cell);
            }
        }
        let pruned = phylo
            .clone()
            .phylo
            .create_tree_without_dead_cells(&cells2keep_ref, 0)
            .unwrap();

        assert_eq!(
            pruned
                .tree
                .to_formatted_newick(phylotree::tree::NewickFormat::OnlyNames)
                .unwrap(),
            "((H,Z)G,(C,E)D)F;"
        );
        let (leaves, leaves_exp) = get_leaves_4test(phylo, &pruned, cells2keep_ref);
        assert_eq!(leaves, leaves_exp);
    }

    #[test]
    fn create_tree_without_dead_cells_remove_a_and_c() {
        let phylo = TreeTest::new();
        let cells2keep: Vec<StemCell> = phylo
            .clone()
            .cell_idx
            .into_iter()
            .map(StemCell::new)
            .collect();
        let mut cells2keep_ref: Vec<&StemCell> = Vec::new();
        for cell in cells2keep.iter() {
            // remove a and c
            if cell.id != 2 && cell.id != 40 {
                cells2keep_ref.push(cell);
            }
        }
        let pruned = phylo
            .clone()
            .phylo
            .create_tree_without_dead_cells(&cells2keep_ref, 0)
            .unwrap();

        assert_eq!(
            pruned
                .tree
                .to_formatted_newick(phylotree::tree::NewickFormat::OnlyNames)
                .unwrap(),
            "((H,Z)G,E)F;"
        );
        let (leaves, leaves_exp) = get_leaves_4test(phylo, &pruned, cells2keep_ref);
        assert_eq!(leaves, leaves_exp);
    }
}
