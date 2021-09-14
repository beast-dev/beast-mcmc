//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by FernFlower decompiler)
//

package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;

public class FixedTreeBranchLengthProvider implements BranchLengthProvider {
    public static final String FIXED_TREE_BRANCHLENGTH_PROVIDER = "FixedTreeBranchLengthProvider";
    private final double[] branchLengths;
    private final double scale;
    private final boolean discrete;
    private final double minBranchlength;
    private final Tree tree;

    public FixedTreeBranchLengthProvider(Tree fixedTree, Tree dataTree, Double scale, double minBranchlength, boolean discrete) {
        this.scale = scale;
        this.discrete = discrete;
        this.minBranchlength = minBranchlength;
        this.tree = fixedTree;
        this.branchLengths = new double[dataTree.getNodeCount()];
        if (this.tree.getNodeCount() != dataTree.getNodeCount()) {
            throw new RuntimeException("The number of nodes in the datatree does not match that in input tree");
        } else {
            Map<String, NodeRef> taxonIdNodeMap = new HashMap();


            for(int i = 0; i < dataTree.getExternalNodeCount(); ++i) {
                NodeRef node = dataTree.getExternalNode(i);
                taxonIdNodeMap.put(dataTree.getNodeTaxon(node).getId(), node);
            }

            for(int i = 0; i < fixedTree.getExternalNodeCount(); ++i) {
                NodeRef node = fixedTree.getExternalNode(i);
                String taxonId = fixedTree.getNodeTaxon(node).getId();
                NodeRef dataNode = taxonIdNodeMap.get(taxonId);
                this.branchLengths[node.getNumber()] = discrete ? (double)Math.round(dataTree.getBranchLength(dataNode) * scale) : dataTree.getBranchLength(dataNode) * scale;
            }

            Map<BitSet, NodeRef> dataTreeMap = this.getBitSetNodeMap(dataTree, dataTree);
            Map<BitSet, NodeRef> treeModelMap = this.getBitSetNodeMap(dataTree, fixedTree);
            HashMap<NodeRef, NodeRef> dataTreeNodeMap = new HashMap();

            for (Entry <BitSet, NodeRef> entry:
                    dataTreeMap.entrySet()) {
                dataTreeNodeMap.put(entry.getValue(), treeModelMap.get(entry.getKey()));
            }


            for(int i = 0; i < dataTree.getInternalNodeCount(); ++i) {
                NodeRef dataNode = dataTree.getInternalNode(i);
                NodeRef node = dataTreeNodeMap.get(dataNode);
                this.branchLengths[node.getNumber()] = discrete ? (double)Math.round(dataTree.getBranchLength(dataNode) * scale) : dataTree.getBranchLength(dataNode) * scale;
            }

        }
    }

    public FixedTreeBranchLengthProvider(Tree tree, Tree dataTree) {
        this(tree, dataTree, 1.0D, 0.0D, true);
    }

    public double getBranchLength(Tree tree, NodeRef node) {
        if (this.tree == tree) {
            return this.branchLengths[node.getNumber()];
        } else {
            throw new RuntimeException("Unrecognized Tree");
        }
    }

    private HashMap<BitSet, NodeRef> getBitSetNodeMap(Tree referenceTree, Tree tree) {
        HashMap<BitSet, NodeRef> map = new HashMap();
        this.addBits(referenceTree, tree, tree.getRoot(), map);
        return map;
    }

    private BitSet addBits(Tree referenceTree, Tree tree, NodeRef node, HashMap map) {
        BitSet bits = new BitSet();
        if (tree.isExternal(node)) {
            String taxonId = tree.getNodeTaxon(node).getId();
            bits.set(referenceTree.getTaxonIndex(taxonId));
        } else {
            for(int i = 0; i < tree.getChildCount(node); ++i) {
                NodeRef node1 = tree.getChild(node, i);
                bits.or(this.addBits(referenceTree, tree, node1, map));
            }
        }

        map.put(bits, node);
        return bits;
    }
}
