package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;


public class ConstrainedTreeBranchLengthProvider  implements BranchLengthProvider {
    public static final String CONSTRAINED_TREE_BRANCHLENGTH_PROVIDER = "ConstrainedTreeBranchLengthProvider";
    public ConstrainedTreeBranchLengthProvider(ConstrainedTreeModel constrainedTreeModel,Tree dataTree,Double scale,double minBranchlength,boolean discrete){
        this.scale = scale;
        this.discrete = discrete;
        this.minBranchlength = minBranchlength;


        externalBranchLengths = new double[dataTree.getExternalNodeCount()];
        cladeBranchLengths = new double[dataTree.getInternalNodeCount()];

        Map<String, NodeRef> taxonIdNodeMap = new HashMap<>();
        for (int i = 0; i < dataTree.getExternalNodeCount(); i++) {
            NodeRef node = dataTree.getExternalNode(i);
            taxonIdNodeMap.put(dataTree.getNodeTaxon(node).getId(), node);
        }
        //set up external branchLengths
        for (int i = 0; i < constrainedTreeModel.getExternalNodeCount(); i++) {
            NodeRef node = constrainedTreeModel.getExternalNode(i);
            String taxonId =  constrainedTreeModel.getNodeTaxon(node).getId();
            NodeRef dataNode = taxonIdNodeMap.get(taxonId);

            externalBranchLengths[node.getNumber()] = discrete? Math.round(dataTree.getBranchLength(dataNode)*scale):dataTree.getBranchLength(dataNode)*scale ;
        }
        // need to make subtrees to nodes in the dataTree Subtrees
        Map<BitSet, NodeRef> dataTreeMap = getBitSetNodeMap(dataTree,dataTree);
        Map<BitSet, NodeRef> treeModelMap = getBitSetNodeMap(dataTree,constrainedTreeModel);
        HashMap <NodeRef, NodeRef> dataTreeNodeMap = new HashMap<>();

        for (Map.Entry<BitSet, NodeRef> entry: dataTreeMap.entrySet()){
            dataTreeNodeMap.put(entry.getValue(), treeModelMap.get(entry.getKey()));
        }

        for(int i=0; i<dataTree.getInternalNodeCount(); i++){
            NodeRef dataNode = dataTree.getInternalNode(i);
            NodeRef constrainedNode = dataTreeNodeMap.get(dataNode);
            cladeBranchLengths[constrainedTreeModel.getSubtreeIndex(constrainedNode)] = discrete? Math.round(dataTree.getBranchLength(dataNode)*scale):dataTree.getBranchLength(dataNode)*scale ;
        }
    }
    public ConstrainedTreeBranchLengthProvider(ConstrainedTreeModel constrainedTreeModel,Tree dataTree){
        this(constrainedTreeModel,dataTree,1.0,0.0,true);
    }

    @Override
    public double getBranchLength(Tree tree, NodeRef node) {
        if (tree.isExternal(node)) {
            return externalBranchLengths[node.getNumber()];
        }
        WrappedSubtree subtree = ((ConstrainedTreeModel) tree).getSubtree(node);
        NodeRef nodeInSubtree = subtree.getUnWrappedNode(((ConstrainedTreeModel) tree).getNodeInWrappedTree(node));
        if (subtree.isRoot(nodeInSubtree)) {
            int subtreeIndex = ((ConstrainedTreeModel) tree).getSubtreeIndex(node);
            return cladeBranchLengths[subtreeIndex];
        }else{
            return minBranchlength;
        }
    }


    /**
     * Gets a HashMap of clade bitsets to nodes in tree. This is useful for comparing the topology of trees
     * @param referenceTree  the tree that will be used to define taxa and tip numbers
     * @param tree the tree for which clades are being defined
     * @return A HashMap with a BitSet of descendent taxa as the key and a node as value
     */
    private HashMap<BitSet, NodeRef> getBitSetNodeMap(Tree referenceTree,Tree  tree) {
        HashMap<BitSet, NodeRef> map = new HashMap<>();
        addBits(referenceTree,tree,tree.getRoot(),map);
        return map;
    }

    /**
     *  A private recursive function used by getBitSetNodeMap
     *  This is modeled after the addClades in CladeSet and getClades in compatibility statistic
     * @param referenceTree  the tree that will be used to define taxa and tip numbers
     * @param tree the tree for which clades are being defined
     * @param node current node
     * @param map map that is being appended to
     */
    private BitSet addBits(Tree referenceTree, Tree tree, NodeRef node, HashMap map) {
        BitSet bits = new BitSet();
        if (tree.isExternal(node)) {
            String taxonId = tree.getNodeTaxon(node).getId();
            bits.set(referenceTree.getTaxonIndex(taxonId));

        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef node1 = tree.getChild(node, i);
                bits.or(addBits(referenceTree,tree, node1, map));
            }
        }

        map.put(bits, node);
        return bits;
    }

    private final double[] cladeBranchLengths;
    private final double[] externalBranchLengths;
    private final double scale;
    private final boolean discrete;
    private final double minBranchlength;
}
