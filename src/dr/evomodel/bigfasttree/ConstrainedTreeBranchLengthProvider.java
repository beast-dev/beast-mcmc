package dr.evomodel.bigfasttree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;

public class ConstrainedTreeBranchLengthProvider implements BranchLengthProvider{

    public ConstrainedTreeBranchLengthProvider(Tree dataTree, Tree referenceTree) {

        // An array where the entry points to the node in the data tree that maps to the node in the time tree
        // Nodes that result from resolving a polytomy in the data tree will point to the polytomy node.
        this.nodeInDataTree = new int[referenceTree.getNodeCount()];

        HashMap<BitSet, NodeRef> dataTreeMap = getBitSetNodeMap(dataTree,dataTree);

        HashMap<BitSet, NodeRef> treeModelMap = getBitSetNodeMap(dataTree,referenceTree);
        // reverse map to be <NodeRef, BitSet>
        HashMap <NodeRef, BitSet> treeModelNodeMap = new HashMap<>();
        for (Map.Entry<BitSet, NodeRef> entry: treeModelMap.entrySet()){
            treeModelNodeMap.put(entry.getValue(), entry.getKey());
        }


        setUpNodeMap(treeModelNodeMap,dataTreeMap,referenceTree.getRoot(),referenceTree,null);

        branchLengths = new double[dataTree.getNodeCount()];

        for (int i = 0; i < dataTree.getNodeCount(); i++) {
            double x = dataTree.getBranchLength(dataTree.getNode(i));
            branchLengths[i] = x;
        }

    }

    /**
     * A private recursive method that sets the map from the treemodel to the data tree. If the data tree is resolved it
     * will be a 1 to 1 map. If there are polytomies in the data tree all inserted nodes will map to the polytomy node.
     * @param treeModelNodeMap
     * @param dataTreeMap
     * @param node
     * @param parentsClade
     */
    private void setUpNodeMap(HashMap <NodeRef, BitSet> treeModelNodeMap, HashMap<BitSet, NodeRef> dataTreeMap,NodeRef node, Tree referenceTree,BitSet parentsClade){
        BitSet clade = treeModelNodeMap.get(node);
        int j = node.getNumber();
        if(!dataTreeMap.containsKey(clade)){
            clade=parentsClade;
        }
        nodeInDataTree[j] = dataTreeMap.get(clade).getNumber();
        for (int i = 0; i <referenceTree.getChildCount(node) ; i++) {
            NodeRef child = referenceTree.getChild(node,i);
            setUpNodeMap(treeModelNodeMap,dataTreeMap,child,referenceTree,clade);
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


    private final int[] nodeInDataTree;
    private final double[] branchLengths;


    @Override
    public double getBranchLength(Tree tree, NodeRef node) {
            int dataNode = nodeInDataTree[node.getNumber()];
            int dataParentNode = nodeInDataTree[tree.getParent(tree.getNode(node.getNumber())).getNumber()];
            if(dataNode==dataParentNode){
                return 0d;
            }else{
                return branchLengths[dataNode];
            }
    }
}
