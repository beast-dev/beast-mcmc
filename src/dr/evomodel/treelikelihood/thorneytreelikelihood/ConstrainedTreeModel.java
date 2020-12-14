package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.*;
import dr.evomodel.tree.TreeModel;

import java.util.*;

/**
 * A tree model that only allows node height operations an/or operations that maintain constrained clades
 */
public class ConstrainedTreeModel extends AbstractWrappedTree {
    //
    // Public stuff
    //

    public static final String CONSTRAINED_TREE_MODEL = "constrainedTreeModel";


    public ConstrainedTreeModel(String name, TreeModel tree, List<WrappedSubtree> subtrees){
        super(name,tree);
        this.subtrees = subtrees;
        nodeToSubtreeMap = new int[getNodeCount()];
        Arrays.fill(nodeToSubtreeMap,-1);

        int visitedNodes = 0;
        for (int i=0; i<subtrees.size();i++) {
            WrappedSubtree subtree = subtrees.get(i);
            for (int j = 0; i < subtree.getNodeCount(); j++) {

                NodeRef node = subtree.getNode(j);
                NodeRef nodeInThisTree = convertSubtreeNodeToConstrainedNode(subtree,node);
        // subtree roots should be external nodes in parent clades. They are stored in the tree where they are the root
                if(subtree.isExternal(node)&&!isExternal(nodeInThisTree)){
                    continue;
                }
                if(nodeToSubtreeMap[nodeInThisTree.getNumber()]!=-1){
                    throw new IllegalArgumentException("The subtrees in a constrained tree must be monophyletic");
                }
                nodeToSubtreeMap[nodeInThisTree.getNumber()]=i;
                visitedNodes++;
            }
        }
        if(visitedNodes!=getNodeCount()){
            throw new IllegalArgumentException("The subtrees in a constrained tree must include all nodes in the tree");
        }
    }
    public ConstrainedTreeModel(TreeModel tree, Tree constraintsTree){
        this(CONSTRAINED_TREE_MODEL, tree, constraintsTree);
    }
    public ConstrainedTreeModel(String name, TreeModel tree, Tree constraintsTree) {
        super(name, tree);

        // need to make subtrees to nodes in the dataTree Subtrees
        Map<BitSet, NodeRef> constraintsTreeMap = getBitSetNodeMap(constraintsTree,constraintsTree);
        Map<BitSet, NodeRef> treeModelMap = getBitSetNodeMap(constraintsTree,tree);

        Map<NodeRef, NodeRef> constraintsNodeToTreeNode = new HashMap<>();
        for (BitSet clade :
                constraintsTreeMap.keySet()) {
            NodeRef constraintsTreeNode = constraintsTreeMap.get(clade);
            NodeRef treeModelNode = treeModelMap.get(clade);
            constraintsNodeToTreeNode.put(constraintsTreeNode, treeModelNode);
        }


        this.subtrees = new ArrayList<>();
        getSubtrees(constraintsTree, constraintsNodeToTreeNode, constraintsTree.getRoot(), this.subtrees);

        nodeToSubtreeMap = new int[getNodeCount()];
        Arrays.fill(nodeToSubtreeMap,-1);

        int visitedNodes = 0;
        for (int i=0; i<subtrees.size();i++) {
            WrappedSubtree subtree = subtrees.get(i);
            for (int j = 0; j < subtree.getNodeCount(); j++) {

                NodeRef node = subtree.getNode(j);
                NodeRef nodeInThisTree = convertSubtreeNodeToConstrainedNode(subtree,node);
                // subtree roots should be external nodes in parent clades. They are stored in the tree where they are the root
                if(subtree.isExternal(node)&&!isExternal(nodeInThisTree)){
                    continue;
                }
                if(nodeToSubtreeMap[nodeInThisTree.getNumber()]!=-1){
                    throw new IllegalArgumentException("The subtrees in a constrained tree must be monophyletic");
                }
                nodeToSubtreeMap[nodeInThisTree.getNumber()]=i;
                visitedNodes++;
            }
        }
        if(visitedNodes!=getNodeCount()){
            throw new IllegalArgumentException("The subtrees in a constrained tree must include all nodes in the tree");
        }

    }

    private WrappedSubtree getSubtrees(Tree contraintsTree,Map<NodeRef,NodeRef> treeMap,NodeRef node,List<WrappedSubtree>subtrees) {
        Set<NodeRef> tips = new HashSet<>();
        List<WrappedSubtree> childSubtrees = new ArrayList<>();
        for (int i = 0; i < contraintsTree.getChildCount(node); i++) {
            NodeRef child = contraintsTree.getChild(node, i);
            tips.add(treeMap.get(child));
            if(!contraintsTree.isExternal(child)){
                childSubtrees.add(getSubtrees(contraintsTree, treeMap, child,subtrees));
            }
        }
        WrappedSubtree subtree = new WrappedSubtree(wrappedTree, tips);
        subtrees.add(subtree);
        for (WrappedSubtree subtendedSubtree:
             childSubtrees) {
            subtendedSubtree.setParentTree(subtree);
        }
        return subtree;
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



    private NodeRef convertSubtreeNodeToConstrainedNode(WrappedSubtree subtree, NodeRef subtreeNode) {
        NodeRef wrappedNode = subtree.getNodeInWrappedTree(subtreeNode);
        return getUnWrappedNode(wrappedNode);
    }

    public WrappedSubtree getSubtree(NodeRef node){
        return subtrees.get(nodeToSubtreeMap[node.getNumber()]);
    }
    public WrappedSubtree getSubtree(int i){
        return subtrees.get(i);
    }

    public int getSubtreeIndex(NodeRef node){
        return nodeToSubtreeMap[node.getNumber()];
    }


    /**
     * Add child to the children of parent.
     *
     * @param parent
     * @param child
     * @throws IllegalArgumentException If child is already a child of parent
     */
    @Override
    public void addChild(NodeRef parent, NodeRef child) {
        throw new UnsupportedOperationException("Can not edit the topology of a constrained tree");
    }

    /**
     * Removes child from the children of parent.
     *
     * @param parent
     * @param child
     * @throws IllegalArgumentException If child is not a child of parent
     */
    @Override
    public void removeChild(NodeRef parent, NodeRef child) {
        throw new UnsupportedOperationException("Can not edit the topology of a constrained tree");

    }

    /**
     * Replace child with another
     *
     * @param node
     * @param child    of node to replace
     * @param newChild replacment child
     */
    @Override
    public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {
        throw new UnsupportedOperationException("Can not edit the topology of a constrained tree");
    }

    /**
     * set the height of the ith node in the tree (where the first n are internal).
     *
     * @param node
     * @param height
     */
    @Override
    public void setNodeHeight(NodeRef node, double height) {
        NodeRef wrappedNode = getNodeInWrappedTree(node);
        wrappedTree.setNodeHeight(wrappedNode,height);
    }

    @Override
    public void setNodeHeightQuietly(NodeRef node, double height) {
        NodeRef wrappedNode = getNodeInWrappedTree(node);
        wrappedTree.setNodeHeightQuietly(wrappedNode,height);
    }

    /**
     * set the rate of the ith node in the tree (where the first n are internal).
     *
     * @param node
     * @param rate
     */
    @Override
    public void setNodeRate(NodeRef node, double rate) {
        NodeRef wrappedNode = getNodeInWrappedTree(node);
        wrappedTree.setNodeRate(wrappedNode,rate);
    }

    /**
     * Sets an named attribute for a given node.
     *
     * @param node  the node whose attribute is being set.
     * @param name  the name of the attribute.
     * @param value the new value of the attribute.
     */
    @Override
    public void setNodeAttribute(NodeRef node, String name, Object value) {
        NodeRef wrappedNode = getNodeInWrappedTree(node);
        wrappedTree.setNodeAttribute(wrappedNode,name,value);
    }


    /**
     * Modifies the current tree by adopting the provided collection of edges
     *
     * @param edges       Edges are provided as index: child number; parent: array entry
     * @param nodeHeights Also sets the node heights to the provided values
     * @param childOrder  Array that contains whether a child node is left or right child
     * @param taxaNames
     */
    @Override
    public void adoptTreeStructure(int[] edges, double[] nodeHeights, int[] childOrder, String[] taxaNames) {
// Nothing to do I think???
    }

    @Override
    public boolean isTreeValid() {
        return wrappedTree.isTreeValid();
    }
    @Override
    public void setMultivariateTrait(NodeRef n, String name, double[] value) {
        wrappedTree.setMultivariateTrait(getNodeInWrappedTree(n), name, value);
    }


    /**
     * Will throw an exception if any nodes have this node as their children.
     *
     * @param root
     */
    @Override
    public void setRoot(NodeRef root) {
        throw new UnsupportedOperationException("Can not directly change the root of a constrained tree");
    }
    /**
     * @return root node of this tree.
     */
    @Override
    public NodeRef getRoot() {
        return getUnWrappedNode(wrappedTree.getRoot());
    }

    /**
     * @param node the node to test if root
     * @return whether the node is the root.
     */
    @Override
    public boolean isRoot(NodeRef node) {
        return node==getUnWrappedNode(wrappedTree.getRoot());
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void acceptState() {

    }

    @Override
    protected void storeStateHook() {

    }

    @Override
    protected void restoreStateHook() {

    }
    private int[] nodeToSubtreeMap;
    private List<WrappedSubtree> subtrees;
}
