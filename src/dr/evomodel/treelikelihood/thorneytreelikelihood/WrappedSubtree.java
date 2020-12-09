package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.TreeModel;

import java.util.Set;

public class WrappedSubtree extends AbstractWrappedTree{

    public static final String WRAPPED_TREE_MODEL = "wrappedTreeModel";
    public WrappedSubtree(String name, TreeModel tree, Set<NodeRef> externalNodes){
        super(name,tree,externalNodes);

        int[] wrappedNodes = new int[externalNodes.size()];
        int i=0;
        for (NodeRef wrappedNode : externalNodes) {
            wrappedNodes[i] = wrappedNode.getNumber();
            i++;
        }

        NodeRef wrappedNode = TreeUtils.getCommonAncestor(tree, wrappedNodes);

        root = getUnWrappedNode(wrappedNode).getNumber();
        equivalentRoots = tree.isRoot(wrappedNode);
    }
    public WrappedSubtree(TreeModel tree, Set<NodeRef> externalNodes){
        this(WRAPPED_TREE_MODEL, tree, externalNodes);
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
        NodeRef wrappedParent = getNodeInWrappedTree(parent);
        NodeRef wrappedChild = getNodeInWrappedTree(child);

        if(isRoot(child)){
            NodeRef rootParent = wrappedTree.getParent(wrappedChild);
            wrappedTree.removeChild(rootParent, wrappedChild);
            subtendingNode = rootParent.getNumber();
        }
        wrappedTree.addChild(wrappedParent,wrappedChild);
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
        NodeRef wrappedParent = getNodeInWrappedTree(parent);
        NodeRef wrappedChild = getNodeInWrappedTree(child);
        wrappedTree.removeChild(wrappedParent,wrappedChild);
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
        NodeRef wrappedParent = getNodeInWrappedTree(node);
        NodeRef wrappedChild = getNodeInWrappedTree(child);
        NodeRef wrappedNewChild = getNodeInWrappedTree(newChild);
        wrappedTree.replaceChild(wrappedParent,wrappedChild,wrappedNewChild);
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
        wrappedTree.setNodeHeight(wrappedNode, height);
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
        wrappedTree.setNodeRate(wrappedNode, rate);
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
        wrappedTree.setNodeAttribute(wrappedNode, name,value);
    }

    @Override
    public void setMultivariateTrait(NodeRef node, String name, double[] value) {
        NodeRef wrappedNode = getNodeInWrappedTree(node);
        wrappedTree.setMultivariateTrait(wrappedNode, name,value);
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
//Nothing to do I think?
    }

    @Override
    public boolean isTreeValid() {
        return wrappedTree.isTreeValid()&& subtendingNode==-1;
    }

    public void setRoot(NodeRef node){
        int oldRoot = this.root;
        root=node.getNumber();
        if(equivalentRoots){
            wrappedTree.setRoot(getNodeInWrappedTree(getNode(root)));
        }else{
            NodeRef oldRootNode = getNodeInWrappedTree(getNode(oldRoot));
            NodeRef newRoodNode = getNodeInWrappedTree(getNode(root));
            parentTree.replaceNode(oldRootNode,newRoodNode);

            wrappedTree.addChild(wrappedTree.getNode(subtendingNode),getNodeInWrappedTree(getNode(root)));
            subtendingNode=-1;
        }

//        throw new UnsupportedOperationException("wrapped trees can not change roots, yet");
    }
    public NodeRef getRoot(){
        return getNode(root);
    }
    public boolean isRoot(NodeRef nodeRef) {
        return nodeRef.getNumber()==root;
    }


    public void setParentTree(WrappedSubtree tree) {
        parentTree=tree;
    }
    @Override
    public void setNodeHeightQuietly(NodeRef node, double height) {
        NodeRef wrappedNode = getNodeInWrappedTree(node);
        wrappedTree.setNodeHeightQuietly(wrappedNode, height);
    }


    protected void storeStateHook(){
        storedRoot = root;
    }
    protected  void restoreStateHook(){
        root = storedRoot;
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void acceptState() {

    }

    private int subtendingNode = -1;
    private final boolean equivalentRoots;
    private int root;
    private int storedRoot = -1;
    private WrappedSubtree parentTree;
}
