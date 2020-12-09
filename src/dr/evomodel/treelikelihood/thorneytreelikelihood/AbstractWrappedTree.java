package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;


import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * This is a tree model that wraps another tree model and only provides access to a monophyletic subclade on that tree.
 * It has it's own Nodes, but relies on the base tree for the topology and data structure. It is up to each implementation
 * how to handle edits to the tree.
 */
public abstract class AbstractWrappedTree extends TreeModel {
    public static final String CONSTRAINED_TREE_MODEL = "constrainedTreeModel";

    public AbstractWrappedTree(String name, TreeModel tree, Set<NodeRef> externalNodes) {
        super(name, tree.isVariable());
        setId(name);

        this.wrappedTree = tree;
        addModel(tree);

        if(externalNodes==null){
            externalNodes = new HashSet<>();
            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                externalNodes.add(tree.getExternalNode(i));
            }
        }


        externalNodeCount = externalNodes.size();
        nodeCount = externalNodeCount*2-1;
        internalNodeCount = nodeCount - externalNodeCount;

        nodes = new Node[nodeCount];
        nodeMap = new int [nodeCount];
        storedNodeMap = new int [nodeCount];
        reverseNodeMap = new int [tree.getNodeCount()];
        storedReverseNodeMap = new int [tree.getNodeCount()];
        Arrays.fill(reverseNodeMap, -1);
        int[] wrappedNodes = new int[externalNodeCount];

        int i = 0;
        for (NodeRef wrappedNode : externalNodes) {
            Node node = new Node(i);
            nodes[i]= node;
            nodeMap[i]= wrappedNode.getNumber();
            reverseNodeMap[wrappedNode.getNumber()] = node.getNumber();
            wrappedNodes[i] = wrappedNode.getNumber();
            i++;
        }

        NodeRef wrappedNode = TreeUtils.getCommonAncestor(tree, wrappedNodes);
        NodeRef wrappedRoot = TreeUtils.getCommonAncestor(tree, wrappedNodes);

        traverseAndSetup(tree, wrappedNode, externalNodes, i);
//        int visitedTips =0;
//        do{
//            if(tree.isExternal(wrappedNode)){
//                if(!externalNodes.contains(wrappedNode)){
//                    throw new IllegalArgumentException("Wrapped paraphyletic wrapped trees are not supported");
//                }
//            }
//            if(!externalNodes.contains(wrappedNode)){
//                nodes[i] = new Node(i);
//                nodeMap[i]= wrappedNode.getNumber();
//                reverseNodeMap[wrappedNode.getNumber()] =  nodes[i].getNumber();
//                i++;
//            }else{
//                visitedTips++;
//                if(visitedTips==externalNodes.size()){
//                    done=true;
//                }
//            }
//            wrappedNode = TreeUtils.postorderSuccessor(tree, wrappedNode);
//
//        }while(!done);


    }
    public AbstractWrappedTree(String name, TreeModel tree) {
        this(name, tree, null);
    }

    private int traverseAndSetup(Tree tree, NodeRef node, Set<NodeRef> externalNodes,int currentIndex) {
        if(tree.isExternal(node)){
            if(!externalNodes.contains(node)){
                throw new IllegalArgumentException("Wrapped paraphyletic wrapped trees are not supported");
            }
        }
        if(!externalNodes.contains(node)){
            nodes[currentIndex] = new Node(currentIndex);
            nodeMap[currentIndex]= node.getNumber();
            reverseNodeMap[node.getNumber()] =  nodes[currentIndex].getNumber();
            currentIndex++;
            //new counter to count decendents;
            for (int i = 0; i < tree.getChildCount(node); i++) {
                currentIndex = traverseAndSetup(tree,tree.getChild(node,i),externalNodes,currentIndex);
            }
        }
        return currentIndex;
    }

    protected NodeRef getNodeInWrappedTree(NodeRef node) {
        return wrappedTree.getNode(nodeMap[node.getNumber()]);
    }

    protected NodeRef getUnWrappedNode(NodeRef node) {
        int index = reverseNodeMap[node.getNumber()];
        if(index==-1){
            return null;
        }
        return nodes[index] ;
    }


    /**
     * set the length of the branch above the ith node in the tree (where the first n are internal).
     *
     * @param node
     * @param length
     */
    @Override
    public void setBranchLength(NodeRef node, double length) {
        throw new UnsupportedOperationException("WrappedTreeModels cannot have branch lengths set");
    }



    @Override
    public double[] getMultivariateNodeTrait(NodeRef node, String name) {
        NodeRef wrappedNode = getNodeInWrappedTree(node);
        return wrappedTree.getMultivariateNodeTrait(wrappedNode, name);
    }


    /**
     * @return a count of the number of nodes (internal + external) in this
     * tree, currently connected from the root node.
     */
    @Override
    public int getNodeCount() {
        return nodeCount;
    }

    /**
     * @param i node index, terminal nodes are first
     * @return the ith node.
     */
    @Override
    public NodeRef getNode(int i) {
        return nodes[i];
    }

    /**
     * @param i index of an internal node
     * @return the ith internal node.
     */
    @Override
    public NodeRef getInternalNode(int i) {
        return nodes[externalNodeCount + i];
    }

    /**
     * @param i the index of an external node
     * @return the ith external node.
     */
    @Override
    public NodeRef getExternalNode(int i) {
        return nodes[i];
    }

    /**
     * @return a count of the number of external nodes (tips) in this
     * tree, currently connected from the root node.
     */
    @Override
    public int getExternalNodeCount() {
        return externalNodeCount;
    }

    /**
     * @return a count of the number of internal nodes in this
     * tree, currently connected from the root node.
     */
    @Override
    public int getInternalNodeCount() {
        return internalNodeCount;
    }

    /**
     * @param node the node to retrieve the taxon of
     * @return the taxon of this node.
     */
    @Override
    public Taxon getNodeTaxon(NodeRef node) {
        NodeRef wrappedNode = getNodeInWrappedTree((Node)node);
        return wrappedTree.getNodeTaxon(wrappedNode);
    }

    /**
     * @param node the node to retrieve height of
     * @return the height of node in the tree.
     */
    @Override
    public double getNodeHeight(NodeRef node) {
        NodeRef wrappedNode = getNodeInWrappedTree((Node)node);
        return wrappedTree.getNodeHeight(wrappedNode);
    }

    /**
     * @param node the node to retrieve the rate of
     * @return the rate of node in the tree.
     */
    @Override
    public double getNodeRate(NodeRef node) {
        NodeRef wrappedNode = getNodeInWrappedTree((Node)node);
        return wrappedTree.getNodeRate(wrappedNode);
    }

    /**
     * @param node the node whose attribute is being fetched.
     * @param name the name of the attribute of interest.
     * @return an object representing the named attributed for the given node.
     */
    @Override
    public Object getNodeAttribute(NodeRef node, String name) {
        NodeRef wrappedNode = getNodeInWrappedTree((Node)node);
        return wrappedTree.getNodeAttribute(wrappedNode,name);
    }

    /**
     * @param node the node whose attribute is being fetched.
     * @return an interator of attribute names available for this node.
     */
    @Override
    public Iterator getNodeAttributeNames(NodeRef node) {
        NodeRef wrappedNode = getNodeInWrappedTree((Node)node);
        return wrappedTree.getNodeAttributeNames(wrappedNode);
    }


    /**
     * @param node the node to test if external
     * @return whether the node is external.
     */
    @Override
    public boolean isExternal(NodeRef node) {
        return node.getNumber()<externalNodeCount;
    }


    /**
     * @param node the node to get child count of
     * @return the number of children of node.
     */
    @Override
    public int getChildCount(NodeRef node) {
        if(isExternal(node)){
            return 0;
        }
        NodeRef wrappedNode = getNodeInWrappedTree(node);
        int count=0;
        for(int i=0;i<wrappedTree.getChildCount(wrappedNode);i++){
            NodeRef wrappedChild = wrappedTree.getChild(wrappedNode, i);
            if(reverseNodeMap[wrappedChild.getNumber()]>-1){
                count++;
            }
        }
        return count;
    }



    /**
     * @param node the node to get jth child of
     * @param j    the index of child to retrieve
     * @return the jth child of node
     */
    @Override
    public NodeRef getChild(NodeRef node, int j) {
        if(isExternal(node)){
            throw new IllegalArgumentException("Attempted to access a child of an external node in a wrapped tree");
        }

        NodeRef wrappedNode = getNodeInWrappedTree(node);
        int currentChild=0;
        NodeRef wrappedJChild=null;
        for(int i=0;i<wrappedTree.getChildCount(wrappedNode);i++){
            NodeRef wrappedChild = wrappedTree.getChild(wrappedNode, i);
            if(reverseNodeMap[wrappedChild.getNumber()]>-1){
                if(currentChild==j){
                    wrappedJChild=wrappedChild;
                    break;
                }
                currentChild++;
            }
        }
        return getUnWrappedNode(wrappedJChild);
    }
    @Override
    public NodeRef getParent(NodeRef node) {
        if(isRoot(node)){
            return null;
        }else{
            NodeRef wrappedNode = getNodeInWrappedTree(node);
            NodeRef wrappedParent = wrappedTree.getParent(wrappedNode);

            return getUnWrappedNode(wrappedParent);
        }
    }



    public AbstractWrappedTree getParentTree(){
        return this.parentTree;
    }
    public void setParentTree(AbstractWrappedTree parentTree){
        this.parentTree = parentTree;
    }


    /**
     * @param taxonIndex
     * @return the ith taxon.
     */
    @Override
    public Taxon getTaxon(int taxonIndex) {
        int wrappedTaxonIndex = nodeMap[taxonIndex];
        return wrappedTree.getTaxon(wrappedTaxonIndex);
    }

    @Override
    public NodeRef[] getNodes() {
        return nodes;
    }

    @Override
    public boolean beginTreeEdit() {
        wrappedTree.beginTreeEdit();
        return super.beginTreeEdit();
    }

    @Override
    public void endTreeEdit() {
        wrappedTree.endTreeEdit();
        super.endTreeEdit();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if(model==wrappedTree){
            TreeChangedEvent treeChangedEvent = (TreeChangedEvent) object;
            NodeRef wrappedNode = treeChangedEvent.getNode();
            NodeRef node = wrappedNode == null ? null : getUnWrappedNode(wrappedNode);
            if(node!=null) {
                pushTreeChangedEvent(new RemappedTreeChangeEvent(treeChangedEvent, node));
            }
        }else{
            throw new RuntimeException("unknown model");
        }
    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     *
     * @param variable
     * @param index
     * @param type
     */
    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState() {
        System.arraycopy(nodeMap,0,storedNodeMap,0,nodeMap.length);
        System.arraycopy(reverseNodeMap,0,storedReverseNodeMap,0,reverseNodeMap.length);
        storeStateHook();
    }
    abstract protected void storeStateHook();

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState() {

        int[] tmp1 = storedNodeMap;
        storedNodeMap= nodeMap;
        nodeMap=tmp1;

        int[] tmp2 = storedReverseNodeMap;
        storedReverseNodeMap= reverseNodeMap;
        reverseNodeMap=tmp2;
        restoreStateHook();
    }

    abstract protected void restoreStateHook();


    protected void replaceNode(NodeRef currentWrappedNode,NodeRef replacement){
        NodeRef node = getUnWrappedNode(currentWrappedNode);
        nodeMap[node.getNumber()] = replacement.getNumber();
        reverseNodeMap[currentWrappedNode.getNumber()]=-1;
        reverseNodeMap[replacement.getNumber()] = node.getNumber();
    }


    //Private to be moved down later
    private  int[] nodeMap;
    private  int[] storedNodeMap;
    private  int[] reverseNodeMap;
    private  int[] storedReverseNodeMap;


    protected final TreeModel wrappedTree;
    private final Node[] nodes;
    private final int nodeCount;
    private final int externalNodeCount;
    private final int internalNodeCount;
    private AbstractWrappedTree parentTree=null;

    private class RemappedTreeChangeEvent implements TreeChangedEvent {

        final private TreeChangedEvent event;
        final private NodeRef node;

        private RemappedTreeChangeEvent(TreeChangedEvent event, NodeRef node) {
            this.event = event;
            this.node = node;
        }

        @Override public int getIndex() { return event.getIndex(); }

        @Override public NodeRef getNode() { return node; }

        @Override public Parameter getParameter() { return event.getParameter(); }

        @Override public boolean isNodeChanged() { return event.isNodeChanged(); }

        @Override public boolean isTreeChanged() { return event.isTreeChanged(); }

        @Override
        public boolean isNodeOrderChanged() { return event.isNodeOrderChanged(); }

        @Override public boolean isNodeParameterChanged() { return event.isNodeParameterChanged(); }

        @Override public boolean isHeightChanged() { return event.isHeightChanged(); }
    }
    private class Node implements NodeRef {

        private Node(int number) {
            this.number = number;
        }

        @Override
        public int getNumber() {
            return number;
        }

        @Override
        public void setNumber(int n) {
            throw new UnsupportedOperationException("Node is immutable");
        }
        final int number;

    }


}
