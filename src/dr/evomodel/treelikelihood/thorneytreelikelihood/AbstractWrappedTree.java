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


import java.util.*;

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
//        storedNodeMap = new int [nodeCount];
        reverseNodeMap = new HashMap<>();
//        storedReverseNodeMap = new  HashMap<>();
        int[] wrappedNodes = new int[externalNodeCount];

        int i = 0;
        for (NodeRef wrappedNode : externalNodes) {
            Node node = new Node(i);
            nodes[i]= node;
            nodeMap[i]= wrappedNode.getNumber();
            reverseNodeMap.put(wrappedNode.getNumber(), node.getNumber());
            wrappedNodes[i] = wrappedNode.getNumber();
            i++;
        }

        NodeRef wrappedNode = TreeUtils.getCommonAncestor(tree, wrappedNodes);

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
            reverseNodeMap.put(node.getNumber(),  nodes[currentIndex].getNumber());
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
        Integer index = reverseNodeMap.get(node.getNumber());
        if(index==null){
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
            if(reverseNodeMap.containsKey(wrappedChild.getNumber())){
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
            if(reverseNodeMap.containsKey(wrappedChild.getNumber())){
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


    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState() {
//        System.arraycopy(nodeMap,0,storedNodeMap,0,nodeMap.length);
//        storedReverseNodeMap = new HashMap<>(reverseNodeMap);
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

//        int[] tmp1 = storedNodeMap;
//        storedNodeMap= nodeMap;
//        nodeMap=tmp1;
//
//        Map<Integer,Integer> tmp2 = storedReverseNodeMap;
//        storedReverseNodeMap= reverseNodeMap;
//        reverseNodeMap=tmp2;
        restoreStateHook();
    }

    abstract protected void restoreStateHook();


    protected void replaceNode(NodeRef currentWrappedNode,NodeRef replacement){
        NodeRef node = getUnWrappedNode(currentWrappedNode);
        nodeMap[node.getNumber()] = replacement.getNumber();
        reverseNodeMap.remove(currentWrappedNode.getNumber());

        reverseNodeMap.put(replacement.getNumber(),node.getNumber());
    }


    //Private to be moved down later
    private  int[] nodeMap;
//    private  int[] storedNodeMap;
    private Map<Integer, Integer> reverseNodeMap;
//    private Map<Integer, Integer> storedReverseNodeMap;


    protected final TreeModel wrappedTree;
    private final Node[] nodes;
    private final int nodeCount;
    private final int externalNodeCount;
    private final int internalNodeCount;
    private AbstractWrappedTree parentTree=null;

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
