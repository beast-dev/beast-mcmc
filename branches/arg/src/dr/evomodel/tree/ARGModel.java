/*
 * ARGModel.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */

package dr.evomodel.tree;

import dr.evolution.tree.*;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.evomodel.treelikelihood.ARGLikelihood;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.math.MathUtils;
import dr.util.Attributable;
import dr.xml.*;
import org.jdom.Document;
import org.jdom.Element;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

/**
 * A model component for trees.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: ARGModel.java,v 1.18.2.4 2006/11/06 01:38:30 msuchard Exp $
 */
public class ARGModel extends AbstractModel
        implements MutableTree, Loggable {

    //
    // Public stuff
    //

    public static final String TREE_MODEL = "argTreeModel";

    public static final String ROOT_HEIGHT = "rootHeight";
    public static final String LEAF_HEIGHT = "leafHeight";

    public static final String NODE_HEIGHTS = "nodeHeights";
    public static final String NODE_RATES = "nodeRates";
    public static final String NODE_TRAITS = "nodeTraits";

    public static final String ROOT_NODE = "rootNode";
    public static final String INTERNAL_NODES = "internalNodes";
    public static final String LEAF_NODES = "leafNodes";
    public static final String TAXON = "taxon";

    public ARGModel(Tree tree) {

        super(TREE_MODEL);
//        System.err.println("constructor for TreeModel");
        partitioningParameters = new VariableSizeCompoundParameter("partitioning");
        //     initialize(tree);
        //}

        //protected void initialize(Tree tree) {
        //System.err.println("init for TreeModel");
        //System.exit(-1);
        // get a rooted version of the tree to clone
        FlexibleTree binaryTree = new FlexibleTree(tree);
        binaryTree.resolveTree();

        // clone the node structure (this will create the individual parameters
        Node node = new Node(binaryTree, binaryTree.getRoot());

        internalNodeCount = binaryTree.getInternalNodeCount();
        externalNodeCount = binaryTree.getExternalNodeCount();

        nodeCount = internalNodeCount + externalNodeCount;

        //nodes = new Node[nodeCount];
        //storedNodes = new Node[nodeCount];
        nodes = new ArrayList<Node>(nodeCount);
        storedNodes = new ArrayList<Node>(nodeCount);
        for (int i = 0; i < nodeCount; i++) {
            nodes.add(null);
            storedNodes.add(null);
        }
        int i = 0;
        int j = externalNodeCount;

        root = node;

        //System.err.println("Going to do postOrder");

        do {
            node = (Node) Tree.Utils.postorderSuccessor(this, node);

            if (node.isExternal()) {
                node.number = i;

                //nodes[i] = node;
                //storedNodes[i] = new Node();
                //storedNodes[i].taxon = node.taxon;
                //storedNodes[i].number = i;
                nodes.set(i, node);
                Node copy = new Node();
                copy.taxon = node.taxon;
                copy.number = i;
                storedNodes.set(i, copy);

                i++;
            } else {
                node.number = j;

                //nodes[j] = node;
                //storedNodes[j] = new Node();
                //storedNodes[j].number = j;
                nodes.set(j, node);
                Node copy = new Node();
                copy.number = j;
                storedNodes.set(j, copy);

                j++;
            }
        } while (node != root);

        //       System.err.println("Succeed in post-order");

        //       ARGTree t = new ARGTree(this,0);
        //       System.err.println(Tree.Utils.uniqueNewick(t, t.getRoot()));
        //       System.err.println(this.toGraphString());
        // System.exit(-1);

    }

    public VariableSizeCompoundParameter getPartitioningParameters() {
        return partitioningParameters;
    }

    protected void setupHeightBounds() {
        for (Node node : nodes) {
            node.setupHeightBounds();
        }
    }


    public static final String GRAPH_ELEMENT = "graph";
    public static final String NODE_ELEMENT = "node";
    public static final String EDGE_ELEMENT = "edge";
    public static final String ID_ATTRIBUTE = "id";
    public static final String EDGE_TO = "source";
    public static final String EDGE_FROM = "target";
    public static final String TAXON_NAME = "taxonName";
    public static final String EDGE_LENGTH = "len";
    public static final String IS_TIP = "isTip";

    private String getNameOfNode(Node node) {
        if (node.taxon == null)
            return "n" + Integer.toString(node.number);
        else
            return node.taxon.getId();
    }

    private Element makeEdge(Node from, Node to) {
        Element edgeElement = new Element(EDGE_ELEMENT);
        edgeElement.setAttribute(EDGE_FROM, getNameOfNode(from));
        edgeElement.setAttribute(EDGE_TO, getNameOfNode(to));
        edgeElement.setAttribute(EDGE_LENGTH, Double.toString(getNodeHeight(from) - getNodeHeight(to)));
        return edgeElement;
    }

    private Element makeNode(Node node) {
        Element nodeElement = new Element(NODE_ELEMENT);
        nodeElement.setAttribute(ID_ATTRIBUTE, getNameOfNode(node));
        if (node.taxon != null) {
            nodeElement.setAttribute(IS_TIP, "true");
//			nodeElement.setAttribute(TAXON_NAME, node.taxon.getId());
//			nodeElement.setAttribute("style","filled");
//			nodeElement.setAttribute("fillcolor","blue");
        }
        return nodeElement;
    }

    public Element toXML() {
        int cnt = 0;
        for (Node node : nodes)
            node.number = cnt++;

        Element graphElement = new Element(GRAPH_ELEMENT);
        graphElement.setAttribute("edgedefault", "directed");

        for (Node node : nodes) {

            graphElement.addContent(makeNode(node));
            // Add edge to left parent if not root
            if (node.leftParent != null) {
                graphElement.addContent(makeEdge(node.leftParent, node));

            }
            // Add edge to right parent if reassortment
            if (node.rightParent != null && node.isReassortment()) {
                graphElement.addContent(makeEdge(node.rightParent, node));
            }

        }

        return graphElement;
    }


    /**
     * Push a tree changed event into the event stack.
     */
    public void pushTreeChangedEvent() {
        pushTreeChangedEvent(new TreeChangedEvent());
    }


    public void pushTreeSizeChangedEvent() {
        pushTreeChangedEvent(new TreeChangedEvent(this));
    }

    /**
     * Push a tree changed event into the event stack.
     */
    public void pushTreeChangedEvent(NodeRef nodeRef) {
        pushTreeChangedEvent(new TreeChangedEvent((Node) nodeRef));
    }

    /**
     * Push a tree changed event into the event stack.
     */
    public void pushTreeChangedEvent(Node node, Parameter parameter, int index) {
        pushTreeChangedEvent(new TreeChangedEvent(node, parameter, index));
    }

    /**
     * Push a tree changed event into the event stack.
     */
    public void pushTreeChangedEvent(TreeChangedEvent event) {
        if (inEdit) {
            treeChangedEvents.add(event);
        } else {
            listenerHelper.fireModelChanged(this, event);
        }
    }


    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // no submodels so nothing to do
    }

    /**
     * Called when a parameter changes.
     */
    public void handleParameterChangedEvent(Parameter parameter, int index) {

        Node node = getNodeOfParameter(parameter);
        pushTreeChangedEvent(node, parameter, index);
    }

    private ArrayList<ARGLikelihood> likelihoodCalculators;
    private int maxNumberOfPartitions;


    public int getNumberOfPartitions() {
        return maxNumberOfPartitions;
    }


    public int addLikelihoodCalculator(ARGLikelihood calc) {
        //int len = 0;
        if (likelihoodCalculators == null) {
            likelihoodCalculators = new ArrayList<ARGLikelihood>();
        }
        likelihoodCalculators.add(calc);
        int len = likelihoodCalculators.size() - 1;
        maxNumberOfPartitions = likelihoodCalculators.size();
        System.err.println("Add calculator for partition #" + len);
        setPartitionRecursively(getRoot(), len);
        return len;
    }

    public int getMaxPartitionNumber() {
        return maxNumberOfPartitions;
    }

    protected List treeChangedEvents = new ArrayList();

    public class TreeChangedEvent {

        Node node;
        Parameter parameter;
        int index;
        boolean size = false;

        public TreeChangedEvent() {
            this(null, null, -1);
        }

        public TreeChangedEvent(ARGModel arg) {
            this(null, null, -1);
            size = true;
        }

        public TreeChangedEvent(Node node) {
            this(node, null, -1);
        }

        public TreeChangedEvent(Node node, Parameter parameter, int index) {
            this.node = node;
            this.parameter = parameter;
            this.index = index;
        }

        public int getIndex() {
            return index;
        }

        public Node getNode() {
            return node;
        }

        public Parameter getParameter() {
            return parameter;
        }

        public boolean isSizeChanged() {
            return size;
        }

        public boolean isNodeChanged() {
            return node != null;
        }

        public boolean isNodeParameterChanged() {
            return parameter != null;
        }

        public boolean isHeightChanged() {
            return parameter == node.heightParameter;
        }

        public boolean isRateChanged() {
            return parameter == node.rateParameter;
        }

        public boolean isTraitChanged() {
            return parameter == node.traitParameter;
        }
    }

    // *****
    // Interface Loggable
    // *****

    public LogColumn[] getColumns() {
        int numColumns = 3;
        //numColumns += this.getMaxPartitionNumber();
        LogColumn[] logColumns = new LogColumn[numColumns + getMaxPartitionNumber()];
        logColumns[0] = new IsReassortmentColumn("isReassortment");
        logColumns[1] = new CountReassortmentColumn("numberReassortments");
        logColumns[2] = new IsRootTooHighColumn("isRootTooHigh");
        for (int i = 0; i < getMaxPartitionNumber(); i++)
            logColumns[numColumns + i] = new ArgTreeHeightColumn("argTreeHeight", this, i);
        return logColumns;
    }

    private class ArgTreeHeightColumn extends NumberColumn {

        private int partition;
        private ARGModel argModel;

        public ArgTreeHeightColumn(String label, ARGModel argModel, int partition) {
            super(label + partition);
            this.argModel = argModel;
            this.partition = partition;
        }

        public double getDoubleValue() {
            ARGTree argTree = new ARGTree(argModel, partition);
            return argTree.getNodeHeight(argTree.getRoot());
//			return  (new ARGTree(
        }
    }


    private class IsReassortmentColumn extends NumberColumn {

        public IsReassortmentColumn(String label) {
            super(label);    //To change body of overridden methods use File | Settings | File Templates.
        }

        public double getDoubleValue() {
            return getReassortmentNodeCount() == 0 ? 0 : 1;
        }
    }

    private class IsRootTooHighColumn extends NumberColumn {

        public IsRootTooHighColumn(String label) {
            super(label);    //To change body of overridden methods use File | Settings | File Templates.
        }

        public double getDoubleValue() {
            return isBifurcationDoublyLinked(getRoot()) ? 1 : 0;
        }
    }

    private class CountReassortmentColumn extends NumberColumn {

        public CountReassortmentColumn(String label) {
            super(label);    //To change body of overridden methods use File | Settings | File Templates.
        }

        public double getDoubleValue() {
            return getReassortmentNodeCount();
        }
    }

    // *****************************************************************
    // Interface Tree
    // *****************************************************************

    /**
     * Return the units that this tree is expressed in.
     */
    public final int getUnits() {
        return units;
    }

    /**
     * Sets the units that this tree is expressed in.
     */
    public final void setUnits(int units) {
        this.units = units;
    }

    /**
     * @return a count of the number of nodes (internal + external) in this
     *         tree.
     */
    public final int getNodeCount() {
        return nodes.size();
    }

    public final boolean hasNodeHeights() {
        return true;
    }

    public final double getNodeHeight(NodeRef node) {

        //System.err.println(Tree.Utils.uniqueNewick(this, node));
        //((Node)node))

        return ((Node) node).getHeight();
    }

    public final double getMinParentNodeHeight(NodeRef nr) {
        Node node = (Node) nr;
        return Math.min(node.leftParent.getHeight(), node.rightParent.getHeight());
    }


    public final double getNodeHeightUpper(NodeRef node) {
        return ((Node) node).heightParameter.getBounds().getUpperLimit(0);
    }

    public final double getNodeHeightLower(NodeRef node) {
        return ((Node) node).heightParameter.getBounds().getLowerLimit(0);
    }


    /**
     * @param node
     * @return the rate parameter associated with this node.
     */
    public final double getNodeRate(NodeRef node) {
        if (!hasRates) {
            return 1.0;
        }
        return ((Node) node).getRate();
    }

    public Object getNodeAttribute(NodeRef node, String name) {
        throw new UnsupportedOperationException("ARGModel does not use NodeAttributes");
    }

    public Iterator getNodeAttributeNames(NodeRef node) {
        throw new UnsupportedOperationException("ARGModel does not use NodeAttributes");
    }

    public double getNodeTrait(NodeRef node) {
        if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
        return ((Node) node).getTrait();
    }

    public final Taxon getNodeTaxon(NodeRef node) {
        return ((Node) node).taxon;
    }

    public final boolean isExternal(NodeRef node) {
        return ((Node) node).isExternal();
    }

    public final boolean isRoot(NodeRef node) {
        return (node == root);
    }

    public final boolean isBifurcation(NodeRef node) {
        return ((Node) node).isBifurcation();
    }

    public final boolean isBifurcationDoublyLinked(NodeRef node) {
        return ((Node) node).isBifurcationDoublyLinked();
    }

    public final boolean isReassortment(NodeRef node) {
        return ((Node) node).isReassortment();
    }

    public final int countReassortmentNodes(NodeRef nr) {
        Node node = (Node) nr;
        int count = node.countReassortmentChild(this);
        //int count = 0;
        return (count / 2);
    }

    public final int getChildCount(NodeRef node) {
        return ((Node) node).getChildCount();
    }

    public final NodeRef getChild(NodeRef node, int i) {
        return ((Node) node).getChild(i);
    }

    public final NodeRef getChild(NodeRef node, int i, int partition) {
        return ((Node) node).getChild(i, partition);
    }

    //public final NodeRef getParent(NodeRef node) { return ((Node)node).parent; }
    public final NodeRef getParent(NodeRef node) {
        Node left = ((Node) node).leftParent;
        Node right = ((Node) node).rightParent;
        if (left == right)
            return left;
        else
            throw new IllegalArgumentException("No single parent for reassorted node");
    }

    public final NodeRef getParent(NodeRef node, int i) {
        if (i == 0)
            return ((Node) node).leftParent;
        if (i == 1)
            return ((Node) node).rightParent;
        throw new IllegalArgumentException("ARGModel.Node can only have two parents");
    }

//    public final NodeRef getParent(NodeRef node, int partition) {
//    	// TODO
//    	return null;
//    }
//    

    public final boolean hasBranchLengths() {
        return true;
    }

    public final double getBranchLength(NodeRef node) {
        NodeRef parent = getParent(node);
        if (parent == null) {
            return 0.0;
        }

        return getNodeHeight(parent) - getNodeHeight(node);
    }

    public final NodeRef getExternalNode(int i) {
        return nodes.get(i);
    }

    public final NodeRef getInternalNode(int i) {
        return nodes.get(i + externalNodeCount);
    }

    public final NodeRef getNode(int i) {
        return nodes.get(i);
    }

    /**
     * Returns the number of external nodes.
     */
    public final int getExternalNodeCount() {
        return externalNodeCount;
    }

    /**
     * Returns the ith internal node.
     */
    public final int getInternalNodeCount() {
        return internalNodeCount;
    }

    public final int getReassortmentNodeCount() {
        int cnt = 0;
        for (Node node : nodes) {
            if (!node.bifurcation)
                cnt++;
        }
        return cnt;
    }

    /**
     * Returns the root node of this tree.
     */
    public final NodeRef getRoot() {
        return root;
    }

    public final NodeRef getRoot(int partition) {
        // TODO
        return null;
    }

    // *****************************************************************
    // Interface MutableTree
    // *****************************************************************

    /**
     * Set a new node as root node.
     */
    public final void setRoot(NodeRef newRoot) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        root = (Node) newRoot;

        // We shouldn't need this because the addChild will already have fired appropriate events.
        // pushTreeChangedEvent();
    }

    public void swapHeightParameters(NodeRef n1, NodeRef n2) {
        Node node1 = (Node) n1;
        Node node2 = (Node) n2;
        double height1 = node1.getHeight();
        double height2 = node2.getHeight();
        Parameter trans = node1.heightParameter;
        node1.heightParameter = node2.heightParameter;
        node2.heightParameter = trans;
        node1.setHeight(height1);
        node2.setHeight(height2);

    }

    public void addChild(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        Node parent = (Node) p;
        Node child = (Node) c;
        parent.doubleAddChild(child);
    }

    public void singleAddChild(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        Node parent = (Node) p;
        Node child = (Node) c;
        parent.singleAddChild(child);
    }

    public void singleAddChildWithOneParent(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        Node parent = (Node) p;
        Node child = (Node) c;
        parent.singleAddChildWithOneParent(child);
    }

    public void doubleAddChild(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        Node parent = (Node) p;
        Node child = (Node) c;
        parent.doubleAddChild(child);
    }

    public void doubleAddChildWithOneParent(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        Node parent = (Node) p;
        Node child = (Node) c;
        parent.doubleAddChildWithOneParent(child);
    }

    public void addChildAsRecombinant(NodeRef p1, NodeRef p2, NodeRef c,
                                      Parameter partitioning) {
//    public void addChildAsRecombinant(NodeRef p1, NodeRef p2, NodeRef c, BitSet bs1, BitSet bs2) {
        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");
        Node parent1 = (Node) p1;
        Node parent2 = (Node) p2;
        Node child = (Node) c;
        //if (parent1.hasChild(child) || parent2.hasChild(child)) throw new IllegalArgumentException("Child already exists in parent");
        //if (parent2.hasChild(child)) throw new IllegalArgumentException("Child already exists in")
        parent1.addChildRecombinant(child, partitioning);
        parent2.addChildRecombinant(child, partitioning);
        //if( parent2.getChildCount() == 1 )
        //	parent2.addChildNoParentConnection(node);
    }

    public void removeChild(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        Node parent = (Node) p;
        Node child = (Node) c;

        parent.doubleRemoveChild(child);
    }

    public void doubleRemoveChild(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        Node parent = (Node) p;
        Node child = (Node) c;

        parent.doubleRemoveChild(child);
    }

    public void singleRemoveChild(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        Node parent = (Node) p;
        Node child = (Node) c;

        parent.singleRemoveChild(child);
    }

    protected Node oldRoot;

    public void beginTreeEdit() {
        if (inEdit) throw new RuntimeException("Alreading in edit transaction mode!");

        oldRoot = root;

        inEdit = true;
    }

    public void endTreeEdit() throws MutableTree.InvalidTreeException {
        if (!inEdit) throw new RuntimeException("Not in edit transaction mode!");

        inEdit = false;

        if (root != oldRoot) {
            swapParameterObjects(oldRoot, root);
        }

        //for (int i =0; i < nodes.length; i++) {
        for (Node node : nodes) {
            if (!node.heightParameter.isWithinBounds()) {
                throw new MutableTree.InvalidTreeException("height parameter out of bounds");
            }
        }
//ystem.err.println("There are "+treeChangedEvents.size()+" events waiting");
//System.exit(-1);
        for (int i = 0; i < treeChangedEvents.size(); i++) {
            listenerHelper.fireModelChanged(this, treeChangedEvents.get(i));
        }
        treeChangedEvents.clear();
    }

    public void setNodeHeight(NodeRef n, double height) {
        ((Node) n).setHeight(height);
    }


    public void setNodeRate(NodeRef n, double rate) {
        if (!hasRates) throw new IllegalArgumentException("Rate parameters have not been created");
        ((Node) n).setRate(rate);

    }

    public void setNodeTrait(NodeRef n, double value) {
        if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
        ((Node) n).setTrait(value);
    }

    public void setBranchLength(NodeRef node, double length) {
        throw new UnsupportedOperationException("ARGModel cannot have branch lengths set");
    }

    public void setNodeAttribute(NodeRef node, String name, Object value) {
        throw new UnsupportedOperationException("ARGModel does not use NodeAttributes");
    }

    // *****************************************************************
    // Interface ModelComponent
    // *****************************************************************

    /**
     * Store current state
     */
    protected void storeState() {
/*		System.err.println("Storing state");
		this.checkBranchSanity();
		System.err.println("sane before operation");*/
        copyNodeStructure(storedNodes);
        //storedRootNumber = storedNodes.indexOf(root.getNumber();
        storedRootNumber = nodes.indexOf(root);
        storedNodeCount = nodeCount;
        storedInternalNodeCount = internalNodeCount;
        addedParameters = null;
        addedPartitioningParameter = null;
        addedNodes = null;
        removedParameters = null;
        removedPartitioningParameter = null;
        removedNodes = null;
        //System.err.println("Stored: "+Tree.Utils.uniqueNewick(this, getRoot()));
        //System.err.println("Stored : "+this.toString());
    }

    /**
     * Restore the stored state
     */
    protected void restoreState() {
//		System.err.println("Restoring state.");
        ArrayList<Node> tmp = storedNodes;
        storedNodes = nodes;
        nodes = tmp;
        root = nodes.get(storedRootNumber);

        nodeCount = storedNodeCount;
        internalNodeCount = storedInternalNodeCount;

        if (addedParameters != null) {
            storedInternalNodeHeights.removeParameter(addedParameters[0]);
            storedInternalNodeHeights.removeParameter(addedParameters[1]);
            removeParameter(addedParameters[0]);
            removeParameter(addedParameters[1]);

            storedInternalAndRootNodeHeights.removeParameter(addedParameters[0]);
            storedInternalAndRootNodeHeights.removeParameter(addedParameters[1]);
        }

        if (addedPartitioningParameter != null) {
            partitioningParameters.removeParameter(addedPartitioningParameter);
            removeParameter(addedPartitioningParameter);
        }
        if (removedParameters != null) {
            storedInternalNodeHeights.addParameter(removedParameters[0]);
            storedInternalNodeHeights.addParameter(removedParameters[1]);
            addParameter(removedParameters[0]);
            addParameter(removedParameters[1]);

            storedInternalAndRootNodeHeights.addParameter(removedParameters[0]);
            storedInternalAndRootNodeHeights.addParameter(removedParameters[1]);


        }

        if (removedPartitioningParameter != null) {
            partitioningParameters.addParameter(removedPartitioningParameter);
            addParameter(removedPartitioningParameter);
        }
        //System.err.println("Restore: "+Tree.Utils.uniqueNewick(this, getRoot()));
        //System.err.println("Restore: "+this.toString());
    }

    /**
     * accept the stored state
     */
    protected void acceptState() {
        //System.err.println("Accepted ARG\n"+this.toGraphString());
    } // nothing to do

    /**
     * Adopt the state of the model component from source.
     */
    protected void adoptState(Model source) {
    }
/*
	public void addNewHeightParameters(Parameter newbie1, Parameter newbie2,
	                                   VariableSizeCompoundParameter internalNodeParameters,
                                       VariableSizeCompoundParameter internalAndRootNodeParameters) {
		addParameter(newbie1);
		addParameter(newbie2);
		addedParameters = new Parameter[2];
		addedParameters[0] = newbie1;
		addedParameters[1] = newbie2;

        storedInternalNodeHeights = internalNodeParameters;
		storedInternalNodeHeights.addParameter(newbie1);
		storedInternalNodeHeights.addParameter(newbie2);

 *//*       storedInternalAndRootNodeHeights = internalAndRootNodeParameters;
		storedInternalAndRootNodeHeights.addParameter(newbie1);
		storedInternalAndRootNodeHeights.addParameter(newbie2);*//*
    }*/

    public void expandARGWithRecombinant(Node newbie1, Node newbie2,
                                         VariableSizeCompoundParameter internalNodeParameters,
                                         VariableSizeCompoundParameter internalAndRootNodeParameters) {
        addParameter(newbie1.heightParameter);
        addParameter(newbie2.heightParameter);
        addParameter(newbie2.partitioning);
        addedParameters = new Parameter[2];
        addedParameters[0] = newbie1.heightParameter;
        addedParameters[1] = newbie2.heightParameter;
        addedPartitioningParameter = newbie2.partitioning;

        storedInternalNodeHeights = internalNodeParameters;
        storedInternalNodeHeights.addParameter(newbie1.heightParameter);
        storedInternalNodeHeights.addParameter(newbie2.heightParameter);

        storedInternalAndRootNodeHeights = internalAndRootNodeParameters;
        storedInternalAndRootNodeHeights.addParameter(newbie1.heightParameter);
        storedInternalAndRootNodeHeights.addParameter(newbie2.heightParameter);

        partitioningParameters.addParameter(newbie2.partitioning);
        nodes.add(newbie1);
        nodes.add(newbie2);
        internalNodeCount += 2;
        sanityNodeCheck(internalNodeParameters);
    }

    public void sanityNodeCheck(VariableSizeCompoundParameter inodes) {
        int len = inodes.getNumParameters();
        for (int i = 0; i < len; i++) {
            Parameter p = inodes.getParameter(i);
            for (int j = 0; j < internalNodeCount; j++) {
                Node node = (Node) getInternalNode(j);
                if (node.heightParameter == p) {
                    if (isRoot(node)) {
                        System.err.println("Root height found in internal nodes");
                        System.exit(-1);
                    }
                }


            }
        }
    }


    public void contractARGWithRecombinant(Node oldie1, Node oldie2,
                                           VariableSizeCompoundParameter internalNodeParameters,
                                           VariableSizeCompoundParameter internalAndRootNodeParameters) {
        removeParameter(oldie1.heightParameter);
        removeParameter(oldie2.heightParameter);
        removeParameter(oldie2.partitioning);
        removedParameters = new Parameter[2];
        removedParameters[0] = oldie1.heightParameter;
        removedParameters[1] = oldie2.heightParameter;
        partitioningParameters.removeParameter(oldie2.partitioning);
        removedPartitioningParameter = oldie2.partitioning;
        storedInternalNodeHeights = internalNodeParameters;
        storedInternalNodeHeights.removeParameter(oldie1.heightParameter);
        storedInternalNodeHeights.removeParameter(oldie2.heightParameter);

        storedInternalAndRootNodeHeights = internalAndRootNodeParameters;
        storedInternalAndRootNodeHeights.removeParameter(oldie1.heightParameter);
        storedInternalAndRootNodeHeights.removeParameter(oldie2.heightParameter);

        nodes.remove(oldie1);
        nodes.remove(oldie2);
        internalNodeCount -= 2;
    }


    /**
     * Copies the node connections from this ARGModel's nodes array to the
     * destination array. Basically it connects up the nodes in destination
     * in the same way as this ARGModel is set up. This method is package
     * private.
     */
    void copyNodeStructure(ArrayList<Node> destination) {

        //if ( nodes.length != destination.length ) {
        //    throw new IllegalArgumentException("Node arrays are of different lengths");
        //}
        while (destination.size() < nodes.size())
            destination.add(new Node());
        while (destination.size() > nodes.size())
            destination.remove(0);
        int n = nodes.size();
        for (int i = 0; i < n; i++) {
            //for( Node node0 : nodes ) {
            Node node0 = nodes.get(i);
            Node node1 = destination.get(i);

            // the parameter values are automatically stored and restored
            // just need to keep the links
            node1.heightParameter = node0.heightParameter;
            node1.rateParameter = node0.rateParameter;
            node1.traitParameter = node0.traitParameter;
            node1.partitioning = node0.partitioning;

            node1.taxon = node0.taxon;
            node1.bifurcation = node0.bifurcation;
            node1.number = node0.number;
            //node1.partitionSet = (BitSet)node0.partitionSet.clone();
//			if (node0.leftPartition != null) {
//				node1.leftPartition = (BitSet) node0.leftPartition.clone();
//			} else {
//				node1.leftPartition = null;
//			}
//			if (node0.rightPartition != null) {
//				node1.rightPartition = (BitSet) node0.rightPartition.clone();
//			} else {
//				node1.rightPartition = null;
//			}
//            


            if (node0.leftParent != null) {
                node1.leftParent = //storedNodes.get(node0.leftParent.getNumber());
                        storedNodes.get(nodes.indexOf(node0.leftParent));
            } else {
                node1.leftParent = null;
            }

            if (node0.rightParent != null) {
                node1.rightParent = //storedNodes.get(node0.rightParent.getNumber());
                        storedNodes.get(nodes.indexOf(node0.rightParent));
            } else {
                node1.rightParent = null;
            }

            if (node0.leftChild != null) {
                node1.leftChild = //storedNodes.get(node0.leftChild.getNumber());
                        storedNodes.get(nodes.indexOf(node0.leftChild));
            } else {
                node1.leftChild = null;
            }

            if (node0.rightChild != null) {
                node1.rightChild = //storedNodes.get(node0.rightChild.getNumber());
                        storedNodes.get(nodes.indexOf(node0.rightChild));
            } else {
                node1.rightChild = null;
            }
        }
    }

    public void setPartitionRecursively(NodeRef nr, int partition) {
        Node node = (Node) nr;
        node.setPartitionRecursively(partition);
    }


    /**
     * @return the number of statistics of this component.
     */
    public int getStatisticCount() {
        return 1;
    }

    /**
     * @return the ith statistic of the component
     */
    public Statistic getStatistic(int i) {
        if (i == 0) return root.heightParameter;
        throw new IllegalArgumentException();
    }

    public String getModelComponentName() {
        return TREE_MODEL;
    }

    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

    /**
     * @return a count of the number of taxa in the list.
     */
    public int getTaxonCount() {
        return getExternalNodeCount();
    }

    /**
     * @return the ith taxon in the list.
     */
    public Taxon getTaxon(int taxonIndex) {
        return ((Node) getExternalNode(taxonIndex)).taxon;
    }

    /**
     * @return the ID of the taxon of the ith external node. If it doesn't have
     *         a taxon, returns the ID of the node itself.
     */
    public String getTaxonId(int taxonIndex) {
        Taxon taxon = getTaxon(taxonIndex);
        if (taxon != null) {
            return taxon.getId();
        } else {
            return null;
        }
    }

    /**
     * returns the index of the taxon with the given id.
     */
    public int getTaxonIndex(String id) {
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            if (getTaxonId(i).equals(id)) return i;
        }
        return -1;
    }

    /**
     * returns the index of the given taxon.
     */
    public int getTaxonIndex(Taxon taxon) {
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            if (getTaxon(i) == taxon) return i;
        }
        return -1;
    }

    /**
     * @param taxonIndex the index of the taxon whose attribute is being fetched.
     * @param name       the name of the attribute of interest.
     * @return an object representing the named attributed for the taxon of the given
     *         external node. If the node doesn't have a taxon then the nodes own attribute
     *         is returned.
     */
    public final Object getTaxonAttribute(int taxonIndex, String name) {
        Taxon taxon = getTaxon(taxonIndex);
        if (taxon != null) {
            return taxon.getAttribute(name);
        }
        return null;
    }

    // **************************************************************
    // MutableTaxonList IMPLEMENTATION
    // **************************************************************

    public int addTaxon(Taxon taxon) {
        throw new IllegalArgumentException("Cannot add taxon to a ARGModel");
    }

    public boolean removeTaxon(Taxon taxon) {
        throw new IllegalArgumentException("Cannot add taxon to a ARGModel");
    }

    public void setTaxonId(int taxonIndex, String id) {
        throw new IllegalArgumentException("Cannot set taxon id in a ARGModel");
    }

    public void setTaxonAttribute(int taxonIndex, String name, Object value) {
        throw new IllegalArgumentException("Cannot set taxon attribute in a ARGModel");
    }

    public void addMutableTreeListener(MutableTreeListener listener) {
    } // Do nothing at the moment

    public void addMutableTaxonListListener(MutableTaxonListListener listener) {
    } // Do nothing at the moment

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

    protected String id = null;

    /**
     * @return the id.
     */
    public String getId() {
        return id;
    }

    /**
     * Sets the id.
     */
    public void setId(String id) {
        this.id = id;
    }

    // **************************************************************
    // Attributable IMPLEMENTATION
    // **************************************************************

    private Attributable.AttributeHelper treeAttributes = null;

    /**
     * Sets an named attribute for this object.
     *
     * @param name  the name of the attribute.
     * @param value the new value of the attribute.
     */
    public void setAttribute(String name, Object value) {
        if (treeAttributes == null)
            treeAttributes = new Attributable.AttributeHelper();
        treeAttributes.setAttribute(name, value);
    }

    /**
     * @param name the name of the attribute of interest.
     * @return an object representing the named attributed for this object.
     */
    public Object getAttribute(String name) {
        if (treeAttributes == null)
            return null;
        else
            return treeAttributes.getAttribute(name);
    }

    /**
     * @return an iterator of the attributes that this object has.
     */
    public Iterator getAttributeNames() {
        if (treeAttributes == null)
            return null;
        else
            return treeAttributes.getAttributeNames();
    }

    /**
     * @return a string containing a newick representation of the tree
     */
    public final String getNewick(int partition) {
        return Tree.Utils.newick(new ARGTree(this, partition));
        //return Tree.Utils.newick(this);
    }

    public void checkBranchSanity() {
        boolean plotted = false;
        for (Node node : nodes) {
            if (!node.isRoot()) {
                double length1 = 0;
                double length2 = 0;
                if (node.leftParent != null)
                    length1 = getNodeHeight(node.leftParent) - getNodeHeight(node);
                if (node.rightParent != null)
                    length2 = getNodeHeight(node.rightParent) - getNodeHeight(node);
                if (String.valueOf(length1).equals("NaN") || String.valueOf(length2).equals("NaN")) {
                    if (!plotted) {
                        System.err.println(toGraphString());
                        plotted = true;
                    }
                    System.err.println("Caught the NaN: node=" + node.number + " (" + node.getHeight() + ") lp=" + node.leftParent.number +
                            " (" + node.leftParent.getHeight() + ") rp=" + node.rightParent.number +
                            " (" + node.rightParent.getHeight() + ")");
                    System.exit(-1);
                }
            }
        }

    }


    /**
     * @return a string containing a newick representation of the tree
     */
    public String toString() {
        //System.err.println("toString");
        //System.err.println(this.maxNumberOfPartitions);
        //System.exit(-1);
        //ARGTree tree = new ARGTree(treeModel, partition);
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < maxNumberOfPartitions; i++) {
            sb.append(i + ": ");
            sb.append(getNewick(i));
        }
        //System.err.println(new String(sb));
        return new String(sb);

    }

    public static final String nullEdge = " -";


    public void appendGraphStringOld(StringBuffer sb) {
        int cnt = 0;

        for (Node node : nodes)
            node.number = cnt++;

        cnt = 0;
        for (Node node : nodes) {
            sb.append(cnt == 0 ? "[" : ",[");
            cnt++;
            sb.append(node.number + ":");

            if (node.leftParent == null)
                sb.append(nullEdge);
            else
                sb.append(" " + node.leftParent.number);

            if (node.rightParent == null)
                sb.append(nullEdge);
            else
                sb.append(" " + node.rightParent.number);

            if (node.leftChild == null)
                sb.append(nullEdge);
            else
                sb.append(" " + node.leftChild.number);

            if (node.rightChild == null)
                sb.append(nullEdge);
            else
                sb.append(" " + node.rightChild.number);
//			sb.append(" " + node.bifurcation);
            if (node.taxon != null)
                sb.append(" " + node.taxon.toString());
//            if (node.leftPartition != null)
//                sb.append(" l");
//            if (node.rightPartition != null)
//                sb.append(" r");
            sb.append("]");
//
//            );
        }
//        sb.append("Root = " + ((Node) getRoot()).number + "\n");
//		sb.append("\n");
    }

    public String toGraphString() {
        //	if( true )
        //		return null;
        int cnt = 1;
        for (Node node : nodes) {
            node.number = cnt;
            cnt++;
        }
        StringBuffer sb = new StringBuffer();
        sb.append("Total length: " + nodes.size() + "\n");
        for (Node node : nodes) {
            sb.append(node.number + ":");

            if (node.leftParent == null)
                sb.append(" 0");
            else
                sb.append(" " + node.leftParent.number);

            if (node.rightParent == null)
                sb.append(" 0");
            else
                sb.append(" " + node.rightParent.number);

            if (node.leftChild == null)
                sb.append(" 0");
            else
                sb.append(" " + node.leftChild.number);

            if (node.rightChild == null)
                sb.append(" 0");
            else
                sb.append(" " + node.rightChild.number);
//			sb.append(" " + node.bifurcation);
            if (node.taxon != null)
                sb.append(" " + node.taxon.toString());
            if (node.partitioning != null)
                sb.append(" p");
            /*		if (node.leftPartition != null)
                   sb.append(" l");
               if (node.rightPartition != null)
                   sb.append(" r");*/
            sb.append("\t" + getNodeHeight(node));
            sb.append("\n");
        }
        sb.append("Root = " + ((Node) getRoot()).number + "\n");
        return new String(sb);
    }


    public Tree getCopy() {
        throw new UnsupportedOperationException("please don't call this function");
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document document) {
        throw new RuntimeException("Not implemented yet");
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TREE_MODEL;
        }

        /**
         * @return a tree object based on the XML element it was passed.
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Tree tree = (Tree) xo.getChild(Tree.class);
            ARGModel treeModel = new ARGModel(tree);

            Logger.getLogger("dr.evomodel").info("Creating the tree model, '" + xo.getId() + "'");

            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof XMLObject) {

                    XMLObject cxo = (XMLObject) xo.getChild(i);

                    if (cxo.getName().equals(ROOT_HEIGHT)) {

                        replaceParameter(cxo, treeModel.getRootHeightParameter());

                    } else if (cxo.getName().equals(LEAF_HEIGHT)) {

                        String taxonName;
                        if (cxo.hasAttribute(TAXON)) {
                            taxonName = cxo.getStringAttribute(TAXON);
                        } else {
                            throw new XMLParseException("taxa element missing from leafHeight element in treeModel element");
                        }

                        int index = treeModel.getTaxonIndex(taxonName);
                        if (index == -1) {
                            throw new XMLParseException("taxon " + taxonName + " not found for leafHeight element in treeModel element");
                        }
                        NodeRef node = treeModel.getExternalNode(index);
                        replaceParameter(cxo, treeModel.getLeafHeightParameter(node));

                    } else if (cxo.getName().equals(NODE_HEIGHTS)) {

                        boolean rootNode = false;
                        boolean internalNodes = false;
                        boolean leafNodes = false;

                        if (cxo.hasAttribute(ROOT_NODE)) {
                            rootNode = cxo.getBooleanAttribute(ROOT_NODE);
                        }

                        if (cxo.hasAttribute(INTERNAL_NODES)) {
                            internalNodes = cxo.getBooleanAttribute(INTERNAL_NODES);
                        }

                        if (cxo.hasAttribute(LEAF_NODES)) {
                            leafNodes = cxo.getBooleanAttribute(LEAF_NODES);
                        }

                        if (!rootNode && !internalNodes && !leafNodes) {
                            throw new XMLParseException("one or more of root, internal or leaf nodes must be selected for the nodeHeights element");
                        }

                        replaceParameter(cxo, treeModel.createNodeHeightsParameter(rootNode, internalNodes, leafNodes));

                    } else if (cxo.getName().equals(NODE_RATES)) {

                        boolean rootNode = false;
                        boolean internalNodes = false;
                        boolean leafNodes = false;

                        if (cxo.hasAttribute(ROOT_NODE)) {
                            rootNode = cxo.getBooleanAttribute(ROOT_NODE);
                        }

                        if (cxo.hasAttribute(INTERNAL_NODES)) {
                            internalNodes = cxo.getBooleanAttribute(INTERNAL_NODES);
                        }

                        if (cxo.hasAttribute(LEAF_NODES)) {
                            leafNodes = cxo.getBooleanAttribute(LEAF_NODES);
                        }

                        //if (rootNode) {
                        //	throw new XMLParseException("root node does not have a rate parameter");
                        //}

                        if (!rootNode && !internalNodes && !leafNodes) {
                            throw new XMLParseException("one or more of root, internal or leaf nodes must be selected for the nodeRates element");
                        }

                        replaceParameter(cxo, treeModel.createNodeRatesParameter(rootNode, internalNodes, leafNodes));

                    } else if (cxo.getName().equals(NODE_TRAITS)) {

                        boolean rootNode = false;
                        boolean internalNodes = false;
                        boolean leafNodes = false;

                        if (cxo.hasAttribute(ROOT_NODE)) {
                            rootNode = cxo.getBooleanAttribute(ROOT_NODE);
                        }

                        if (cxo.hasAttribute(INTERNAL_NODES)) {
                            internalNodes = cxo.getBooleanAttribute(INTERNAL_NODES);
                        }

                        if (cxo.hasAttribute(LEAF_NODES)) {
                            leafNodes = cxo.getBooleanAttribute(LEAF_NODES);
                        }

                        if (!rootNode && !internalNodes && !leafNodes) {
                            throw new XMLParseException("one or more of root, internal or leaf nodes must be selected for the nodeTraits element");
                        }

                        replaceParameter(cxo, treeModel.createNodeTraitsParameter(rootNode, internalNodes, leafNodes));

                    } else {
                        throw new XMLParseException("illegal child element in " + getParserName() + ": " + cxo.getName());
                    }

                } else if (xo.getChild(i) instanceof Tree) {
                    // do nothing - already handled
                } else {
                    throw new XMLParseException("illegal child element in  " + getParserName() + ": " + xo.getChildName(i) + " " + xo.getChild(i));
                }
            }

            treeModel.setupHeightBounds();

            Logger.getLogger("dr.evomodel").info("  initial tree topology = " + Tree.Utils.uniqueNewick(treeModel, treeModel.getRoot()));
            return treeModel;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a model of the tree. The tree model includes and attributes of the nodes " +
                    "including the age (or <i>height</i>) and the rate of evolution at each node in the tree.";
        }

        public String getExample() {
            return
                    "<!-- the tree model as special sockets for attaching parameters to various aspects of the tree     -->\n" +
                            "<!-- The treeModel below shows the standard setup with a parameter associated with the root height -->\n" +
                            "<!-- a parameter associated with the internal node heights (minus the root height) and             -->\n" +
                            "<!-- a parameter associates with all the internal node heights                                     -->\n" +
                            "<!-- Notice that these parameters are overlapping                                                  -->\n" +
                            "<!-- The parameters are subsequently used in operators to propose changes to the tree node heights -->\n" +
                            "<treeModel id=\"treeModel1\">\n" +
                            "	<tree idref=\"startingTree\"/>\n" +
                            "	<rootHeight>\n" +
                            "		<parameter id=\"treeModel1.rootHeight\"/>\n" +
                            "	</rootHeight>\n" +
                            "	<nodeHeights internalNodes=\"true\" rootNode=\"false\">\n" +
                            "		<parameter id=\"treeModel1.internalNodeHeights\"/>\n" +
                            "	</nodeHeights>\n" +
                            "	<nodeHeights internalNodes=\"true\" rootNode=\"true\">\n" +
                            "		<parameter id=\"treeModel1.allInternalNodeHeights\"/>\n" +
                            "	</nodeHeights>\n" +
                            "</treeModel>";

        }

        public Class getReturnType() {
            return ARGModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(Tree.class),
                new ElementRule(ROOT_HEIGHT, Parameter.class, "A parameter definition with id only (cannot be a reference!)", false),
                new ElementRule(NODE_HEIGHTS,
                        new XMLSyntaxRule[]{
                                AttributeRule.newBooleanRule(ROOT_NODE, true, "If true the root height is included in the parameter"),
                                AttributeRule.newBooleanRule(INTERNAL_NODES, true, "If true the internal node heights (minus the root) are included in the parameter"),
                                new ElementRule(Parameter.class, "A parameter definition with id only (cannot be a reference!)")
                        }, 1, Integer.MAX_VALUE)
        };

        public Parameter getParameter(XMLObject xo) throws XMLParseException {

            int paramCount = 0;
            Parameter param = null;
            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof Parameter) {
                    param = (Parameter) xo.getChild(i);
                    paramCount += 1;
                }
            }

            if (paramCount == 0) {
                throw new XMLParseException("no parameter element in treeModel " + xo.getName() + " element");
            } else if (paramCount > 1) {
                throw new XMLParseException("More than one parameter element in treeModel " + xo.getName() + " element");
            }


            return param;
        }

        public void replaceParameter(XMLObject xo, Parameter newParam) throws XMLParseException {

            for (int i = 0; i < xo.getChildCount(); i++) {

                if (xo.getChild(i) instanceof Parameter) {

                    XMLObject rxo = null;
                    Object obj = xo.getRawChild(i);

                    if (obj instanceof Reference) {
                        rxo = ((Reference) obj).getReferenceObject();
                    } else if (obj instanceof XMLObject) {
                        rxo = (XMLObject) obj;
                    } else {
                        throw new XMLParseException("object reference not available");
                    }

                    if (rxo.getChildCount() > 0) {
                        throw new XMLParseException("No child elements allowed in parameter element.");
                    }

                    if (rxo.hasAttribute(XMLParser.IDREF)) {
                        throw new XMLParseException("References to " + xo.getName() + " parameters are not allowed in treeModel.");
                    }

                    if (rxo.hasAttribute(ParameterParser.VALUE)) {
                        throw new XMLParseException("Parameters in " + xo.getName() + " have values set automatically.");
                    }

                    if (rxo.hasAttribute(ParameterParser.UPPER)) {
                        throw new XMLParseException("Parameters in " + xo.getName() + " have bounds set automatically.");
                    }

                    if (rxo.hasAttribute(ParameterParser.LOWER)) {
                        throw new XMLParseException("Parameters in " + xo.getName() + " have bounds set automatically.");
                    }

                    if (rxo.hasAttribute(XMLParser.ID)) {

                        newParam.setId(rxo.getStringAttribute(XMLParser.ID));
                    }

                    rxo.setNativeObject(newParam);

                    return;
                }
            }
        }
    };

    // ***********************************************************************
    // Private methods
    // ***********************************************************************

    /**
     * @return the node that this parameter is a member of
     */
    protected Node getNodeOfParameter(Parameter parameter) {

        if (parameter == null) throw new IllegalArgumentException("Parameter is null!");

        // for (int i =0; i < nodes.length; i++) {
        for (Node node : nodes) {
            if (node.heightParameter == parameter) {
                return node;
            }
            if (hasRates && node.rateParameter == parameter) {
                return node;
            }
            if (hasTraits && node.traitParameter == parameter) {
                return node;
            }
        }
        throw new RuntimeException("Parameter not found in any nodes:" + parameter.getId());
    }

    /**
     * Get the root height parameter. Is private because it can only be called by the XMLParser
     */
    protected Parameter getRootHeightParameter() {

        return root.heightParameter;
    }

    /**
     * @return the relevant node height parameter. Is private because it can only be called by the XMLParser
     */
    protected Parameter createNodeHeightsParameter(boolean rootNode, boolean internalNodes, boolean leafNodes) {

        if (!rootNode && !internalNodes && !leafNodes) {
            throw new IllegalArgumentException("At least one of rootNode, internalNodes or leafNodes must be true");
        }

        VariableSizeCompoundParameter parameter = new VariableSizeCompoundParameter("nodeHeights");
//        System.err.println("Constructed nodeHeights");
        for (int i = externalNodeCount; i < nodeCount; i++) {
            Node node = nodes.get(i);
            if ((rootNode && node == root) || (internalNodes && node != root)) {
                parameter.addParameter(node.heightParameter);
            }
        }

        if (leafNodes) {
            for (int i = 0; i < externalNodeCount; i++) {
                parameter.addParameter(nodes.get(i).heightParameter);
            }
        }

        return parameter;
    }

    protected Parameter getLeafHeightParameter(NodeRef nr) {
        Node node = (Node) nr;
        if (!isExternal(node)) {
            throw new RuntimeException("only root and leaves can be used with setNodeHeightParameter");
        }

        return nodes.get(nodes.indexOf(node)).heightParameter;
    }

    /**
     * @return the relevant node rate parameter. Is private because it can only be called by the XMLParser
     */
    protected Parameter createNodeRatesParameter(boolean rootNode, boolean internalNodes, boolean leafNodes) {

        if (!rootNode && !internalNodes && !leafNodes) {
            throw new IllegalArgumentException("At least one of rootNode, internalNodes or leafNodes must be true");
        }

        CompoundParameter parameter = new CompoundParameter("nodeRates");

        hasRates = true;

        for (int i = externalNodeCount; i < nodeCount; i++) {
            Node node = nodes.get(i);
            node.createRateParameter();
            if ((rootNode && node == root) || (internalNodes && node != root)) {
                parameter.addParameter(node.rateParameter);
            }
        }

        for (int i = 0; i < externalNodeCount; i++) {
            Node node = nodes.get(i);
            node.createRateParameter();
            if (leafNodes) {
                parameter.addParameter(node.rateParameter);
            }
        }

        return parameter;
    }

    /**
     * Create a node traits parameter. Is private because it can only be called by the XMLParser
     */
    public Parameter createNodeTraitsParameter(boolean rootNode, boolean internalNodes, boolean leafNodes) {

        if (!rootNode && !internalNodes && !leafNodes) {
            throw new IllegalArgumentException("At least one of rootNode, internalNodes or leafNodes must be true");
        }

        CompoundParameter parameter = new CompoundParameter("nodeTraits");

        hasTraits = true;

        for (int i = externalNodeCount; i < nodeCount; i++) {
            Node node = nodes.get(i);
            node.createTraitParameter();
            if ((rootNode && node == root) || (internalNodes && node != root)) {
                parameter.addParameter(node.traitParameter);
            }
        }

        for (int i = 0; i < externalNodeCount; i++) {
            Node node = nodes.get(i);
            node.createTraitParameter();
            if (leafNodes) {
                parameter.addParameter(node.traitParameter);
            }
        }

        return parameter;
    }

    /**
     * This method swaps the parameter objects of the two nodes
     * but maintains the values in each node.
     * This method is used to ensure that root node of the tree
     * always has the same parameter object.
     */
    private void swapParameterObjects(Node n1, Node n2) {

        double height1 = n1.getHeight();
        double height2 = n2.getHeight();

        double rate1 = 1.0, rate2 = 1.0;
        double trait1 = 0.0, trait2 = 0.0;

        if (hasRates) {
            rate1 = n1.getRate();
            rate2 = n2.getRate();
        }

        if (hasTraits) {
            trait1 = n1.getTrait();
            trait2 = n2.getTrait();
        }

        Parameter temp = n1.heightParameter;
        n1.heightParameter = n2.heightParameter;
        n2.heightParameter = temp;

        if (hasRates) {
            temp = n1.rateParameter;
            n1.rateParameter = n2.rateParameter;
            n2.rateParameter = temp;
        }

        if (hasTraits) {
            temp = n1.traitParameter;
            n1.traitParameter = n2.traitParameter;
            n2.traitParameter = temp;
        }

        n1.heightParameter.setParameterValueQuietly(0, height1);
        n2.heightParameter.setParameterValueQuietly(0, height2);

        if (hasRates) {
            n1.rateParameter.setParameterValueQuietly(0, rate1);
            n2.rateParameter.setParameterValueQuietly(0, rate2);
        }

        if (hasTraits) {
            n1.traitParameter.setParameterValueQuietly(0, trait1);
            n2.traitParameter.setParameterValueQuietly(0, trait2);
        }

    }


    // **************************************************************
    // Private inner classes
    // **************************************************************
    public class Node implements NodeRef {

        public Node leftParent, rightParent;
        public Node leftChild, rightChild;
        public int number;
        public Parameter heightParameter;
        public Parameter rateParameter = null;
        public Parameter traitParameter = null;
        public Taxon taxon = null;

//		public BitSet leftPartition = null;
//		public BitSet rightPartition = null;
        //       public Node dupSister = null;
        //       public Node linkSister = null;
        //       public Node dupParent = null;

        //       public Node leftParent;
        //       public Node rightParent;

        //       public int leftPartition;
        //       public int rightPartition;

        public Parameter partitioning;

        public boolean bifurcation = true;

        public int countReassortmentChild(Tree tree) {
            //int cnt = 0;
            if (isExternal())
                return 0;
//        	if( leftChild == null ) {
//        		System.err.println("left is null");
//        		System.err.println("is reassort = "+isReassortment());
//        		System.err.println("right child = "+Tree.Utils.uniqueNewick(tree, rightChild));
//        	}
//        	if( rightChild == null ) {
//        		System.err.println("right is null");
//        		System.err.println("is reassort = "+isReassortment());
//        		System.err.println("left child = "+Tree.Utils.uniqueNewick(tree, leftChild));
//        	}
            if (isReassortment()) {
                return 1 + leftChild.countReassortmentChild(tree);
            } else if (isBifurcationDoublyLinked()) {
                return leftChild.countReassortmentChild(tree);
            } else {
                return leftChild.countReassortmentChild(tree) + rightChild.countReassortmentChild(tree);
            }
        }

        public Node() {
            leftParent = rightParent = null;
            leftChild = rightChild = null;
            heightParameter = null;
            number = 0;
            taxon = null;
//			leftPartition = null;
//			rightPartition = null;
            partitioning = null;
        }

        /**
         * Constructor to build an ARG into a bifurcating tree
         *
         * @param node
         */
        public Node(Node node) {
            leftParent = rightParent = null;
            leftChild = rightChild = null;
            heightParameter = node.heightParameter;
            taxon = node.taxon;
            number = node.number;
            //	if( node.isReassortment() ) {
            //Node parent = leftParent;
            //parent.removeChild(this);
            //	return new Node(leftChild);
            //	} else {
            if (node.leftChild != null) {
                if (node.leftChild.isReassortment())
                    singleAddChild(new Node(node.leftChild.leftChild));
                else
                    singleAddChild(new Node(node.leftChild));
            }
            if (node.rightChild != null) {
                if (node.rightChild.isReassortment())
                    singleAddChild(new Node(node.rightChild.leftChild));
                else
                    singleAddChild(new Node(node.rightChild));
            }
            //	}
        }

        /**
         * constructor used to clone a node and all children with no reassortments
         */
        public Node(Tree tree, NodeRef node) {
            leftParent = rightParent = null;
            leftChild = rightChild = null;
            //leftPartition = new BitSet();
            //leftPartition.set(0);
//			leftPartition = null;
            //rightPartition = new BitSet();
            //rightPartition.set(0);
//			rightPartition = null;
            partitioning = null;
            heightParameter = new Parameter.Default(tree.getNodeHeight(node));
            addParameter(heightParameter);

            number = node.getNumber();
            taxon = tree.getNodeTaxon(node);

            for (int i = 0; i < tree.getChildCount(node); i++) {
                singleAddChild(new Node(tree, tree.getChild(node, i)));
            }
//            System.err.println("Built initial tree");
            //System.exit(-1);
        }

//       public Node(Node node, int partition)
//       {
//            leftParent = rightParent = null;
//            leftChild = rightChild = null;
//            leftPartition = node.leftPartition;
//            rightPartition = node.rightPartition;
//            //leftPartition = rightPartition = null;
//            heightParameter = node.heightParameter;
//            number = node.number;
//            taxon = node.taxon;
//            bifurcation = node.bifurcation;
//  //          System.err.println("Examinging "+number);
//            Node left;
//            if( node.leftChild != null ) {
//            	Node left = node.getChild(0,partition);
//            }
//            Node right;
//            if( node.rightChild != null ) {
//            	Node right = node.getChild(1,partition);
//            }
//            
//            	if( left != null ) {
//  //          		System.err.println("Adding "+number+"->"+left.number);
//            		singleAddChild(new Node(left,partition));
//            	}
//            }
//           
//            	if( right != null ) {
//  //          		System.err.println("Adding "+number+"->"+right.number);
//            		singleAddChild(new Node(right,partition));
//            	}
//            }
//        }

        public Node(Node inode, int partition, ArrayList<Node> nodes) {
            leftParent = rightParent = null;
            leftChild = rightChild = null;
            Node node = inode;
            while (node.isBifurcationDoublyLinked()) {
                node = node.leftChild.leftChild;
                //        	System.err.println("Does this do anything?");
            }
            //else
            //	node = inode;
            //if( node.bifurcation ) {
            //leftPartition = node.leftPartition;
            //rightPartition = node.rightPartition;
            //leftPartition = rightPartition = null;
            heightParameter = node.heightParameter;
            number = node.number;
            //rightParent = leftParent = inode.leftParent;
            //rightParent = inode;
            //nodes.add(this);
            //number = nodes.size();
            taxon = node.taxon;
            bifurcation = true;
            //          System.err.println("Examinging "+number);
            if (node.isExternal())
                nodes.add(this);
            else {
                Node left, right;
                left = node.getChild(0, partition);
                right = node.getChild(1, partition);
                if (left != null || right != null) {
                    if (left != null) {
                        //leftChild = new Node(left,partition,nodes);
                        //this.leftChild
                        //          		System.err.println("Adding "+number+"->"+left.number);
                        singleAddChild(new Node(left, partition, nodes));
                    }
                    //}}
                    //f( node.rightChild != null ) {
                    //Node right = node.getChild(1,partition);
                    if (right != null) {
                        //          		System.err.println("Adding "+number+"->"+right.number);
                        singleAddChild(new Node(right, partition, nodes));
                        //rightChild = new Node(right,partition,nodes);
                    }
                    //rightParent = leftParent = inode.leftParent;
                    nodes.add(this);
                }
            }
        }

        /*		public void setPartitionRecursively(int partition) {
              if (leftChild != null) {
                  if (leftPartition != null)
                      leftPartition.set(partition);
                  leftChild.setPartitionRecursively(partition);
              }
              if (rightChild != null) {
                  if (leftPartition != null)
                      rightPartition.set(partition);
                  rightChild.setPartitionRecursively(partition);
              }

          }*/
        public void setPartitionRecursively(int partition) {
            boolean onLeft = MathUtils.nextBoolean();
            if (leftChild != null) {
                //if( leftPartition != null )
                //	leftPartition.set(partition);
                if (partitioning != null && onLeft)
                    partitioning.setParameterValue(partition, 1);
                leftChild.setPartitionRecursively(partition);
            }
            if (rightChild != null) {
                //if( leftPartition != null )
                //	rightPartition.set(partition);
                if (partitioning != null && !onLeft)
                    partitioning.setParameterValue(partition, 1);
                rightChild.setPartitionRecursively(partition);
            }

        }

        /*      public Node findPartitionTreeRoot(int partition) {
              if (leftPartition.get(partition) && rightPartition.get(partition))
                  return this;
              if (leftPartition.get(partition))
                  return leftChild.findPartitionTreeRoot(partition);
              if (rightPartition.get(partition))
                  return rightChild.findPartitionTreeRoot(partition);
              throw new IllegalArgumentException("Partition " + partition + " never found in " + this.toString());
          }*/

//        private boolean doesBifurcate(int partition) {
//        	if( leftChild.parent != null ) {
//        		if( rightChild.parent != null )
//        			return true;
//        		if( (rightChild.leftParent == this) ||
//        			(rightChild.rightParent == this) )
//        			return true;
//        	}
//        	if( rightChild.parent != null ) {
//        		if( (leftChild.leftParent == this) ||
//        			(leftChild.rightParent == this) )
//        			return true;
//        	}
//        	return false;
//        }
//        
//        private boolean isRecombinantParent(Node parent) {
//        	if( leftParent == parent || rightParent == parent )
//        		return true;
//        	return false;
//        }
//        
//        public boolean isBifurcatingOrExternal(int partition) {
//        	// TODO protect against null errors at tips
//        	if( isExternal() )
//        		return true;
//        	if( leftChild.parent !=null && rightChild.parent !=null ) // standard case
//        		return true;
//        	if( leftChild.parent !=null &&
//        			(	(rightChild.leftParent == this && rightChild.leftPartition == partition) ||
//        				(rightChild.rightParent == this && rightChild.rightPartition == partition)
//        			) ) return true;
//        	if( rightChild.parent !=null &&
//        			(	(leftChild.leftParent == this && leftChild.leftPartition == partition) ||
//        				(leftChild.rightParent == this && leftChild.rightPartition == partition)
//        			) ) return true;
//        	return false;
//        }
//        
//        public Node straightDescendent(int partition) {
//        	// TODO protect against two direct children
//        	if( leftChild.parent != null )
//        		return leftChild;
//        	else {
//        		if( (leftChild.leftParent == this && leftChild.leftPartition == partition ) ||
//        			(leftChild.rightParent == this && leftChild.rightPartition == partition) )
//        			return leftChild;
//        	}
//        	if( rightChild.parent != null )
//        		return rightChild;
//        	else {
//        		if( (rightChild.leftParent == this && rightChild.leftPartition == partition) ||
//        		    (rightChild.rightParent == this && rightChild.rightPartition == partition) )
//        			return rightChild;
//        	}
//        	throw new IllegalArgumentException("No straight descendent found.");
//        }
//        
//        /** constructor used to clone a node and all children for a particular partition */
//        public Node(Node n, int partition) {
//        	Node node = n;
//        	parent = leftParent = rightParent = null;
//        	leftChild = rightChild = null;
//        	if( !node.isBifurcatingOrExternal(partition) )
//        		node = node.straightDescendent(partition);
//        	heightParameter = node.heightParameter;
//        	taxon = node.taxon;
//        	//boolean tip = true;
//        	Node lc = node.leftChild;
//        	Node rc = node.rightChild;
//         	if( lc != null ) {
//        		// LeftChild exists.  
//        		// Does lc bifurcate or do we skip for this partition
//  //      		if( lc.isBifurcatingOrExternal(partition) ) {
//        			addChild(new Node(lc, partition));
//        			//System.err.println(Tree.Utils.newick(lc)+" is bifurcating with P = "+partition);
//  //      		} else {
//  //      			addChild(new Node(lc.straightDescendent(partition),partition));
//        			//System.err.println(lc.straightDescendent(partition)+" is descendent.");
//  //      		}
//         	}
//        	if( rc != null ) {
//           		// RightChild exists.  
//        		// Does rc bifurcate or do we skip for this partition
//  //      		if( rc.isBifurcatingOrExternal(partition) ) {
//        			addChild(new Node(rc, partition));
//        			//System.err.println(rc.toString()+" is bifurcating with P = "+partition);
//  //      		} else {
//  //      			addChild(new Node(rc.straightDescendent(partition),partition));
//        			//System.err.println(rc.straightDescendent(partition)+" is descendent.");
//  //      		}
//         	}
//  		
//        }
//        
//        
//        /** constructor used to clone a subtree without duplicating height parameters
//         * 
//         *
//         */
//        public Node(Tree tree, Node node, int[] bits) {
//        	parent = null;
//        	leftChild = rightChild = null;
//        	heightParameter = node.heightParameter;
//        	linkSister = node;
//        	linkSister.linkSister = this;
//        	//for(int i=0; i < tree.getChildCount())
//        	boolean tip = true;
//        	if( node.leftChild != null ) {
//        		addChild(new Node(tree, node.leftChild, bits));
//        		tip = false;
//        	}
//        	if( node.rightChild != null ) {
//        		addChild(new Node(tree, node.rightChild, bits));
//        		tip = false;
//        	}
//        	if( tip ) {
//        		taxon = node.taxon;
//        		partitionSet = new BitSet();
//        		int len = bits.length;
//        		for(int i=0; i<len; i++) 
//        			partitionSet.set(bits[i]);
//        	}
//        }
//        
//        public final void setDupParent(Node parent) {
//        	this.dupParent = parent;
//        	if( leftChild != null )
//        		leftChild.setDupParent(parent);
//        	if( rightChild != null )
//        		rightChild.setDupParent(parent);
//        }
//        
//        public final void clearLinkSister() {
//        	this.linkSister = null;
//        	if( leftChild != null )
//        		leftChild.clearLinkSister();
//        	if( rightChild != null )
//        		rightChild.clearLinkSister();
//        }
//        
//        public final void clearDupParent() {
//        	this.dupParent = null;
//        	if( leftChild != null )
//        		leftChild.clearDupParent();
//        	if( rightChild != null )
//        		rightChild.clearDupParent();
//        }
//        

/*		public final void recursiveSetPartition(int partition) {
			leftPartition.set(partition);
			rightPartition.set(partition);
			if (leftChild != null)
				leftChild.recursiveSetPartition(partition);
			if (rightChild != null)
				rightChild.recursiveSetPartition(partition);
		}*/


        public void stripOutDeadEnds() {
            if (leftChild != null)
                leftChild.stripOutDeadEnds();
            if (rightChild != null && rightChild != leftChild)
                rightChild.stripOutDeadEnds();
            if (taxon == null && leftChild == null && rightChild == null)
                leftParent.doubleRemoveChild(this);
        }

        public Node stripOutSingleChildNodes(Node cRoot) {
            //Node rtn = cRoot;
            int childCount = getChildCount();
            if (childCount == 0) {
                return cRoot;
            }
            if (childCount == 2) {
                if (hasEqualChildren()) {
                    return leftChild.stripOutSingleChildNodes(leftChild);
                }
                leftChild.stripOutSingleChildNodes(cRoot);
                rightChild.stripOutSingleChildNodes(cRoot);
                return cRoot;
            }
            if (isRoot()) {
                if (leftChild != null) {
                    leftChild.leftParent = leftChild.rightParent = null;
                    return leftChild.stripOutSingleChildNodes(leftChild);
                } else {
                    rightChild.leftParent = rightChild.rightParent = null;
                    return rightChild.stripOutSingleChildNodes(rightChild);
                }
            }
            //       	System.err.println("Unlinking "+number);
            Node parent = leftParent;
            Node child = leftChild;
            if (child == null)
                child = rightChild;
            parent.doubleRemoveChild(this);
            doubleRemoveChild(child);
            parent.singleAddChild(child);
            return child.stripOutSingleChildNodes(cRoot);
        }

        public final void setupHeightBounds() {
            heightParameter.addBounds(new NodeHeightBounds(heightParameter));
        }

        public final void createRateParameter() {
            if (rateParameter == null) {
                rateParameter = new Parameter.Default(1.0);
                if (isRoot()) {
                    rateParameter.setId("root.rate");
                } else if (isExternal()) {
                    rateParameter.setId(getTaxonId(getNumber()) + ".rate");
                } else {
                    rateParameter.setId("node" + getNumber() + ".rate");
                }
                rateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
                addParameter(rateParameter);
            }
        }

        public final void createTraitParameter() {
            if (traitParameter == null) {
                traitParameter = new Parameter.Default(1.0);
                if (isRoot()) {
                    traitParameter.setId("root.trait");
                } else if (isExternal()) {
                    traitParameter.setId(getTaxonId(getNumber()) + ".trait");
                } else {
                    traitParameter.setId("node" + getNumber() + ".trait");
                }
                rateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

                addParameter(traitParameter);
            }
        }

        public final double getHeight() {
            return heightParameter.getParameterValue(0);
        }

        public final double getRate() {
            return rateParameter.getParameterValue(0);
        }

        public final double getTrait() {
            return traitParameter.getParameterValue(0);
        }

        public final void setHeight(double height) {
            heightParameter.setParameterValue(0, height);
        }

        public final void setRate(double rate) {
            //System.out.println("Rate set for parameter " + rateParameter.getParameterName());
            rateParameter.setParameterValue(0, rate);
        }

        public final void setTrait(double trait) {
            //System.out.println("Trait set for parameter " + traitParameter.getParameterName());
            traitParameter.setParameterValue(0, trait);
        }

        public int getNumber() {
            return number;
        }

        /**
         * Returns the number of children this node has.
         */
        public final int getChildCount() {
            int n = 0;
            if (leftChild != null) n++;
            if (rightChild != null && bifurcation) n++;
            return n;
        }

        public final int getParentCount() {
            int n = 0;
            if (leftParent != null) n++;
            if (rightParent != null) n++;
            return n;
        }

        public Node getChild(int n) {
            if (n == 0) return leftChild;
            if (n == 1) return rightChild;
            throw new IllegalArgumentException("ARGModel.Nodes can only have up to 2 children");
        }

//        private Node getChild(int n, int partition) {
//        	//if( isExternal() )
//        	//	return null;
//        	if (n == 0) { // Handle left side
//        		Node left = leftChild;
//        		if( left.isExternal() )
//        			return left;
//        		Node grandLeft = left.leftChild;
//        		Node grandRight = left.rightChild;
//        		boolean grandLeftValid = false;
//        		boolean grandRightValid = true;
//        		if( (grandLeft.leftParent == left && grandLeft.partitionSet.get(partition)) || 
//        			(grandLeft.rightParent == left && !grandLeft.partitionSet.get(partition)) )
//        			grandLeftValid = true;
//        		if( (grandRight.leftParent == left && grandRight.partitionSet.get(partition)) ||
//        			(grandRight.rightParent == left && !grandRight.partitionSet.get(partition)) )
//        			grandRightValid = true;
//        		if( grandLeftValid && grandRightValid )
//        			return left;
//        		if( grandLeftValid )
//        			return grandLeft;
//        		else
//        			return grandRight;
//        	}
//          	if (n == 1) { // Handle right side
//        		Node right = rightChild;
//        		if( right.isExternal() )
//        			return right;
//        		Node grandLeft = right.leftChild;
//        		Node grandRight = right.rightChild;
//        		boolean grandLeftValid = false;
//        		boolean grandRightValid = true;
//        		if( (grandLeft.leftParent == right && grandLeft.partitionSet.get(partition)) || 
//        			(grandLeft.rightParent == right && !grandLeft.partitionSet.get(partition)) )
//        			grandLeftValid = true;
//        		if( (grandRight.leftParent == right && grandRight.partitionSet.get(partition)) ||
//        			(grandRight.rightParent == right && !grandRight.partitionSet.get(partition)) )
//        			grandRightValid = true;
//        		if( grandLeftValid && grandRightValid )
//        			return right;
//        		if( grandLeftValid )
//        			return grandLeft;
//        		else
//        			return grandRight;
//        	}
//           	throw new IllegalArgumentException("ARGModel.Nodes can only have 2 children");      
//        }
//        

/*
		private boolean isRecombinantLinkedToMe(Node node, int partition) {
			if (node.leftParent == this && node.leftPartition.get(partition))
				return true;
			if (node.rightParent == this && node.rightPartition.get(partition))
				return true;
			return false;
		}
*/

        private boolean isBifurcationDoublyLinked() {
            return bifurcation && (leftChild == rightChild) && leftChild != null;
        }

        private boolean recombinantIsLinked(Node parent, int partition) {
            boolean left = leftParent == parent;
            boolean right = rightParent == parent;
            //       	System.err.println("Testing links on "+number);
//        	if( left && right ) {
//        		System.err.println("Doubly linked recombinant.");
//        		if( side == 0 ) {
//        			System.err.println("Going left");
//        			return true;
//        		}
//        		else {
//        			System.err.println("Not going right");
//        			return false;
//        		}
//        	//		return false;
//        	}

/*            if (leftPartition == null) {
                System.err.println("no left Partition");
                //System.exit(-1);
            }
            if (rightPartition == null) {
                System.err.println("no right Partition");
                //System.exit(-1);
            }
            if (left && leftPartition.get(partition))
                return true;
            if (right && rightPartition.get(partition))
                return true;

            return false;*/

            final double partitionSide = partitioning.getParameterValue(partition);
            if (left && partitionSide == 0)
                return true;
            if (right && partitionSide == 1)
                return true;
            return false;


        }


        public boolean checkForNullRights() {
            if (isExternal())
                return false;
            if (rightChild == null)
                return true;
            else
                return rightChild.checkForNullRights() || leftChild.checkForNullRights();
        }


        private Node findNextTreeNode(Node parent, int partition) { // searches down the ARG for the next bifurcation or tip
            if (isExternal())
                return this;
            Node next = this;
            while (next.isReassortment()) {
                if (recombinantIsLinked(parent, partition))
                    return leftChild.findNextTreeNode(parent, partition);
                else
                    return null;
            }
            //if( leftChild.findNextTreeMode())
            next = leftChild.findNextTreeNode(parent, partition);
            if (next == null)
                next = rightChild.findNextTreeNode(parent, partition);
            if (next == null) // TODO Error check for removal
                throw new IllegalArgumentException("Can't find next tree node.");
            return next;
        }


        private Node getChild(int n, int partition) { // Assuming an acyclic bifurcating tree
            Node child = null;
            if (n == 0)  // Handle left side
                child = leftChild;
            if (n == 1) // Handle right side
                child = rightChild;
            if (child.isExternal())
                return child;
            if (child.isBifurcationDoublyLinked()) {
                //System.err.println("Passing double from "+number+"->"+child.leftChild.number);
                return child.leftChild.getChild(0, partition);
            }
            if (child.isReassortment()) {
                if (child.recombinantIsLinked(this, partition))
                    return child.getChild(0, partition);
                else
                    return null;
            }
            if (child.leftChild.isReassortment() && child.rightChild.isReassortment()) {
                if (child.leftChild.recombinantIsLinked(child, partition))
                    return child;
                if (child.rightChild.recombinantIsLinked(child, partition))
                    return child;
                return null;
            }
            return child;
        }

//        		if( child.isExternal() )
//        			return child;
//        		//if( left.isReassortment() )
//        		//while( left.isReassortment() )
//        		//	left = left.leftChild; // Just pass through
//        		//if( left.isExternal() )
//        		//	return left;
//        		//if( left.isReassortment() )
//        		//	;
//        		Node grandLeft = child.leftChild.findNextTreeNode(child.leftChild, partition);
//        		Node grandRight = child.rightChild;
//        		boolean grandLeftValid = false;
//        		boolean grandRightValid = false;
//        		while( ! (grandLeftValid || grandRightValid) ) {
//        			
//        		}
//        			left.leftPartition.get(partition);
//        		boolean grandRightValid = left.rightPartition.get(partition);
//        		if( grandLeftValid && grandRightValid )
//        			return left;
//        		if( grandLeftValid )
//        			return left.leftChild;
//        		if( grandRightValid )
//        			return left.rightChild;
//        		throw new IllegalArgumentException("Partition not found");
//        	}
//          	if (n == 1) { // Handle right side
//          		if( !rightPartition.get(partition) )
//          			throw new IllegalArgumentException("Partition not found");
//        		Node right = rightChild;
//        		if( right.isExternal() )
//        			return right;
//        		boolean grandLeftValid = right.leftPartition.get(partition);
//        		boolean grandRightValid = right.rightPartition.get(partition);
//        		if( grandLeftValid && grandRightValid )
//        			return right;
//        		if( grandLeftValid )
//        			return right.leftChild;
//        		if( grandRightValid )
//        			return right.rightChild;
//        		throw new IllegalArgumentException("Partition not found");
//        	}
//           	throw new IllegalArgumentException("ARGModel.Nodes can only have 2 children");      
//        }

        public Node getParent(int n) {
            if (n == 0) return leftParent;
            if (n == 1) return rightParent;
            throw new IllegalArgumentException("ARGModel.Nodes can only have 2 parents");
        }

        public boolean hasChild(Node node) {
            return (leftChild == node || rightChild == node);
        }

        /**
         * add new child node
         *
         * @param node new child node
         */
        public void singleAddChild(Node node) {
            if (leftChild == null) {
                leftChild = node;
            } else if (rightChild == null) {
                rightChild = node;
            } else {
                throw new IllegalArgumentException("ARGModel.Nodes can only have 2 children");
            }
            //if( node.leftParent == null )
            if (node.leftParent == null)
                node.leftParent = this;
            if (node.rightParent == null)
                node.rightParent = this;
        }

        public void singleAddChildWithOneParent(Node node) {
            if (leftChild == null) {
                leftChild = node;
            } else if (rightChild == null) {
                rightChild = node;
            } else {
                throw new IllegalArgumentException("ARGModel.Node " + number + " can only have 2 children");
            }
            //if( node.leftParent == null )
            if (node.leftParent == null) {
                node.leftParent = this;
            } else if (node.rightParent == null) {
                node.rightParent = this;
            } else {
                throw new IllegalArgumentException("ARGModel.Nodes can only have 2 parents");
            }
        }

        public void doubleAddChild(Node node) {
            if (leftChild == null) {
                leftChild = node;
            }
            if (rightChild == null) {
                rightChild = node;
            }// else {
            //  throw new IllegalArgumentException("ARGModel.Nodes can only have 2 children");
            //}
            //if( node.leftParent == null )
            if (node.leftParent == null)
                node.leftParent = this;
            if (node.rightParent == null)
                node.rightParent = this;
        }

        public void doubleAddChildWithOneParent(Node node) {
            if (leftChild == null) {
                leftChild = node;
            }
            if (rightChild == null) {
                rightChild = node;
            }// else {
            //  throw new IllegalArgumentException("ARGModel.Nodes can only have 2 children");
            //}
            if (node.leftParent == null) {
                node.leftParent = this;
            } else if (node.rightParent == null) {
                node.rightParent = this;
            } else {
                throw new IllegalArgumentException("ARGModel.Nodes can only have 2 parents");
            }
        }

        public void addChildNoParentConnection(Node node) {
            if (leftChild == null)
                leftChild = node;
            else if (rightChild == null)
                rightChild = node;
            else
                throw new IllegalArgumentException("Nodes can only have 2 children.");
        }

        //public void addChild()

        //		public void addChildRecombinant(Node node, BitSet partition) {
        public void addChildRecombinant(Node node, Parameter partition) {
            //if( leftChild == null && rightChild == null ) {
            //	leftChild = rightChild = node;
            //} else
            //if( leftChild == null && rightChild == null ) {
            //	System.err.println("yep");
            //	//System.exit(-1);
            //	leftChild = rightChild = node;
            //}
            if (leftChild == null)
                leftChild = node;
            //leftPartition = partition;
            if (rightChild == null)
                rightChild = node;
            //rightPartition = partition;
            //} else {
            //	throw new IllegalArgumentException("Nodes can only have 2 children.");
            //}
//        	node.parent = null;
            if (node.leftParent == null) {
                node.leftParent = this;
//				node.leftPartition = partition;
                node.partitioning = partition;
            } else if (node.rightParent == null) {
                node.rightParent = this;
//				node.rightPartition = partition;
                node.partitioning = partition;
            } else {
                throw new IllegalArgumentException("Recombinant nodes can only have 2 parents.");
            }
        }

        /**
         * remove child
         *
         * @param node child to be removed
         */
        public Node doubleRemoveChild(Node node) {
            boolean found = false;
            if (leftChild == node) {
                leftChild = null;
                //leftPartition = null;
                found = true;
            }
            if (rightChild == node) {
                rightChild = null;
                //rightPartition = null;
                found = true;
            }
            if (!found)
                throw new IllegalArgumentException("Unknown child node");
            if (node.leftParent == this)
                node.leftParent = null;
            //node.leftPartition = null;
            if (node.rightParent == this)
                node.rightParent = null;
            //node.rightPartition = null;
            return node;
        }

        public Node singleRemoveChild(Node node) {
            //boolean found = false;
            if (leftChild == node) {
                leftChild = null;
                //leftPartition = null;
                //found = true;
            } else if (rightChild == node) {
                rightChild = null;
                //rightPartition = null;
                //found = true;
            }
            //if( !found )
            else
                throw new IllegalArgumentException("Unknown child node");
            if (node.bifurcation) {
                node.leftParent = node.rightParent = null;
                return null;
            }
            if (node.leftParent == this)
                node.leftParent = null;
                //node.leftPartition = null;
            else if (node.rightParent == this)
                node.rightParent = null;
            //node.rightPartition = null;
            return node;
        }

        /**
         * remove child
         *
         * @param n number of child to be removed
         */
        public Node removeChild(int n) {
            Node node;
            if (n == 0) {
                node = leftChild;
                leftChild = null;
            } else if (n == 1) {
                node = rightChild;
                rightChild = null;
            } else {
                throw new IllegalArgumentException("TreeModel.Nodes can only have 2 children");
            }
            if (node.leftParent == this)
                node.leftParent = null;
            if (node.rightParent == this)
                node.rightParent = null;
            return node;
        }

        public boolean hasChildren() {
            return (leftChild != null || rightChild != null);
        }

        public boolean isExternal() {
            return !hasChildren();
        }

        public boolean isRoot() {
            return (leftParent == null && rightParent == null);
        }

        public boolean hasEqualChildren() {
            return (leftChild == rightChild);
        }

        public boolean isBifurcation() {
            return bifurcation;
        }

        //public boolean isReassortment() { return hasChildren() && (leftChild == rightChild); }
        public boolean isReassortment() {
            return !bifurcation;
        }

        public String toString() {
            return taxon.getId();
        }
    }

    /**
     * This class provides bounds for parameters that represent a node height
     * in this tree model.
     */
    private class NodeHeightBounds implements Bounds {

        public NodeHeightBounds(Parameter parameter) {
            nodeHeightParameter = parameter;
        }

        public double getUpperLimit(int i) {
            // I think only upper bounds are of concern with linked subtrees
            // because everything below has only one parameter
            // TODO -- check this!
            Node node = getNodeOfParameter(nodeHeightParameter);
            // Returns the first node in nodes[] with this height parameter

            if (node.isRoot()) {
                return Double.POSITIVE_INFINITY;
//                return 10.0;
            } else {
                if (node.leftParent == null) {
                    System.err.println("leftParent of " + node.number + " is null");
                }
                if (node.rightParent == null) {
                    System.err.println("rightParent of " + node.number + " is null");
                }
                return Math.min(node.leftParent.getHeight(), node.rightParent.getHeight());
            }
        }

        public double getLowerLimit(int i) {

            Node node = getNodeOfParameter(nodeHeightParameter);
            // System.err.println("Is node recombinant? "+node.isReassortment());
            if (node.isExternal()) {
                return 0.0;
            } else {
                if (node.leftChild == null)
                    System.err.println("Node " + node.number + " has null leftChild");
                if (node.rightChild == null)
                    System.err.println("Node " + node.number + " has null rightChild");
                //System.err.println(node.number+" "+(node.leftChild==null)+" "+(node.rightChild==null));
                return Math.max(node.leftChild.getHeight(), node.rightChild.getHeight());
            }
        }

        public int getBoundsDimension() {
            return 1;
        }


        private Parameter nodeHeightParameter = null;
    }

    // ***********************************************************************
    // Private members
    // ***********************************************************************


    /**
     * root node
     */
    protected Node root = null;
    protected int storedRootNumber;

    /**
     * list of internal nodes (including root)
     */
    //protected Node[] nodes = null;
    //protected Node[] storedNodes = null;
    protected ArrayList<Node> nodes = null;
    protected ArrayList<Node> storedNodes = null;

    /**
     * number of nodes (including root and tips)
     */
    protected int nodeCount;
    protected int storedNodeCount;

    /**
     * number of external nodes
     */
    protected int externalNodeCount;

    /**
     * number of internal nodes (including root)
     */
    protected int internalNodeCount;
    protected int storedInternalNodeCount;

    /**
     * holds the units of the trees branches.
     */
    private int units = SUBSTITUTIONS;

    protected boolean inEdit = false;

    private boolean hasRates = false;
    private boolean hasTraits = false;

    protected Parameter[] addedParameters = null;
    protected Parameter[] removedParameters = null;
    protected Parameter addedPartitioningParameter = null;
    protected Parameter removedPartitioningParameter = null;
    protected VariableSizeCompoundParameter partitioningParameters;

    protected VariableSizeCompoundParameter storedInternalNodeHeights;
    protected VariableSizeCompoundParameter storedInternalAndRootNodeHeights;
    protected Node[] addedNodes = null;
    protected Node[] removedNodes = null;
}
