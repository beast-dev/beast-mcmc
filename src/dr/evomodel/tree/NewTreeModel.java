/*
 * TreeModel.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.tree;

import dr.evolution.tree.*;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.inference.model.*;
import dr.util.Attributable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.*;

/**
 * A model component for trees.
 * This is a cleaned up and simplified version of the original TreeModel and is intended to replace it.
 *
 * Before it can be used, considerable refactoring is required to switch from using node rates and traits
 * to TreeParameterModel
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeModel.java,v 1.129 2006/01/05 17:55:47 rambaut Exp $
 */
public class NewTreeModel extends AbstractModel implements MutableTree, Citable {

    //
    // Public stuff
    //

    public static final String TREE_MODEL = "treeModel";

    private static final boolean TEST_NODE_BOUNDS = false;

    public NewTreeModel(String name) {
        super(name);
        nodeCount = 0;
        externalNodeCount = 0;
        internalNodeCount = 0;
        isTreeRandom = true;
    }

    public NewTreeModel(Tree tree) {
        this(TREE_MODEL, tree, false, false, false);
    }

    public NewTreeModel(String id, Tree tree) { this(id, tree, false, false); }

    public NewTreeModel(String id, Tree tree, boolean fixHeights, boolean fixTree) {

        this(id, tree, false, fixHeights, fixTree);
        setId(id);
    }

    /* New constructor that copies the attributes of Tree tree into the new TreeModel
     * Useful for constructing a TreeModel from a NEXUS file entry
     */

    public NewTreeModel(String name, Tree tree, boolean copyAttributes, boolean fixHeights, boolean fixTree) {

        super(name);

        // get a rooted version of the tree to clone
        FlexibleTree binaryTree = new FlexibleTree(tree, copyAttributes);
        binaryTree.resolveTree();

        // adjust the heights to be compatible with the tip dates and perturb
        // any zero branches.
        if (!fixHeights) {
            MutableTree.Utils.correctHeightsForTips(binaryTree);
        }

        this.isTreeRandom = !fixTree;

        // clone the node structure (this will create the individual parameters)
        Node node = new Node(binaryTree, binaryTree.getRoot());

        internalNodeCount = binaryTree.getInternalNodeCount();
        externalNodeCount = binaryTree.getExternalNodeCount();

        nodeCount = internalNodeCount + externalNodeCount;

        nodes = new Node[nodeCount];
        storedNodes = new Node[nodeCount];

        int i = 0;
        int j = externalNodeCount;

        root = node;

        do {
            node = (Node) TreeUtils.postorderSuccessor(this, node);

            if (node.isExternal()) {
                node.number = i;

                nodes[i] = node;
                storedNodes[i] = new Node();
                storedNodes[i].taxon = node.taxon;
                storedNodes[i].number = i;

                i++;
            } else {
                node.number = j;

                nodes[j] = node;
                storedNodes[j] = new Node();
                storedNodes[j].number = j;

                j++;
            }
        } while (node != root);

        // must be done here to allow programmatic running of BEAST
        setupHeightBounds();
    }


    boolean heightBoundsSetup = false;

    public void setupHeightBounds() {

        if (heightBoundsSetup) {
            throw new IllegalArgumentException("Node height bounds set up twice");
        }

        for (int i = 0; i < nodeCount; i++) {
            nodes[i].setupHeightBounds();
        }

        heightBoundsSetup = true;
    }

    /**
     * Push a tree changed event into the event stack.
     */
    public void pushTreeChangedEvent() {
        pushTreeChangedEvent(new TreeChangedEvent());
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

        if (!isTreeRandom) throw new IllegalStateException("Attempting state change in fixed tree");

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
    public void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        final Node node = getNodeOfParameter((Parameter) variable);
        if (type == Parameter.ChangeType.ALL_VALUES_CHANGED) {
            //this signals events where values in all dimensions of a parameter is changed.
            pushTreeChangedEvent(new TreeChangedEvent(node, (Parameter) variable, TreeChangedEvent.CHANGE_IN_ALL_INTERNAL_NODES));
        } else {
            pushTreeChangedEvent(node, (Parameter) variable, index);
        }
    }


    private final List<TreeChangedEvent> treeChangedEvents = new ArrayList<TreeChangedEvent>();

    public boolean inTreeEdit() {
        return inEdit;
    }

    public class TreeChangedEvent implements dr.evomodel.tree.TreeChangedEvent {
        static final int CHANGE_IN_ALL_INTERNAL_NODES = -2;

        final Node node;
        final Parameter parameter;
        final int index;
        final boolean nodeOrderChanged;

        public TreeChangedEvent() {
            this(null, null, -1,false);
        }

        public TreeChangedEvent(Node node) {
            this(node, null, -1,false);
        }
        public TreeChangedEvent(Node node, Parameter parameter, int index) {
            this(node, parameter, -index,false);
        }

        public TreeChangedEvent(Node node, Parameter parameter, int index,boolean nodeOrderChanged) {
            this.node = node;
            this.parameter = parameter;
            this.index = index;
            this.nodeOrderChanged = nodeOrderChanged;
        }

        @Override
        public int getIndex() {
            return index;
        }

        @Override
        public Node getNode() {
            return node;
        }

        public Parameter getParameter() {
            return parameter;
        }

        public boolean isTreeChanged() {
            return parameter == null;
        }

        public boolean isNodeChanged() {
            return node != null;
        }
        public boolean isNodeOrderChanged(){return nodeOrderChanged;};


        public boolean isNodeParameterChanged() {
            return parameter != null;
        }

        public boolean isHeightChanged() {
            return parameter == node.heightParameter;
        }

        public boolean isOnlyHeightChanged() { return false; }

        public boolean areAllInternalHeightsChanged() {
            if (parameter != null) {
                return parameter == node.heightParameter && index == CHANGE_IN_ALL_INTERNAL_NODES;
            }
            return false;
        }

    }

    // *****************************************************************
    // Interface Tree
    // *****************************************************************

    /**
     * Return the units that this tree is expressed in.
     */
    public Type getUnits() {
        return units;
    }

    /**
     * Sets the units that this tree is expressed in.
     */
    public void setUnits(Type units) {
        this.units = units;
    }

    /**
     * @return a count of the number of nodes (internal + external) in this
     *         tree.
     */
    public int getNodeCount() {
        return nodeCount;
    }

    public boolean hasNodeHeights() {
        return true;
    }

    public double getNodeHeight(NodeRef node) {
        return ((Node) node).getHeight();
    }

    public final double getNodeHeightUpper(NodeRef node) {
        return ((Node) node).heightParameter.getBounds().getUpperLimit(0);
    }

    public final double getNodeHeightLower(NodeRef node) {
        return ((Node) node).heightParameter.getBounds().getLowerLimit(0);
    }


    public Object getNodeAttribute(NodeRef node, String name) {
        throw new UnsupportedOperationException("getNodeAttribute is not used for TreeModel");
    }

    public Iterator getNodeAttributeNames(NodeRef node) {
        throw new UnsupportedOperationException("getNodeAttribute is not used for TreeModel");
    }

    public Taxon getNodeTaxon(NodeRef node) {
        return ((Node) node).taxon;
    }

    public void setNodeTaxon(NodeRef node, Taxon taxon) {
        ((Node) node).taxon = taxon;
    }

    public boolean isExternal(NodeRef node) {
        return ((Node) node).isExternal();
    }

    public boolean isRoot(NodeRef node) {
        return (node == root);
    }

    public int getChildCount(NodeRef node) {
        return ((Node) node).getChildCount();
    }

    public NodeRef getChild(NodeRef node, int i) {
        return ((Node) node).getChild(i);
    }

    public NodeRef getParent(NodeRef node) {
        return ((Node) node).parent;
    }

    public boolean hasBranchLengths() {
        return true;
    }

    public double getBranchLength(NodeRef node) {
        NodeRef parent = getParent(node);
        if (parent == null) {
            return 0.0;
        }

        return getNodeHeight(parent) - getNodeHeight(node);
    }

    @Override
    public double getNodeRate(NodeRef node) {
        throw new UnsupportedOperationException("getNodeRate is deprecated");
    }

    @Override
    public void setNodeRate(NodeRef node, double height) {
        throw new UnsupportedOperationException("setNodeRate is deprecated");
    }


    public NodeRef getExternalNode(int i) {
        return nodes[i];
    }

    public NodeRef getInternalNode(int i) {
        return nodes[i + externalNodeCount];
    }

    public NodeRef getNode(int i) {
        return nodes[i];
    }

    public NodeRef[] getNodes() {
        return nodes;
    }

    /**
     * Returns the number of external nodes.
     */
    public int getExternalNodeCount() {
        return externalNodeCount;
    }

    /**
     * Returns the ith internal node.
     */
    public int getInternalNodeCount() {
        return internalNodeCount;
    }

    /**
     * Returns the root node of this tree.
     */
    public NodeRef getRoot() {
        return root;
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
        pushTreeChangedEvent(root);
    }

    public void addChild(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        Node parent = (Node) p;
        Node child = (Node) c;
        if (parent.hasChild(child)) throw new IllegalArgumentException("Child already exists in parent");

        parent.addChild(child);
        pushTreeChangedEvent(parent);
    }

    public void removeChild(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        Node parent = (Node) p;
        Node child = (Node) c;

        parent.removeChild(child);
    }

    public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {
        throw new RuntimeException("Unimplemented");
    }

    private Node oldRoot;

    public boolean beginTreeEdit() {
        if (inEdit) throw new RuntimeException("Alreading in edit transaction mode!");

        oldRoot = root;

        inEdit = true;

        return false;
    }

    public void endTreeEdit() {
        if (!inEdit) throw new RuntimeException("Not in edit transaction mode!");

        inEdit = false;

        if (TEST_NODE_BOUNDS) {
            try {
                checkTreeIsValid();
            } catch (InvalidTreeException ite) {
                throw new RuntimeException(ite.getMessage());
            }
        }

        for (TreeChangedEvent treeChangedEvent : treeChangedEvents) {
            listenerHelper.fireModelChanged(this, treeChangedEvent);
        }
        treeChangedEvents.clear();
    }

    public void checkTreeIsValid() throws MutableTree.InvalidTreeException {
        for (Node node : nodes) {
            if (!node.heightParameter.isWithinBounds()) {
                throw new InvalidTreeException("height parameter out of bounds");
            }
        }
    }

    public void setNodeHeight(NodeRef n, double height) {
        ((Node) n).setHeight(height);
    }

    public void setBranchLength(NodeRef node, double length) {
        throw new UnsupportedOperationException("TreeModel cannot have branch lengths set");
    }

    public void setNodeAttribute(NodeRef node, String name, Object value) {
        throw new UnsupportedOperationException("TreeModel does not use NodeAttributes");
    }

    // *****************************************************************
    // Interface ModelComponent
    // *****************************************************************

    /**
     * Store current state
     */
    protected void storeState() {

        copyNodeStructure(storedNodes);
        storedRootNumber = root.getNumber();

    }

    /**
     * Restore the stored state
     */
    protected void restoreState() {

        Node[] tmp = storedNodes;
        storedNodes = nodes;
        nodes = tmp;

        root = nodes[storedRootNumber];
    }

    /**
     * accept the stored state
     */
    protected void acceptState() {
    } // nothing to do

    /**
     * Copies the node connections from this TreeModel's nodes array to the
     * destination array. Basically it connects up the nodes in destination
     * in the same way as this TreeModel is set up. This method is package
     * private.
     */
    private void copyNodeStructure(Node[] destination) {

        if (nodes.length != destination.length) {
            throw new IllegalArgumentException("Node arrays are of different lengths");
        }

        for (int i = 0, n = nodes.length; i < n; i++) {
            Node node0 = nodes[i];
            Node node1 = destination[i];

            // the parameter values are automatically stored and restored
            // just need to keep the links
            node1.heightParameter = node0.heightParameter;

            if (node0.parent != null) {
                node1.parent = storedNodes[node0.parent.getNumber()];
            } else {
                node1.parent = null;
            }

            if (node0.leftChild != null) {
                node1.leftChild = storedNodes[node0.leftChild.getNumber()];
            } else {
                node1.leftChild = null;
            }

            if (node0.rightChild != null) {
                node1.rightChild = storedNodes[node0.rightChild.getNumber()];
            } else {
                node1.rightChild = null;
            }
        }
    }

    /**
     * Copies a different tree into the current treeModel. Needs to reconnect
     * the existing internal and external nodes, taking into account that the
     * node numbers of the external nodes may differ between the two trees.
     */
    public void adoptTreeStructure(Tree donor) {

        /*System.err.println("internalNodeCount: " + this.internalNodeCount);
          System.err.println("externalNodeCount: " + this.externalNodeCount);
          for (int i = 0; i < this.nodeCount; i++) {
              System.err.println(nodes[i]);
          }*/

        //first remove all the child nodes of the internal nodes
        for (int i = this.externalNodeCount; i < this.nodeCount; i++) {
            int childCount = nodes[i].getChildCount();
            for (int j = 0; j < childCount; j++) {
                nodes[i].removeChild(j);
            }
        }

        // set-up nodes in this.nodes[] to mirror connectedness in donor via a simple recursion on donor.getRoot()
        addNodeStructure(donor, donor.getRoot());

        //Tree donor has no rates nor traits, only heights

    }

    /**
     * Modifies the current tree by adopting the provided collection of edges
     * @param edges Edges are provided as index: child number; parent: array entry
     * @param nodeHeights Also sets the node heights to the provided values
     * @param childOrder Array that contains whether a child node is left or right child
     */
    public void adoptTreeStructure(int[] edges, double[] nodeHeights, int[] childOrder) {

        if (this.nodeCount != edges.length) {
            throw new RuntimeException("Incorrect number of edges provided: " + edges.length + " versus " + this.nodeCount + " nodes.");
        }

        //first remove all the child nodes of the internal nodes
        for (int i = this.externalNodeCount; i < this.nodeCount; i++) {
            int childCount = nodes[i].getChildCount();
            for (int j = 0; j < childCount; j++) {
                nodes[i].removeChild(j);
            }
        }

        //set the node heights
        for (int i = 0; i < nodeHeights.length; i++) {
            setNodeHeight(nodes[i], nodeHeights[i]);
        }

        int newRootIndex = -1;
        //now add the parent-child links again to ALL the nodes
        for (int i = 0; i < edges.length; i++) {
            if (edges[i] != -1) {
                nodes[edges[i]].addChild(nodes[i]);
            } else {
                newRootIndex = i;
            }
        }

        //not possible to determine correct ordering of child nodes in the loop where they're being assigned
        //hence perform possible swaps in a separate loop
        for (int i = 0; i < edges.length; i++) {
            if (edges[i] != -1) {
                if (childOrder[i] == 0 && nodes[edges[i]].getChild(0) != nodes[i]) {
                    //swap child nodes
                    Node childOne = nodes[edges[i]].removeChild(0);
                    Node childTwo = nodes[edges[i]].removeChild(1);
                    nodes[edges[i]].addChild(childTwo);
                    nodes[edges[i]].addChild(childOne);
                }
            }
        }

        this.setRoot(nodes[newRootIndex]);

    }

    /**
     * Recursive algorithm to copy a proposed tree structure into the current treeModel.
     */
    private void addNodeStructure(Tree donorTree, NodeRef donorNode) {

        NodeRef acceptorNode = null;
        if (donorTree.isExternal(donorNode)) {
            //external nodes can have different numbers between both trees
            acceptorNode = this.nodes[this.getTaxonIndex(donorTree.getTaxonId(donorNode.getNumber()))];
        } else {
            //not really important for internal nodes
            acceptorNode = this.nodes[donorNode.getNumber()];
        }

        setNodeHeight(acceptorNode, donorTree.getNodeHeight(donorNode));

        //removing all child nodes up front currently works
        //((Node)acceptorNode).leftChild = null;
        //((Node)acceptorNode).rightChild = null;
        /*int nrChildren = getChildCount(acceptorNode);
          for (int i = 0; i < nrChildren; i++) {
              this.removeChild(acceptorNode, this.getChild(acceptorNode, i));
          }*/

        for (int i = 0; i < donorTree.getChildCount(donorNode); i++) {
            //add a check when the added child is an external node
            if (donorTree.isExternal(donorTree.getChild(donorNode, i))) {
                addChild(acceptorNode, this.nodes[this.getTaxonIndex(donorTree.getTaxonId(donorTree.getChild(donorNode, i).getNumber()))]);
            } else {
                addChild(acceptorNode, this.nodes[donorTree.getChild(donorNode, i).getNumber()]);
            }
        }

        pushTreeChangedEvent(acceptorNode);

        if (!donorTree.isExternal(donorNode)) {
            for (int i = 0; i < donorTree.getChildCount(donorNode); i++) {
                addNodeStructure(donorTree, donorTree.getChild(donorNode, i));
            }
        }

    }

    /**
     * @return the number of statistics of this component.
     */
    public int getStatisticCount() {
        return super.getStatisticCount() + 1;
    }

    /**
     * @return the ith statistic of the component
     */
    public Statistic getStatistic(int i) {
        if (i == super.getStatisticCount()) return root.heightParameter;
        return super.getStatistic(i);
    }

//    public String getModelComponentName() {
//        return TREE_MODEL;
//    }

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

    public List<Taxon> asList() {
        List<Taxon> taxa = new ArrayList<Taxon>();
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            taxa.add(getTaxon(i));
        }
        return taxa;
    }

    public Iterator<Taxon> iterator() {
        return new Iterator<Taxon>() {
            private int index = -1;

            public boolean hasNext() {
                return index < getTaxonCount() - 1;
            }

            public Taxon next() {
                index++;
                return getTaxon(index);
            }

            public void remove() { /* do nothing */ }
        };
    }

    /**
     * @param taxonIndex the index of the taxon whose attribute is being fetched.
     * @param name       the name of the attribute of interest.
     * @return an object representing the named attributed for the taxon of the given
     *         external node. If the node doesn't have a taxon then the nodes own attribute
     *         is returned.
     */
    public Object getTaxonAttribute(int taxonIndex, String name) {
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
        throw new IllegalArgumentException("Cannot add taxon to a TreeModel");
    }

    public boolean removeTaxon(Taxon taxon) {
        throw new IllegalArgumentException("Cannot add taxon to a TreeModel");
    }

    public void setTaxonId(int taxonIndex, String id) {
        throw new IllegalArgumentException("Cannot set taxon id in a TreeModel");
    }

    public void setTaxonAttribute(int taxonIndex, String name, Object value) {
        throw new IllegalArgumentException("Cannot set taxon attribute in a TreeModel");
    }

    public void addMutableTreeListener(MutableTreeListener listener) {
    } // Do nothing at the moment

    public void addMutableTaxonListListener(MutableTaxonListListener listener) {
    } // Do nothing at the moment

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

    private String id = null;

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
    public Iterator<String> getAttributeNames() {
        if (treeAttributes == null)
            return null;
        else
            return treeAttributes.getAttributeNames();
    }

    /**
     * @return a string containing a newick representation of the tree
     */
    public final String getNewick() {
        return TreeUtils.newick(this);
    }

    /**
     * @return a string containing a newick representation of the tree
     */
    public String toString() {
        return getNewick();
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

    // ***********************************************************************
    // Private methods
    // ***********************************************************************

    /**
     * @return the node that this parameter is a member of
     */
    public Node getNodeOfParameter(Parameter parameter) {

        if (parameter == null) throw new IllegalArgumentException("Parameter is null!");

        for (Node node : nodes) {
            if (node.heightParameter == parameter) {
                return node;
            }
        }

        throw new RuntimeException("Parameter not found in any nodes:" + parameter.getId() + " " + parameter.hashCode());
        // assume it is a trait parameter and return null
//		return null;
    }

    /**
     * Get the root height parameter. Is private because it can only be called by the XMLParser
     */
    public Parameter getRootHeightParameter() {

        return root.heightParameter;
    }

    /**
     * @return the relevant node height parameter. Is private because it can only be called by the XMLParser
     */
    public Parameter createNodeHeightsParameter(boolean rootNode, boolean internalNodes, boolean leafNodes) {

        if (!rootNode) {
            throw new UnsupportedOperationException("A node height parameter without the root is not currently implemented");
        }

        checkValidFlags(rootNode, internalNodes, leafNodes);

        CompoundParameter parameter = new CompoundParameter("nodeHeights(" + getId() + ")");

        for (int i = externalNodeCount; i < nodeCount; i++) {
            if ((rootNode && nodes[i] == root) || (internalNodes && nodes[i] != root)) {
                parameter.addParameter(nodes[i].heightParameter);
            }
        }

        if (leafNodes) {
            for (int i = 0; i < externalNodeCount; i++) {
                parameter.addParameter(nodes[i].heightParameter);
            }
        }

        return parameter;
    }

    public Parameter getLeafHeightParameter(NodeRef node) {

        if (!isExternal(node)) {
            throw new RuntimeException("only leaves can be used with getLeafHeightParameter");
        }

        isTipDateSampled = true;

        return nodes[node.getNumber()].heightParameter;
    }

    private void checkValidFlags(boolean rootNode, boolean internalNodes, boolean leafNodes) {
        if (!rootNode && !internalNodes && !leafNodes) {
            throw new IllegalArgumentException("At least one of rootNode, internalNodes or leafNodes must be true");
        }
    }

    // **************************************************************
    // Private inner classes
    // **************************************************************

    public class Node implements NodeRef {

        public Node parent;
        public Node leftChild, rightChild;
        private int number;
        public Parameter heightParameter;
        public Taxon taxon = null;

        public Node() {
            parent = null;
            leftChild = rightChild = null;
            heightParameter = null;
            number = 0;
            taxon = null;
        }

        /**
         * constructor used to clone a node and all children
         */
        public Node(Tree tree, NodeRef node) {
            parent = null;
            leftChild = rightChild = null;

            heightParameter = new Parameter.Default(tree.getNodeHeight(node));
            addVariable(heightParameter);

            number = node.getNumber();
            taxon = tree.getNodeTaxon(node);
            heightParameter.setId(getId() + "." + number);
            for (int i = 0; i < tree.getChildCount(node); i++) {
                addChild(new Node(tree, tree.getChild(node, i)));
            }
        }

        public final void setupHeightBounds() {
            heightParameter.addBounds(new NodeHeightBounds(heightParameter));
        }

        private void setParameterValues(Parameter parameter, int dim, double[] initialValues) {
            if (initialValues != null && initialValues.length > 0) {
                for (int i = 0; i < dim; i++) {
                    if (initialValues.length == dim) {
                        parameter.setParameterValue(i, initialValues[i]);
                    } else {
                        parameter.setParameterValue(i, initialValues[0]);
                    }
                }
            }
        }

        public final double getHeight() {
            return heightParameter.getParameterValue(0);
        }

        public final void setHeight(double height) {
            heightParameter.setParameterValue(0, height);
        }

        public int getNumber() {
            return number;
        }

        public void setNumber(int n) {
            number = n;
        }

        /**
         * Returns the number of children this node has.
         */
        public final int getChildCount() {
            int n = 0;
            if (leftChild != null) n++;
            if (rightChild != null) n++;
            return n;
        }

        public Node getChild(int n) {
            if (n == 0) return leftChild;
            if (n == 1) return rightChild;
            throw new IllegalArgumentException("TreeModel.Nodes can only have 2 children");
        }

        public boolean hasChild(Node node) {
            return (leftChild == node || rightChild == node);
        }

        /**
         * add new child node
         *
         * @param node new child node
         */
        public void addChild(Node node) {
            if (leftChild == null) {
                leftChild = node;
            } else if (rightChild == null) {
                rightChild = node;
            } else {
                throw new IllegalArgumentException("TreeModel.Nodes can only have 2 children");
            }
            node.parent = this;
        }

        /**
         * remove child
         *
         * @param node child to be removed
         */
        public Node removeChild(Node node) {
            if (leftChild == node) {
                leftChild = null;
            } else if (rightChild == node) {
                rightChild = null;
            } else {
                throw new IllegalArgumentException("Unknown child node");
            }
            node.parent = null;
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
            node.parent = null;
            return node;
        }

        public boolean hasNoChildren() {
            return (leftChild == null && rightChild == null);
        }

        public boolean isExternal() {
            return hasNoChildren();
        }

        public boolean isRoot() {
            return (parent == null);
        }

        public String toString() {
            return "node " + number + ", height=" + getHeight() + (taxon != null ? ": " + taxon.getId() : "");
        }
    }

    /**
     * This class provides bounds for parameters that represent a node height
     * in this tree model.
     */
    private class NodeHeightBounds implements Bounds<Double> {

        public NodeHeightBounds(Parameter parameter) {
            nodeHeightParameter = parameter;
        }

        public Double getUpperLimit(int i) {

            Node node = getNodeOfParameter(nodeHeightParameter);
            if (node.isRoot()) {
                return Double.POSITIVE_INFINITY;
            } else {
                return node.parent.getHeight();
            }
        }

        public Double getLowerLimit(int i) {

            Node node = getNodeOfParameter(nodeHeightParameter);
            if (node.isExternal()) {
                return 0.0;
            } else {
                return Math.max(node.leftChild.getHeight(), node.rightChild.getHeight());
            }
        }

        public int getBoundsDimension() {
            return 1;
        }


        private Parameter nodeHeightParameter = null;
    }

    // ***********************************************************************
    // Interface: Keywordable
    // ***********************************************************************

    @Override
    public void addKeyword(String keyword) {
        keywords.add(keyword);
    }

    @Override
    public List<String> getKeywords() {
        return keywords;
    }

    private final List<String> keywords = new ArrayList<String>();

    // ***********************************************************************
    // Private members
    // ***********************************************************************


    /**
     * root node
     */
    private Node root = null;
    private int storedRootNumber;

    /**
     * list of internal nodes (including root)
     */
    private Node[] nodes = null;
    private Node[] storedNodes = null;

    /**
     * number of nodes (including root and tips)
     */
    private final int nodeCount;

    /**
     * number of external nodes
     */
    private final int externalNodeCount;

    /**
     * number of internal nodes (including root)
     */
    private final int internalNodeCount;

    /**
     * holds the units of the trees branches.
     */
    private Type units = Type.SUBSTITUTIONS;

    private boolean inEdit = false;

    private boolean isTipDateSampled = false;
    private final boolean isTreeRandom;

    public boolean isTipDateSampled() {
        return isTipDateSampled;
    }

    @Override
    public boolean isVariable() {
        return isTreeRandom;
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Sampling tip dates model";
    }

    @Override
    public List<Citation> getCitations() {
        if (isTipDateSampled()) {
            return Arrays.asList(new Citation(
                            new Author[]{
                                    new Author("B", "Shapiro"),
                                    new Author("SYW", "Ho"),
                                    new Author("AJ", "Drummond"),
                                    new Author("MA", "Suchard"),
                                    new Author("OG", "Pybus"),
                                    new Author("A", "Rambaut"),
                            },
                            "A Bayesian phylogenetic method to estimate unknown sequence ages",
                            2010,
                            "Mol Biol Evol",
                            28,
                            879, 887,
                            "10.1093/molbev/msq262"
                    ),
                    new Citation(
                            new Author[]{
                                    new Author("AJ", "Drummond"),
                            },
                            "PhD Thesis",
                            2002,
                            "University of Auckland",
                            ""
                    ));
        } else {
            return  Collections.EMPTY_LIST;
        }
    }

}
