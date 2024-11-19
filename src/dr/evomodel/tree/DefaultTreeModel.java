/*
 * DefaultTreeModel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.tree;

import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.inference.model.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * A model component for trees.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class DefaultTreeModel extends TreeModel {

    //
    // Public stuff
    //

    public static final String TREE_MODEL = "treeModel";

    private static final boolean TEST_NODE_BOUNDS = false;

    public DefaultTreeModel(String name) {
        super(name, true);
        nodeCount = 0;
        externalNodeCount = 0;
        internalNodeCount = 0;
    }

    public DefaultTreeModel(Tree tree) {
        this(TREE_MODEL, tree, false, false, false);
    }

    public DefaultTreeModel(String id, Tree tree) { this(id, tree, false, false); }

    public DefaultTreeModel(String id, Tree tree, boolean fixHeights, boolean fixTree) {

        this(id, tree, false, fixHeights, fixTree);
        setId(id);
    }

    /* New constructor that copies the attributes of Tree tree into the new TreeModel
     * Useful for constructing a TreeModel from a NEXUS file entry
     */

    public DefaultTreeModel(String name, Tree tree, boolean copyAttributes, boolean fixHeights, boolean fixTree) {

        super(name, !fixTree);

        // get a rooted version of the tree to clone
        FlexibleTree binaryTree = new FlexibleTree(tree, copyAttributes);
        binaryTree.resolveTree();

        // adjust the heights to be compatible with the tip dates and perturb
        // any zero branches.
        if (!fixHeights) {
            MutableTree.Utils.correctHeightsForTips(binaryTree);
        }

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
    public void pushTreeChangedEvent(Node node, Parameter parameter) {
        pushTreeChangedEvent(TreeChangedEvent.create(node, parameter, parameter == node.heightParameter));
    }

    /**
     * Push a tree changed event into the event stack.
     */
    public void pushTreeChangedEvent(Node node, Parameter parameter, int index) {
        pushTreeChangedEvent(TreeChangedEvent.create(node, parameter, index, parameter == node.heightParameter));
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
            pushTreeChangedEvent(node, (Parameter) variable);
        } else {
            pushTreeChangedEvent(node, (Parameter) variable, index);
        }
    }

    @Override
    public boolean beginTreeEdit() {
        oldRoot = root;
        return super.beginTreeEdit();
    }

    @Override
    public void endTreeEdit() {
        super.endTreeEdit();
        if (root != oldRoot) {
            swapParameterObjects(oldRoot, root);
        }

    }

    public boolean hasRates() {
        return hasRates;
    }

    // *****************************************************************
    // Interface Tree
    // *****************************************************************

    /**
     * @return a count of the number of nodes (internal + external) in this
     *         tree.
     */
    public int getNodeCount() {
        return nodeCount;
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


    /**
     * @param node
     * @return the rate parameter associated with this node.
     */
    public double getNodeRate(NodeRef node) {
        if (!hasRates) {
            return 1.0;
        }
        return ((Node) node).getRate();
    }

    public Object getNodeAttribute(NodeRef node, String name) {

        if (name.equals("rate")) {
            return getNodeRate(node);
        }

        return null;
    }

    public Iterator getNodeAttributeNames(NodeRef node) {
        return new Iterator() {

            int i = 0;
            String[] attributes = {"rate"};

            public boolean hasNext() {
                return i < attributes.length;
            }

            public Object next() {
                return attributes[i++];
            }

            public void remove() {
                throw new UnsupportedOperationException("can't remove from this iterator!");
            }
        };
    }

    public boolean hasNodeTraits() {
        return hasTraits;
    }

    public Map<String, Parameter> getTraitMap(NodeRef node) {
        if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
        return ((Node) node).getTraitMap();
    }

    public double getNodeTrait(NodeRef node, String name) {
        if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
        return ((Node) node).getTrait(name);
    }

    public Parameter getNodeTraitParameter(NodeRef node, String name) {
        if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
        return ((Node) node).getTraitParameter(name);
    }

    public double[] getMultivariateNodeTrait(NodeRef node, String name) {
        if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
        return ((Node) node).getMultivariateTrait(name);
    }

    public final void swapAllTraits(NodeRef node1, NodeRef node2) {
        if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
        swapAllTraits((Node) node1, (Node) node2);
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

    public void setTipDateSampled(boolean tipDateSampled) {
        isTipDateSampled = tipDateSampled;
    }

    @Override
    public boolean isTipDateSampled() {
        return isTipDateSampled;
    }


    // *****************************************************************
    // Interface MutableTree
    // *****************************************************************

    /**
     * Set a new node as root node.
     */
    public void setRoot(NodeRef newRoot) {

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

    @Override
    public boolean isTreeValid() {
        for (Node node : nodes) {
            if (!node.heightParameter.isWithinBounds()) {
                return false;
            }
        }
        return true;
    }

    public void setNodeHeight(NodeRef n, double height) {
        ((Node) n).setHeight(height);
    }

    @Override
    public void setNodeHeightQuietly(NodeRef n, double height) {
        ((Node) n).setHeightQuietly(height);
    }

    public void setNodeRate(NodeRef n, double rate) {
        if (!hasRates) throw new IllegalArgumentException("Rate parameters have not been created");
        ((Node) n).setRate(rate);

    }

    public void setNodeTrait(NodeRef n, String name, double value) {
        if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
        ((Node) n).setTrait(name, value);
    }

    public void setMultivariateTrait(NodeRef n, String name, double[] value) {
        if (!hasTraits) throw new IllegalArgumentException("Trait parameters have not been created");
        ((Node) n).setMultivariateTrait(name, value);
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

        remapParameterNodes();

    }

    private void remapParameterNodes() {
        for (Node node : nodes) {
            parameterNodeMap.put(node.heightParameter, node);
            if (hasRates) {
                parameterNodeMap.put(node.rateParameter, node);
            }
            if (hasTraits) {
                for (Parameter trait : node.traitParameters.values()) {
                    parameterNodeMap.put(trait, node);
                }
            }
        }
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
            node1.rateParameter = node0.rateParameter;
            node1.traitParameters = node0.traitParameters;

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

        //first remove all the child nodes of the internal nodes
        for (int i = this.externalNodeCount; i < this.nodeCount; i++) {
            int childCount = nodes[i].getChildCount();
            for (int j = 0; j < childCount; j++) {
                nodes[i].removeChild(j);
            }
        }

        // set-up nodes in this.nodes[] to mirror connectedness in donor via a simple recursion on donor.getRoot()
        addNodeStructure(donor, donor.getRoot());

        remapParameterNodes();
    }

    /**
     * Modifies the current tree by adopting the provided collection of edges
     * @param edges Edges are provided as index: child number; parent: array entry
     * @param nodeHeights Also sets the node heights to the provided values
     * @param childOrder Array that contains whether a child node is left or right child
     */
    public void adoptTreeStructure(int[] edges, double[] nodeHeights, int[] childOrder, String[] taxaNames) {

        int[] nodeMap = createNodeMap(taxaNames);

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

        //start with setting the external node heights
        for (int i = 0; i < this.getExternalNodeCount(); i++) {
            this.setNodeHeight(this.getExternalNode(nodeMap[i]), nodeHeights[i]);
        }
        //set the internal node heights
        for (int i = 0; i < (this.getExternalNodeCount() - 1); i++) {
            //No just restart counting, will fix later on in the code by adding additionalTaxa variable
            this.setNodeHeight(this.getInternalNode(i), nodeHeights[this.getExternalNodeCount() + i]);
        }

        int newRootIndex = -1;
        //now add the parent-child links again to ALL the nodes
        for (int i = 0; i < edges.length; i++) {
            if (edges[i] != -1) {
                //make distinction between external nodes and internal nodes
                if (i < this.getExternalNodeCount()) {
                    //external node
                    this.addChild(this.getNode(edges[i]), this.getExternalNode(nodeMap[i]));
                    System.out.println("external: " + edges[i] + " > " + nodeMap[i]);
                } else {
                    //internal node
                    this.addChild(this.getNode(edges[i]), this.getNode(i));
                    System.out.println("internal: " + edges[i] + " > " + i);
                }
            } else {
                newRootIndex = i;
            }
        }

        //not possible to determine correct ordering of child nodes in the loop where they're being assigned
        //hence perform possible swaps in a separate loop

        for (int i = 0; i < edges.length; i++) {
            if (edges[i] != -1) {
                if(i < this.externalNodeCount) {
                    if (childOrder[i] == 0 && nodes[edges[i]].getChild(0) != nodes[nodeMap[i]]) {
                        //swap child nodes
                        Node childOne = nodes[edges[i]].removeChild(0);
                        Node childTwo = nodes[edges[i]].removeChild(1);
                        nodes[edges[i]].addChild(childTwo);
                        nodes[edges[i]].addChild(childOne);
                    }
                }else{
                    if (childOrder[i] == 0 && nodes[edges[i]].getChild(0) != nodes[i]) {
                        //swap child nodes
                        Node childOne = nodes[edges[i]].removeChild(0);
                        Node childTwo = nodes[edges[i]].removeChild(1);
                        nodes[edges[i]].addChild(childTwo);
                        nodes[edges[i]].addChild(childOne);
                    }
                }
            }

        }

        remapParameterNodes();
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

    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

    /**
     * @return the ith taxon in the list.
     */
    public Taxon getTaxon(int taxonIndex) {
        return ((Node) getExternalNode(taxonIndex)).taxon;
    }

    // **************************************************************
    // MutableTaxonList IMPLEMENTATION
    // **************************************************************

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

    // **************************************************************
    // Attributable IMPLEMENTATION
    // **************************************************************

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    // ***********************************************************************
    // Private methods
    // ***********************************************************************

    /**
     * @return the node that this parameter is a member of
     */
    public Node getNodeOfParameter(Parameter parameter) {
        return parameterNodeMap.get(parameter);
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

        if (!rootNode && !internalNodes && !leafNodes) {
            throw new IllegalArgumentException("At least one of rootNode, internalNodes or leafNodes must be true");
        }

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

        setTipDateSampled(true);

        return nodes[node.getNumber()].heightParameter;
    }

    /**
     * @return the relevant node rate parameter. Is private because it can only be called by the XMLParser
     */
    public Parameter createNodeRatesParameter(double[] initialValues, boolean rootNode, boolean internalNodes, boolean leafNodes) {

        if (!rootNode && !internalNodes && !leafNodes) {
            throw new IllegalArgumentException("At least one of rootNode, internalNodes or leafNodes must be true");
        }

        CompoundParameter parameter = new CompoundParameter("nodeRates(" + getId() + ")");

        hasRates = true;

        for (int i = externalNodeCount; i < nodeCount; i++) {
            nodes[i].createRateParameter(initialValues);
            if ((rootNode && nodes[i] == root) || (internalNodes && nodes[i] != root)) {
                parameter.addParameter(nodes[i].rateParameter);
            }
        }

        for (int i = 0; i < externalNodeCount; i++) {
            nodes[i].createRateParameter(initialValues);
            if (leafNodes) {
                parameter.addParameter(nodes[i].rateParameter);
            }
        }

        return parameter;
    }

    public Parameter createNodeTraitsParameter(String name, double[] initialValues) {
        return createNodeTraitsParameter(name, initialValues.length,
                initialValues, true, true, true, true);
    }

    /**
     * Create a node traits parameter. Is private because it can only be called by the XMLParser
     */
    public Parameter createNodeTraitsParameter(String name, int dim, double[] initialValues,
                                               boolean rootNode, boolean internalNodes,
                                               boolean leafNodes, boolean firesTreeEvents) {

        checkValidFlags(rootNode, internalNodes, leafNodes);

        CompoundParameter parameter = new CompoundParameter(name);

        hasTraits = true;

        for (int i = externalNodeCount; i < nodeCount; i++) {
            nodes[i].createTraitParameter(name, dim, initialValues, firesTreeEvents);
            if ((rootNode && nodes[i] == root) || (internalNodes && nodes[i] != root)) {
                parameter.addParameter(nodes[i].getTraitParameter(name));
            }
        }

        for (int i = 0; i < externalNodeCount; i++) {
            nodes[i].createTraitParameter(name, dim, initialValues, firesTreeEvents);
            if (leafNodes) {
                parameter.addParameter(nodes[i].getTraitParameter(name));
            }
        }

        return parameter;
    }

    public Parameter createNodeTraitsParameterAsMatrix(String name, int dim, double[] initialValues,
                                                       boolean rootNode, boolean internalNodes,
                                                       boolean leafNodes, boolean firesTreeEvents,
                                                       boolean signalComponents) {

        checkValidFlags(rootNode, internalNodes, leafNodes);

        final int rowDim = dim;
        final int colDim = (rootNode ? 1 : 0)
                + (internalNodes ? internalNodeCount - 1 : 0)
                + (leafNodes ? externalNodeCount : 0);

        FastMatrixParameter parameter = new FastMatrixParameter(name, rowDim, colDim, 0.0, signalComponents);
        parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
                rowDim * colDim));

        hasTraits = true;

        int parameterIndex = 0;
        for (int i = externalNodeCount; i < nodeCount; i++) {
            if ((rootNode && nodes[i] == root) || (internalNodes && nodes[i] != root)) {
                nodes[i].addTraitParameter(name, parameter.getParameter(parameterIndex), initialValues, firesTreeEvents);
                ++parameterIndex;
            }
        }

        for (int i = 0; i < externalNodeCount; i++) {
            if (leafNodes) {
                nodes[i].addTraitParameter(name, parameter.getParameter(parameterIndex), initialValues, firesTreeEvents);
                ++parameterIndex;
            }
        }

        return parameter;
    }

    private void checkValidFlags(boolean rootNode, boolean internalNodes, boolean leafNodes) {
        if (!rootNode && !internalNodes && !leafNodes) {
            throw new IllegalArgumentException("At least one of rootNode, internalNodes or leafNodes must be true");
        }
    }

    private void swapAllTraits(Node n1, Node n2) {

        for (Map.Entry<String, Parameter> entry : n1.traitParameters.entrySet()) {
            Parameter p1 = n1.traitParameters.get(entry.getKey());
            Parameter p2 = n2.traitParameters.get(entry.getKey());
            final int dim = p1.getDimension();
            for (int i = 0; i < dim; i++) {
                double transfer = p1.getParameterValue(i);
                p1.setParameterValue(i, p2.getParameterValue(i));
                p2.setParameterValue(i, transfer);
            }
            parameterNodeMap.put(p1, n1);
            parameterNodeMap.put(p2, n2);
        }

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

        if (hasRates) {
            rate1 = n1.getRate();
            rate2 = n2.getRate();
        }

        // swap all trait parameters

        if (hasTraits) {
            Map<String, Parameter> traits1 = new HashMap<String, Parameter>();
            Map<String, Parameter> traits2 = new HashMap<String, Parameter>();

            traits1.putAll(n1.traitParameters);
            traits2.putAll(n2.traitParameters);

            Map<String, Parameter> temp = n1.traitParameters;
            n1.traitParameters = n2.traitParameters;
            n2.traitParameters = temp;

            for (Parameter trait : n1.traitParameters.values()) {
                parameterNodeMap.put(trait, n1);
            }
            for (Parameter trait : n2.traitParameters.values()) {
                parameterNodeMap.put(trait, n2);
            }

            for (Map.Entry<String, Parameter> entry : traits1.entrySet()) {
                n1.traitParameters.get(entry.getKey()).setParameterValueQuietly(0, entry.getValue().getParameterValue(0));
            }
            for (Map.Entry<String, Parameter> entry : traits2.entrySet()) {
                n2.traitParameters.get(entry.getKey()).setParameterValueQuietly(0, entry.getValue().getParameterValue(0));
            }
        }

        Parameter temp = n1.heightParameter;
        n1.heightParameter = n2.heightParameter;
        n2.heightParameter = temp;
        parameterNodeMap.put(n1.heightParameter, n1);
        parameterNodeMap.put(n2.heightParameter, n2);

        if (hasRates) {
            temp = n1.rateParameter;
            n1.rateParameter = n2.rateParameter;
            n2.rateParameter = temp;
            parameterNodeMap.put(n1.rateParameter, n1);
            parameterNodeMap.put(n2.rateParameter, n2);
        }

        n1.heightParameter.setParameterValueQuietly(0, height1);
        n2.heightParameter.setParameterValueQuietly(0, height2);

        if (hasRates) {
            n1.rateParameter.setParameterValueQuietly(0, rate1);
            n2.rateParameter.setParameterValueQuietly(0, rate2);
        }
    }

    // **************************************************************
    // Private inner classes
    // **************************************************************

    protected class Node implements NodeRef {

        public Node parent;
        public Node leftChild, rightChild;
        private int number;
        public Parameter heightParameter;
        public Parameter rateParameter = null;
        //public Parameter traitParameter = null;
        public Taxon taxon = null;

        Map<String, Parameter> traitParameters = new HashMap<String, Parameter>();

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
            parameterNodeMap.put(heightParameter, this);

            number = node.getNumber();
            taxon = tree.getNodeTaxon(node);
            heightParameter.setId("" + number);
            for (int i = 0; i < tree.getChildCount(node); i++) {
                addChild(new Node(tree, tree.getChild(node, i)));
            }
        }

        public final void setupHeightBounds() {
            heightParameter.addBounds(new NodeHeightBounds(heightParameter));
        }

        public final void createRateParameter(double[] initialValues) {
            if (rateParameter == null) {
                if (initialValues != null) {
                    rateParameter = new Parameter.Default(initialValues[0]);
                } else {
                    rateParameter = new Parameter.Default(1.0);
                }
                setParameterId("rate", rateParameter);
                rateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
                parameterNodeMap.put(rateParameter, this);
                addVariable(rateParameter);
            }
        }

        public final void addTraitParameter(String name, Parameter trait, double[] initialValues, boolean firesTreeEvents) {
            if (!traitParameters.containsKey(name)) {
                setParameterId(name, trait);
                setParameterValues(trait, trait.getDimension(), initialValues);

                traitParameters.put(name, trait);

                if (firesTreeEvents) {
                    addVariable(trait);
                }
            }
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

        public final void createTraitParameter(String name, int dim, double[] initialValues, boolean firesTreeEvents) {

            if (!traitParameters.containsKey(name)) {

                Parameter trait = new Parameter.Default(dim);
                setParameterId(name, trait);
                trait.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, dim));

                setParameterValues(trait, dim, initialValues);

                traitParameters.put(name, trait);
                parameterNodeMap.put(trait, this);

                if (firesTreeEvents) {
                    addVariable(trait);
                }
            }
        }

        private void setParameterId(String name, Parameter trait) {
            if (isRoot()) {
                trait.setId("root." + name);
            } else if (isExternal()) {
                trait.setId(getTaxonId(getNumber()) + "." + name);
            } else {
                trait.setId("node" + getNumber() + "." + name);
            }
        }

        public final Parameter getHeightParameter() { return heightParameter; }

        public final double getHeight() {
            return heightParameter.getParameterValue(0);
        }

        public final double getRate() {
            return rateParameter.getParameterValue(0);
        }

        public final double getTrait(String name) {
            return traitParameters.get(name).getParameterValue(0);
        }

        public final double[] getMultivariateTrait(String name) {
            return traitParameters.get(name).getParameterValues();
        }

        public final Map<String, Parameter> getTraitMap() {
            return traitParameters;
        }

        public final void setHeight(double height) {
            heightParameter.setParameterValue(0, height);
        }

        public final void setHeightQuietly(double height) {
            heightParameter.setParameterValueQuietly(0, height);
        }

        public final void setRate(double rate) {
            //System.out.println("Rate set for parameter " + rateParameter.getParameterName());
            rateParameter.setParameterValue(0, rate);
        }

        public final void setTrait(String name, double trait) {
            //System.out.println("Trait set for parameter " + traitParameter.getParameterName());
            traitParameters.get(name).setParameterValue(0, trait);
        }

        public final void setMultivariateTrait(String name, double[] trait) {
            int dim = trait.length;
            for (int i = 0; i < dim; i++)
                traitParameters.get(name).setParameterValue(i, trait[i]);
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

        public Parameter getTraitParameter(String name) {
            return traitParameters.get(name);
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

    private boolean hasRates = false;
    private boolean hasTraits = false;
    private boolean isTipDateSampled = false;

    private final Map<Parameter, Node> parameterNodeMap = new HashMap<>();

}
