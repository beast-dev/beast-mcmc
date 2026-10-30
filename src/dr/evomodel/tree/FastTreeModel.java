/*
 * FastTreeModel.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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

import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;

/**
 * A tree model for large trees that behaves exactly as DefaultTreeModel (the same node numbering,
 * height parameters, bounds and events, so the same MCMC output), but keeps the tree in flat arrays
 * and stores it copy-on-write, so storing, restoring and accepting cost in proportion to the number
 * of nodes that changed. The node height parameters are views on the arrays: they hold no values and
 * are not variables of the model, so they are not visited when the model is stored.
 * <p>
 * Node rates and node traits are not supported (use DefaultTreeModel).
 *
 * @author Andrew Rambaut
 */
public class FastTreeModel extends TreeModel {

    public FastTreeModel(String id, Tree tree, boolean fixHeights, boolean fixTree) {
        this(id, tree, false, fixHeights, fixTree);
        setId(id);
    }

    public FastTreeModel(String name, Tree tree, boolean copyAttributes, boolean fixHeights, boolean fixTree) {

        super(name, !fixTree);

        // get a rooted version of the tree to clone (as DefaultTreeModel)
        FlexibleTree binaryTree = new FlexibleTree(tree, copyAttributes);
        binaryTree.resolveTree();

        // adjust the heights to be compatible with the tip dates and perturb
        // any zero branches.
        if (!fixHeights) {
            MutableTree.Utils.correctHeightsForTips(binaryTree);
        }

        internalNodeCount = binaryTree.getInternalNodeCount();
        externalNodeCount = binaryTree.getExternalNodeCount();
        nodeCount = internalNodeCount + externalNodeCount;

        nodes = new Node[nodeCount];
        heights = new double[nodeCount];
        parent = new int[nodeCount];
        child0 = new int[nodeCount];
        child1 = new int[nodeCount];
        nodeSlot = new int[nodeCount];
        slotNode = new int[nodeCount];
        slotParameters = new SlotHeightParameter[nodeCount];

        // number the nodes as DefaultTreeModel does: tips then internal nodes, each in post-order
        int[] number = new int[binaryTree.getNodeCount()];
        int i = 0;
        int j = externalNodeCount;
        NodeRef binaryRoot = binaryTree.getRoot();
        NodeRef node = binaryRoot;
        do {
            node = TreeUtils.postorderSuccessor(binaryTree, node);
            if (binaryTree.isExternal(node)) {
                number[node.getNumber()] = i;
                i++;
            } else {
                number[node.getNumber()] = j;
                j++;
            }
        } while (node != binaryRoot);

        for (int k = 0; k < binaryTree.getNodeCount(); k++) {
            NodeRef binaryNode = binaryTree.getNode(k);
            int n = number[binaryNode.getNumber()];

            nodes[n] = new Node(n, binaryTree.getNodeTaxon(binaryNode));
            heights[n] = binaryTree.getNodeHeight(binaryNode);
            parent[n] = binaryTree.isRoot(binaryNode) ? -1 : number[binaryTree.getParent(binaryNode).getNumber()];
            int childCount = binaryTree.getChildCount(binaryNode);
            child0[n] = childCount > 0 ? number[binaryTree.getChild(binaryNode, 0).getNumber()] : -1;
            child1[n] = childCount > 1 ? number[binaryTree.getChild(binaryNode, 1).getNumber()] : -1;

            // each node starts with its own height parameter (slot), named as in DefaultTreeModel
            nodeSlot[n] = n;
            slotNode[n] = n;
            slotParameters[n] = new SlotHeightParameter(n);
            slotParameters[n].setId("" + binaryNode.getNumber());
            // the tree is the first listener of each height parameter, as in DefaultTreeModel (where
            // they are added with addVariable), but they are not variables of this model
            slotParameters[n].addParameterListener(this);
        }
        root = number[binaryRoot.getNumber()];

        savedHeight = new double[nodeCount];
        savedParent = new int[nodeCount];
        savedChild0 = new int[nodeCount];
        savedChild1 = new int[nodeCount];
        savedSlot = new int[nodeCount];
        isChanged = new boolean[nodeCount];
        changedNodes = new int[nodeCount];

        // must be done here to allow programmatic running of BEAST
        setupHeightBounds();
    }

    private void setupHeightBounds() {
        for (int n = 0; n < nodeCount; n++) {
            slotParameters[n].addBounds(new NodeHeightBounds(n));
        }
    }

    // *****************************************************************
    // Height parameters
    // *****************************************************************

    public Parameter getRootHeightParameter() {
        return slotParameters[nodeSlot[root]];
    }

    /**
     * @return a parameter of the heights of the selected nodes, in the same order as
     * DefaultTreeModel.createNodeHeightsParameter
     */
    public Parameter createNodeHeightsParameter(boolean rootNode, boolean internalNodes, boolean leafNodes) {

        if (!rootNode && !internalNodes && !leafNodes) {
            throw new IllegalArgumentException("At least one of rootNode, internalNodes or leafNodes must be true");
        }

        NodeHeightsParameter parameter = new NodeHeightsParameter("nodeHeights(" + getId() + ")");

        for (int n = externalNodeCount; n < nodeCount; n++) {
            if ((rootNode && n == root) || (internalNodes && n != root)) {
                parameter.addParameter(slotParameters[nodeSlot[n]]);
            }
        }

        if (leafNodes) {
            for (int n = 0; n < externalNodeCount; n++) {
                parameter.addParameter(slotParameters[nodeSlot[n]]);
            }
        }

        return parameter;
    }

    /**
     * @return a parameter of the given leaves' heights (as built by TreeModelParser for leafHeights)
     */
    public CompoundParameter createLeafHeightsParameter(String name) {
        return new NodeHeightsParameter(name);
    }

    public Parameter getLeafHeightParameter(NodeRef node) {

        if (!isExternal(node)) {
            throw new RuntimeException("only leaves can be used with getLeafHeightParameter");
        }

        setTipDateSampled(true);

        return slotParameters[nodeSlot[node.getNumber()]];
    }

    /**
     * Called when a node height parameter changes (this tree is its first listener).
     */
    public void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        final SlotHeightParameter parameter = (SlotHeightParameter) variable;
        final Node node = nodes[slotNode[parameter.slot]];
        if (type == Parameter.ChangeType.ALL_VALUES_CHANGED) {
            //this signals events where values in all dimensions of a parameter is changed.
            pushTreeChangedEvent(TreeChangedEvent.create(node, parameter, true));
        } else {
            pushTreeChangedEvent(TreeChangedEvent.create(node, parameter, index, true));
        }
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // no submodels so nothing to do
    }

    // *****************************************************************
    // Edits
    // *****************************************************************

    @Override
    public boolean beginTreeEdit() {
        oldRoot = root;
        return super.beginTreeEdit();
    }

    @Override
    public void endTreeEdit() {
        super.endTreeEdit();
        if (root != oldRoot) {
            // keep the root height parameter at the root (as DefaultTreeModel swaps parameter objects)
            swapSlots(oldRoot, root);
        }
    }

    private void swapSlots(int n1, int n2) {
        saveNode(n1);
        saveNode(n2);
        int slot1 = nodeSlot[n1];
        int slot2 = nodeSlot[n2];
        nodeSlot[n1] = slot2;
        nodeSlot[n2] = slot1;
        slotNode[slot2] = n1;
        slotNode[slot1] = n2;
    }

    // *****************************************************************
    // Interface Tree
    // *****************************************************************

    public int getNodeCount() {
        return nodeCount;
    }

    public double getNodeHeight(NodeRef node) {
        return heights[node.getNumber()];
    }

    public final double getNodeHeightUpper(NodeRef node) {
        return slotParameters[nodeSlot[node.getNumber()]].getBounds().getUpperLimit(0);
    }

    public final double getNodeHeightLower(NodeRef node) {
        return slotParameters[nodeSlot[node.getNumber()]].getBounds().getLowerLimit(0);
    }

    /**
     * @return 1.0: there are no node rates (as DefaultTreeModel without node rates; generic tree
     * copying asks for them)
     */
    public double getNodeRate(NodeRef node) {
        return 1.0;
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
            final String[] attributes = {"rate"};

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

    public Taxon getNodeTaxon(NodeRef node) {
        return ((Node) node).taxon;
    }

    public boolean isExternal(NodeRef node) {
        final int n = node.getNumber();
        return child0[n] == -1 && child1[n] == -1;
    }

    public boolean isRoot(NodeRef node) {
        return node.getNumber() == root;
    }

    public int getChildCount(NodeRef node) {
        final int n = node.getNumber();
        int count = 0;
        if (child0[n] != -1) count++;
        if (child1[n] != -1) count++;
        return count;
    }

    public NodeRef getChild(NodeRef node, int i) {
        final int c;
        if (i == 0) {
            c = child0[node.getNumber()];
        } else if (i == 1) {
            c = child1[node.getNumber()];
        } else {
            throw new IllegalArgumentException("TreeModel.Nodes can only have 2 children");
        }
        return c == -1 ? null : nodes[c];
    }

    public NodeRef getParent(NodeRef node) {
        final int p = parent[node.getNumber()];
        return p == -1 ? null : nodes[p];
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

    public int getExternalNodeCount() {
        return externalNodeCount;
    }

    public int getInternalNodeCount() {
        return internalNodeCount;
    }

    public NodeRef getRoot() {
        return nodes[root];
    }

    public void setTipDateSampled(boolean tipDateSampled) {
        isTipDateSampled = tipDateSampled;
    }

    @Override
    public boolean isTipDateSampled() {
        return isTipDateSampled;
    }

    public Taxon getTaxon(int taxonIndex) {
        return nodes[taxonIndex].taxon;
    }

    // *****************************************************************
    // Interface MutableTree
    // *****************************************************************

    public void setRoot(NodeRef newRoot) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        root = newRoot.getNumber();

        // We shouldn't need this because the addChild will already have fired appropriate events.
        pushTreeChangedEvent(nodes[root]);
    }

    public void addChild(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        final int n = p.getNumber();
        final int child = c.getNumber();
        if (child0[n] == child || child1[n] == child) throw new IllegalArgumentException("Child already exists in parent");

        addChild(n, child);
        pushTreeChangedEvent(nodes[n]);
    }

    public void removeChild(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        removeChild(p.getNumber(), c.getNumber());
    }

    private void addChild(int n, int child) {
        saveNode(n);
        saveNode(child);
        if (child0[n] == -1) {
            child0[n] = child;
        } else if (child1[n] == -1) {
            child1[n] = child;
        } else {
            throw new IllegalArgumentException("TreeModel.Nodes can only have 2 children");
        }
        parent[child] = n;
    }

    private void removeChild(int n, int child) {
        saveNode(n);
        saveNode(child);
        if (child0[n] == child) {
            child0[n] = -1;
        } else if (child1[n] == child) {
            child1[n] = -1;
        } else {
            throw new IllegalArgumentException("Unknown child node");
        }
        parent[child] = -1;
    }

    /**
     * Removes the i'th child of node n (with no checks or events).
     */
    private void removeChildAt(int n, int i) {
        final int child = (i == 0 ? child0[n] : child1[n]);
        saveNode(n);
        saveNode(child);
        if (i == 0) {
            child0[n] = -1;
        } else {
            child1[n] = -1;
        }
        parent[child] = -1;
    }

    public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {
        throw new RuntimeException("Unimplemented");
    }

    @Override
    public boolean isTreeValid() {
        for (int n = 0; n < nodeCount; n++) {
            if (!slotParameters[nodeSlot[n]].isWithinBounds()) {
                return false;
            }
        }
        return true;
    }

    public void setNodeHeight(NodeRef n, double height) {
        slotParameters[nodeSlot[n.getNumber()]].setParameterValue(0, height);
    }

    @Override
    public void setNodeHeightQuietly(NodeRef n, double height) {
        slotParameters[nodeSlot[n.getNumber()]].setParameterValueQuietly(0, height);
    }

    public void setBranchLength(NodeRef node, double length) {
        throw new UnsupportedOperationException("TreeModel cannot have branch lengths set");
    }

    public void setNodeAttribute(NodeRef node, String name, Object value) {
        throw new UnsupportedOperationException("TreeModel does not use NodeAttributes");
    }

    // *****************************************************************
    // Node rates and traits: not supported
    // *****************************************************************

    public void setNodeRate(NodeRef n, double rate) {
        throw new UnsupportedOperationException(UNSUPPORTED);
    }

    public Parameter createNodeRatesParameter(double[] initialValues, boolean rootNode, boolean internalNodes, boolean leafNodes) {
        throw new UnsupportedOperationException(UNSUPPORTED);
    }

    public Parameter createNodeTraitsParameter(String name, int dim, double[] initialValues,
                                               boolean rootNode, boolean internalNodes,
                                               boolean leafNodes, boolean firesTreeEvents) {
        throw new UnsupportedOperationException(UNSUPPORTED);
    }

    public boolean hasRates() {
        return false;
    }

    public boolean hasNodeTraits() {
        return false;
    }

    public Map<String, Parameter> getTraitMap(NodeRef node) {
        throw new UnsupportedOperationException(UNSUPPORTED);
    }

    public double getNodeTrait(NodeRef node, String name) {
        throw new UnsupportedOperationException(UNSUPPORTED);
    }

    public Parameter getNodeTraitParameter(NodeRef node, String name) {
        throw new UnsupportedOperationException(UNSUPPORTED);
    }

    public double[] getMultivariateNodeTrait(NodeRef node, String name) {
        throw new UnsupportedOperationException(UNSUPPORTED);
    }

    public void setNodeTrait(NodeRef n, String name, double value) {
        throw new UnsupportedOperationException(UNSUPPORTED);
    }

    public void setMultivariateTrait(NodeRef n, String name, double[] value) {
        throw new UnsupportedOperationException(UNSUPPORTED);
    }

    public void setNodeTaxon(NodeRef node, Taxon taxon) {
        throw new UnsupportedOperationException("FastTreeModel does not allow taxa to be changed");
    }

    private static final String UNSUPPORTED = "FastTreeModel does not support node rates or node traits: use a treeModel without implementation=\"fast\"";

    // *****************************************************************
    // Interface ModelComponent: copy-on-write store
    // *****************************************************************

    /**
     * Saves a node's height, links and height parameter slot the first time it changes after a store.
     */
    private void saveNode(int n) {
        if (!isChanged[n]) {
            savedHeight[n] = heights[n];
            savedParent[n] = parent[n];
            savedChild0[n] = child0[n];
            savedChild1[n] = child1[n];
            savedSlot[n] = nodeSlot[n];
            isChanged[n] = true;
            changedNodes[changedNodeCount] = n;
            changedNodeCount++;
        }
    }

    protected void storeState() {
        for (int i = 0; i < changedNodeCount; i++) {
            isChanged[changedNodes[i]] = false;
        }
        changedNodeCount = 0;
        storedRoot = root;
    }

    protected void restoreState() {
        for (int i = 0; i < changedNodeCount; i++) {
            final int n = changedNodes[i];
            heights[n] = savedHeight[n];
            parent[n] = savedParent[n];
            child0[n] = savedChild0[n];
            child1[n] = savedChild1[n];
            nodeSlot[n] = savedSlot[n];
        }
        for (int i = 0; i < changedNodeCount; i++) {
            final int n = changedNodes[i];
            slotNode[nodeSlot[n]] = n;
            isChanged[n] = false;
        }
        changedNodeCount = 0;
        root = storedRoot;
    }

    protected void acceptState() {
    } // nothing to do

    // *****************************************************************
    // Adopting a tree (e.g. from a checkpoint)
    // *****************************************************************

    /**
     * Copies a different tree into this tree model (as DefaultTreeModel.adoptTreeStructure(Tree)).
     */
    public void adoptTreeStructure(Tree donor) {

        //first remove all the child nodes of the internal nodes
        for (int n = this.externalNodeCount; n < this.nodeCount; n++) {
            int childCount = getChildCount(nodes[n]);
            for (int j = 0; j < childCount; j++) {
                removeChildAt(n, j);
            }
        }

        // set-up nodes to mirror connectedness in donor via a simple recursion on donor.getRoot()
        addNodeStructure(donor, donor.getRoot());
    }

    private void addNodeStructure(Tree donorTree, NodeRef donorNode) {

        NodeRef acceptorNode;
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
     * Modifies the current tree by adopting the provided collection of edges (as
     * DefaultTreeModel.adoptTreeStructure).
     */
    public void adoptTreeStructure(int[] edges, double[] nodeHeights, int[] childOrder, String[] taxaNames) {

        int[] nodeMap = createNodeMap(taxaNames);

        if (this.nodeCount != edges.length) {
            throw new RuntimeException("Incorrect number of edges provided: " + edges.length + " versus " + this.nodeCount + " nodes.");
        }

        //first remove all the child nodes of the internal nodes
        for (int n = this.externalNodeCount; n < this.nodeCount; n++) {
            int childCount = getChildCount(nodes[n]);
            for (int j = 0; j < childCount; j++) {
                removeChildAt(n, j);
            }
        }

        //start with setting the external node heights
        for (int i = 0; i < this.getExternalNodeCount(); i++) {
            this.setNodeHeight(this.getExternalNode(nodeMap[i]), nodeHeights[i]);
        }
        //set the internal node heights
        for (int i = 0; i < (this.getExternalNodeCount() - 1); i++) {
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
                int child = i < this.externalNodeCount ? nodeMap[i] : i;
                if (childOrder[i] == 0 && child0[edges[i]] != child) {
                    //swap child nodes
                    int childOne = child0[edges[i]];
                    int childTwo = child1[edges[i]];
                    removeChildAt(edges[i], 0);
                    removeChildAt(edges[i], 1);
                    addChild(edges[i], childTwo);
                    addChild(edges[i], childOne);
                }
            }
        }

        this.setRoot(nodes[newRootIndex]);
    }

    // *****************************************************************
    // Statistics
    // *****************************************************************

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
        if (i == super.getStatisticCount()) return getRootHeightParameter();
        return super.getStatistic(i);
    }

    // *****************************************************************
    // Inner classes
    // *****************************************************************

    /**
     * A node: just a number (and a taxon for tips); the tree is held in the arrays.
     */
    private class Node implements NodeRef {
        private final int number;
        private final Taxon taxon;

        Node(int number, Taxon taxon) {
            this.number = number;
            this.taxon = taxon;
        }

        public int getNumber() {
            return number;
        }

        public void setNumber(int n) {
            throw new UnsupportedOperationException("FastTreeModel nodes cannot be renumbered");
        }

        public String toString() {
            return "node " + number + ", height=" + heights[number] + (taxon != null ? ": " + taxon.getId() : "");
        }
    }

    /**
     * The height parameter of one slot (one of DefaultTreeModel's node height parameter objects): a
     * view of the height of the node the slot is at. It holds no value and stores nothing (the tree
     * stores the heights).
     */
    private class SlotHeightParameter extends Parameter.Abstract {

        final int slot;
        private Bounds<Double> bounds = null;

        SlotHeightParameter(int slot) {
            this.slot = slot;
        }

        public double getParameterValue(int dim) {
            return heights[slotNode[slot]];
        }

        public void setParameterValue(int dim, double value) {
            setHeight(value);
            fireParameterChangedEvent(dim, Parameter.ChangeType.VALUE_CHANGED);
        }

        public void setParameterValueQuietly(int dim, double value) {
            setHeight(value);
        }

        public void setParameterValueNotifyChangedAll(int dim, double value) {
            setHeight(value);
            fireParameterChangedEvent(-1, Parameter.ChangeType.ALL_VALUES_CHANGED);
        }

        private void setHeight(double value) {
            final int n = slotNode[slot];
            saveNode(n);
            heights[n] = value;
        }

        public String getParameterName() {
            return getId();
        }

        public void addBounds(Bounds<Double> boundary) {
            // as Parameter.Default
            if (bounds == null) {
                bounds = boundary;
            } else {
                if (!(bounds instanceof HeightIntersectionBounds)) {
                    HeightIntersectionBounds newBounds = new HeightIntersectionBounds();
                    newBounds.addBounds(bounds);
                    bounds = newBounds;
                }

                ((HeightIntersectionBounds) bounds).addBounds(boundary);
            }
        }

        public Bounds<Double> getBounds() {
            if (bounds == null) {
                throw new NullPointerException(getParameterName() + " parameter: Bounds not set");
            }
            return bounds;
        }

        public void addDimension(int index, double value) {
            throw new UnsupportedOperationException("node height parameters have one dimension");
        }

        public double removeDimension(int index) {
            throw new UnsupportedOperationException("node height parameters have one dimension");
        }

        protected void storeValues() {
            // the tree stores the heights
        }

        protected void restoreValues() {
            // the tree restores the heights
        }

        protected void acceptValues() {
        }

        protected void adoptValues(Parameter source) {
            setHeight(source.getParameterValue(0));
        }
    }

    /**
     * The intersection of several bounds on a height, as dr.inference.model.IntersectionBounds (whose
     * constructor is not public): the highest lower limit and the lowest upper limit.
     */
    private static class HeightIntersectionBounds implements Bounds<Double> {

        private final java.util.List<Bounds<Double>> bounds = new java.util.ArrayList<Bounds<Double>>();

        void addBounds(Bounds<Double> boundary) {
            if (boundary.getBoundsDimension() != 1) {
                throw new IllegalArgumentException("Incorrect dimension of bounds, expected 1 but received " +
                        boundary.getBoundsDimension());
            }
            bounds.add(boundary);
        }

        public Double getLowerLimit(int index) {
            double lower = Double.NEGATIVE_INFINITY;
            for (Bounds<Double> boundary : bounds) {
                if (boundary.getLowerLimit(index) > lower) {
                    lower = boundary.getLowerLimit(index);
                }
            }
            return lower;
        }

        public Double getUpperLimit(int index) {
            double upper = Double.POSITIVE_INFINITY;
            for (Bounds<Double> boundary : bounds) {
                if (boundary.getUpperLimit(index) < upper) {
                    upper = boundary.getUpperLimit(index);
                }
            }
            return upper;
        }

        public int getBoundsDimension() {
            return 1;
        }
    }

    /**
     * Bounds of a slot's height: the height of the parent of its node above (none for the root), the
     * higher of its children below (0 for a tip), as DefaultTreeModel.NodeHeightBounds.
     */
    private class NodeHeightBounds implements Bounds<Double> {

        private final int slot;

        NodeHeightBounds(int slot) {
            this.slot = slot;
        }

        public Double getUpperLimit(int i) {
            final int n = slotNode[slot];
            if (parent[n] == -1) {
                return Double.POSITIVE_INFINITY;
            } else {
                return heights[parent[n]];
            }
        }

        public Double getLowerLimit(int i) {
            final int n = slotNode[slot];
            if (child0[n] == -1 && child1[n] == -1) {
                return 0.0;
            } else {
                return Math.max(heights[child0[n]], heights[child1[n]]);
            }
        }

        public int getBoundsDimension() {
            return 1;
        }
    }

    /**
     * A compound parameter of slot height parameters that finds the index of a changed slot directly
     * (rather than by searching the components).
     */
    private class NodeHeightsParameter extends CompoundParameter {

        private final int[] slotIndex = new int[nodeCount];

        NodeHeightsParameter(String name) {
            super(name);
            Arrays.fill(slotIndex, -1);
        }

        @Override
        public void addParameter(Parameter param) {
            super.addParameter(param);
            if (param instanceof SlotHeightParameter && slotIndex[((SlotHeightParameter) param).slot] == -1) {
                // all components have one dimension, so the offset is the position
                slotIndex[((SlotHeightParameter) param).slot] = getParameterCount() - 1;
            } else {
                indexed = false;
            }
        }

        @Override
        public void removeParameter(Parameter param) {
            super.removeParameter(param);
            indexed = false;
        }

        @Override
        public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            if (indexed && variable instanceof SlotHeightParameter) {
                if (!doNotPropagateChangeUp) {
                    final int k = slotIndex[((SlotHeightParameter) variable).slot];
                    if (k != -1) {
                        fireParameterChangedEvent(index == -1 ? -1 : k + index, type);
                    }
                }
            } else {
                super.variableChangedEvent(variable, index, type);
            }
        }

        private boolean indexed = true;
    }

    // ***********************************************************************
    // Private members
    // ***********************************************************************

    private final int nodeCount;
    private final int externalNodeCount;
    private final int internalNodeCount;

    private final Node[] nodes;

    // the tree
    private final double[] heights;
    private final int[] parent;
    private final int[] child0;
    private final int[] child1;
    private int root;
    private int oldRoot;

    // the height parameter slot at each node, and the node of each slot
    private final int[] nodeSlot;
    private final int[] slotNode;
    private final SlotHeightParameter[] slotParameters;

    // copy-on-write store: the nodes changed since the last store, as they were then
    private final double[] savedHeight;
    private final int[] savedParent;
    private final int[] savedChild0;
    private final int[] savedChild1;
    private final int[] savedSlot;
    private final boolean[] isChanged;
    private final int[] changedNodes;
    private int changedNodeCount = 0;
    private int storedRoot;

    private boolean isTipDateSampled = false;
}
