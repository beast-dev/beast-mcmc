/*
 * TreeModel.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evolution.util.Taxon;
import dr.inference.model.*;

import java.util.*;

/**
 * A model component for trees.
 *
 * @author Andrew Rambaut
 * @version $Id:$
 */
public class BigFastTreeModel extends AbstractTreeModel {

    //
    // Public stuff
    //

    public static final String BIG_FAST_TREE_MODEL = "bigFastTreeModel";

    private static final boolean TEST_NODE_BOUNDS = false;

    public BigFastTreeModel(Tree tree) {
        this(BIG_FAST_TREE_MODEL, tree, false, false);
    }

    public BigFastTreeModel(String name, Tree tree) { this(name, tree, false, false); }

    /* New constructor that copies the attributes of Tree tree into the new TreeModel
      * Useful for constructing a TreeModel from a NEXUS file entry
      */

    public BigFastTreeModel(String name, Tree tree, boolean fixHeights, boolean fixTree) {

        super(name, !fixTree);
        setId(name);

        // get a rooted version of the tree to clone
        FlexibleTree binaryTree = new FlexibleTree(tree);
        binaryTree.resolveTree();

        // adjust the heights to be compatible with the tip dates and perturb
        // any zero branches.
        if (!fixHeights) {
            MutableTree.Utils.correctHeightsForTips(binaryTree);
        }

        internalNodeCount = binaryTree.getInternalNodeCount();
        externalNodeCount = binaryTree.getExternalNodeCount();

        nodeCount = internalNodeCount + externalNodeCount;

        nodes = new NodeRef[nodeCount];

        edges = new int[nodeCount * 3];
        storedEdges = new int[nodeCount * 3];

        heights = new double[nodeCount];
        storedHeights = new double[nodeCount];

        int i = 0;
        int j = externalNodeCount;

        boolean done = false;
        do {
            NodeRef node = TreeUtils.postorderSuccessor(binaryTree, binaryTree.getRoot());

            if (binaryTree.isExternal(node)) {
                nodes[i] = new Node(i, tree.getNodeTaxon(node));
                edges[(j * 3)] = -1; // parent
                edges[(i * 3) + 1] = -1; // child 1
                edges[(i * 3) + 2] = -1; // child 2

                i++;
            } else {
                if (binaryTree.isRoot(node)) {
                    done = true;
                }
                nodes[j] = new Node(j);

                nodes[j] = node;
                edges[(j * 3)] = -1; // parent
                edges[(j * 3) + 1] = -1; // child 1
                edges[(j * 3) + 2] = -1; // child 2

                j++;
            }
        } while (!done);
        
        root = j;

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


    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // no submodels so nothing to do
    }

    /**
     * Called when a parameter changes.
     */
    public void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
//        final Node node = getNodeOfParameter((Parameter) variable);
//        if (type == Parameter.ChangeType.ALL_VALUES_CHANGED) {
//            //this signals events where values in all dimensions of a parameter is changed.
//            pushTreeChangedEvent(new TreeChangedEvent(node, (Parameter) variable, TreeChangedEvent.CHANGE_IN_ALL_INTERNAL_NODES));
//        } else {
//            pushTreeChangedEvent(node, (Parameter) variable, index);
//        }
    }


    public boolean inTreeEdit() {
        return inEdit;
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
        throw new UnsupportedOperationException("getNodeRate not available in BigFastTreeModel");
    }

    public Object getNodeAttribute(NodeRef node, String name) {
        throw new UnsupportedOperationException("getNodeAttribute not available in BigFastTreeModel");
    }

    public Iterator getNodeAttributeNames(NodeRef node) {
        throw new UnsupportedOperationException("getNodeAttributeNames not available in BigFastTreeModel");
    }

    public double[] getMultivariateNodeTrait(NodeRef node, String name) {
        throw new UnsupportedOperationException("getMultivariateNodeTrait not available in BigFastTreeModel");
    }

    public Taxon getNodeTaxon(NodeRef node) {
        return ((Node) node).taxon;
    }

    public boolean isExternal(NodeRef node) {
        return ((Node) node).isExternal();
    }

    public boolean isRoot(NodeRef node) {
        return ((Node) node).number == root;
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
        return nodes[root];
    }

    @Override
    public boolean isTipDateSampled() {
        return false;
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

    public boolean beginTreeEdit() {
        if (inEdit) throw new RuntimeException("Already in edit transaction mode!");

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

        for (dr.evomodel.tree.TreeChangedEvent treeChangedEvent : treeChangedEvents) {
            listenerHelper.fireModelChanged(this, treeChangedEvent);
        }
        treeChangedEvents.clear();
    }

    public void checkTreeIsValid() throws InvalidTreeException {
        for (Node node : nodes) {
            if (!node.heightParameter.isWithinBounds()) {
                throw new InvalidTreeException("height parameter out of bounds");
            }
        }
    }

    public void setNodeHeight(NodeRef n, double height) {
        ((Node) n).setHeight(height);
    }

    public void setNodeHeightQuietly(NodeRef n, double height) {
        ((Node) n).setHeightQuietly(height);
    }

    public void setNodeRate(NodeRef n, double rate) {
        throw new UnsupportedOperationException("BigFastTreeModel cannot have node rates set");
    }

    public void setNodeTrait(NodeRef n, String name, double value) {
        throw new UnsupportedOperationException("BigFastTreeModel cannot have traits set");
    }

    public void setMultivariateTrait(NodeRef n, String name, double[] value) {
        throw new UnsupportedOperationException("BigFastTreeModel cannot have traits set");
    }

    public void setBranchLength(NodeRef node, double length) {
        throw new UnsupportedOperationException("BigFastTreeModel cannot have branch lengths set");
    }

    public void setNodeAttribute(NodeRef node, String name, Object value) {
        throw new UnsupportedOperationException("BigFastTreeModel does not use NodeAttributes");
    }

    // *****************************************************************
    // Interface ModelComponent
    // *****************************************************************

    /**
     * Store current state
     */
    protected void storeState() {
        System.arraycopy(edges, 0, storedEdges, 0, edges.length);
        System.arraycopy(heights, 0, storedHeights, 0, heights.length);

        storedRoot = root;

    }

    /**
     * Restore the stored state
     */
    protected void restoreState() {

        int[] tmp = edges;
        storedEdges = edges;
        edges = tmp;

        double[] tmp2 = heights;
        storedHeights = heights;
        heights = tmp2;

        root = storedRoot;
    }

    /**
     * accept the stored state
     */
    protected void acceptState() {
    } // nothing to do

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

    // **************************************************************
    // Private inner classes
    // **************************************************************

    private class Node implements NodeRef {
        private Node(int number) {
            this(number, null);
        }

        private Node(int number, Taxon taxon) {
            this.number = number;
            this.taxon = taxon;
        }

        @Override
        public int getNumber() {
            return number;
        }

        @Override
        public void setNumber(int n) {
            throw new UnsupportedOperationException("Node is immutable");
        }

        public Taxon getTaxon() {
            return taxon;
        }

        final Taxon taxon;
        final int number;
    }
    // ***********************************************************************
    // Private members
    // ***********************************************************************


    /**
     * root node
     */
    private int root = -1;
    private int storedRoot;

    /**
     * list of internal nodes (including root)
     */
    private int[] edges = null;
    private int[] storedEdges = null;

    private double[] heights = null;
    private double[] storedHeights = null;

    private NodeRef[] nodes;

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

    private boolean inEdit = false;


}
