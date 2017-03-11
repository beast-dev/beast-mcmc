/*
 * SimpleTree.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evolution.tree;

import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.util.Attributable;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * A data structure for binary rooted trees.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: SimpleTree.java,v 1.65 2006/08/23 10:46:33 rambaut Exp $
 */
public class SimpleTree implements MutableTree {

    /** Constructor tree with no nodes. Use adoptNodes to add some nodes. */
    public SimpleTree() {
        root = null;
    }

    /** clone constructor */
    public SimpleTree(Tree tree) {

        setUnits(tree.getUnits());

        root = new SimpleNode(tree, tree.getRoot());

        nodeCount = tree.getNodeCount();
        internalNodeCount = tree.getInternalNodeCount();
        externalNodeCount = tree.getExternalNodeCount();

        nodes = new SimpleNode[nodeCount];

        SimpleNode node = root;
        do {
            node = (SimpleNode) TreeUtils.postorderSuccessor(this, node);
            if ((node.getNumber() >= externalNodeCount && node.isExternal()) ||
                (node.getNumber() < externalNodeCount && !node.isExternal())) {
                throw new RuntimeException("Error cloning tree: node numbers are incompatible");
            }
            nodes[node.getNumber()] = node;
        } while (node != root);
    }

    /** clone constructor */
    public SimpleTree(SimpleNode root) {
        adoptNodes(root);
    }

    /**
     * @return a copy of this tree
     */
    public Tree getCopy() {
        return new SimpleTree(this);
    }

    /**
     * Adopt a node hierarchy as its own. Only called by the SimpleTree(SimpleNode, TaxonList).
     * This creates the node list and stores the nodes in post-traversal order.
     */
    protected void adoptNodes(SimpleNode node) {

        if (inEdit) throw new RuntimeException("Mustn't be in an edit transaction to call this method!");

        internalNodeCount = 0;
        externalNodeCount = 0;

        root = node;

        do {
            node = (SimpleNode) TreeUtils.postorderSuccessor(this, node);
            if (node.isExternal()) {
                externalNodeCount++;
            } else
                internalNodeCount++;

        } while(node != root);

        nodeCount = internalNodeCount + externalNodeCount;

        nodes = new SimpleNode[nodeCount];

        node = root;
        int i = 0;
        int j = externalNodeCount;

        do {
            node = (SimpleNode) TreeUtils.postorderSuccessor(this, node);
            if (node.isExternal()) {
                node.setNumber(i);
                nodes[i] = node;
                i++;
            } else {
                node.setNumber(j);
                nodes[j] = node;
                j++;
            }
        } while(node != root);
    }

    /**
     * Return the units that this tree is expressed in.
     */
    public final Type getUnits() {
        return units;
    }

    /**
     * Sets the units that this tree is expressed in.
     */
    public final void setUnits(Type units) {
        this.units = units;
    }

    /**
     * @return a count of the number of nodes (internal + external) in this
     * tree.
     */
    public int getNodeCount() {
        return nodeCount;
    }

    public boolean hasNodeHeights() { return true; }
    public double getNodeHeight(NodeRef node) { return ((SimpleNode)node).getHeight(); }
    public double getNodeRate(NodeRef node) { return ((SimpleNode)node).getRate(); }
    public Taxon getNodeTaxon(NodeRef node) { return ((SimpleNode)node).getTaxon(); }
    public int getChildCount(NodeRef node) { return ((SimpleNode)node).getChildCount(); }
    public boolean isExternal(NodeRef node) { return ((SimpleNode)node).getChildCount() == 0; }
    public boolean isRoot(NodeRef node) { return (node == root); }
    public NodeRef getChild(NodeRef node, int i) { return ((SimpleNode)node).getChild(i); }
    public NodeRef getParent(NodeRef node) { return ((SimpleNode)node).getParent(); }

    public boolean hasBranchLengths() { return true; }
    public double getBranchLength(NodeRef node) {
        NodeRef parent = getParent(node);
        if (parent == null) {
            return 0.0;
        }

        return getNodeHeight(parent) - getNodeHeight(node);
    }
    public void setBranchLength(NodeRef node, double length) {
        throw new UnsupportedOperationException("SimpleTree cannot have branch lengths set... use FlexibleTree");
    }


    public final SimpleNode getExternalNode(int i) {
        return nodes[i];
    }

    public final SimpleNode getInternalNode(int i) {
        return nodes[i+externalNodeCount];
    }

    public final NodeRef getNode(int i) {
        return nodes[i];
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

    /**
     * Returns the root node of this tree.
     */
    public final NodeRef getRoot() {
        return root;
    }

    /**
     * Set a new node as root node.
     */
    public final void setRoot(NodeRef r) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        if (!(r instanceof SimpleNode)) { throw new IllegalArgumentException(); }
        root = (SimpleNode)r;
        root.setParent(null);
    }


    /**
     * @return the height of the root node.
     */
    public final double getRootHeight() {
        return root.getHeight();
    }

    /**
     * Set the height of the root node.
     */
    public final void setRootHeight(double height) {
        root.setHeight(height);
        fireTreeChanged();
    }


    public void addChild(NodeRef p, NodeRef c) {
        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");
        SimpleNode parent = (SimpleNode)p;
        SimpleNode child = (SimpleNode)c;
        if (parent.hasChild(child)) throw new IllegalArgumentException("Child already existists in parent");

        parent.addChild(child);
    }

    public void removeChild(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        SimpleNode parent = (SimpleNode)p;
        SimpleNode child = (SimpleNode)c;

        parent.removeChild(child);
    }

    public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {
        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");
        SimpleNode parent = (SimpleNode)node;
        parent.replaceChild((SimpleNode)child, (SimpleNode)newChild);
    }

    public boolean beginTreeEdit() {
        boolean r = inEdit;
        inEdit = true;
        return r;
    }

    public void endTreeEdit() {
        inEdit = false;
        fireTreeChanged();
    }

    public void setNodeHeight(NodeRef n, double height) {
        SimpleNode node = (SimpleNode)n;
        node.setHeight(height);
        fireTreeChanged();
    }

    public void setNodeRate(NodeRef n, double rate) {
        SimpleNode node = (SimpleNode)n;
        node.setRate(rate);
        fireTreeChanged();
    }

    /**
     * Sets an named attribute for a given node.
     * @param node the node whose attribute is being set.
     * @param name the name of the attribute.
     * @param value the new value of the attribute.
     */
    public void setNodeAttribute(NodeRef node, String name, Object value) {
        ((SimpleNode)node).setAttribute(name, value);
        fireTreeChanged();
    }

    /**
     * @return an object representing the named attributed for the given node.
     * @param node the node whose attribute is being fetched.
     * @param name the name of the attribute of interest.
     */
    public Object getNodeAttribute(NodeRef node, String name) {
        return ((SimpleNode)node).getAttribute(name);
    }

    /**
     * @return an interator of attribute names available for this node.
     * @return a key set of attribute names available for this node.
     */
    public Iterator getNodeAttributeNames(NodeRef node) {
        return ((SimpleNode)node).getAttributeNames();
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
        return getExternalNode(taxonIndex).getTaxon();
    }

    /**
     * @return the ID of the taxon of the ith external node. If it doesn't have
     * a taxon, returns the ID of the node itself.
     */
    public String getTaxonId(int taxonIndex) {
        Taxon taxon = getTaxon(taxonIndex);
        if (taxon != null)
            return taxon.getId();
        else
            return getExternalNode(taxonIndex).getId();
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
                index ++;
                return getTaxon(index);
            }

            public void remove() { /* do nothing */ }
        };
    }

    /**
     * @return an object representing the named attributed for the taxon of the given
     * external node. If the node doesn't have a taxon then the nodes own attribute
     * is returned.
     * @param taxonIndex the index of the taxon whose attribute is being fetched.
     * @param name the name of the attribute of interest.
     */
    public Object getTaxonAttribute(int taxonIndex, String name) {
        Taxon taxon = getTaxon(taxonIndex);
        if (taxon != null)
            return taxon.getAttribute(name);
        else
            return getExternalNode(taxonIndex).getAttribute(name);
    }

    // **************************************************************
    // MutableTaxonList IMPLEMENTATION
    // **************************************************************

    public int addTaxon(Taxon taxon) { throw new IllegalArgumentException("Cannot add taxon to a MutableTree"); }
    public boolean removeTaxon(Taxon taxon) { throw new IllegalArgumentException("Cannot add taxon to a MutableTree"); }

    /**
     * Sets the ID of the taxon of the ith external node. If it doesn't have
     * a taxon, sets the ID of the node itself.
     */
    public void setTaxonId(int taxonIndex, String id) {
        Taxon taxon = getTaxon(taxonIndex);
        if (taxon != null)
            taxon.setId(id);
        else
            getExternalNode(taxonIndex).setId(id);

        fireTreeChanged();
        fireTaxaChanged();
    }

    /**
     * Sets an named attribute for the taxon of a given external node. If the node
     * doesn't have a taxon then the attribute is added to the node itself.
     * @param taxonIndex the index of the taxon whose attribute is being set.
     * @param name the name of the attribute.
     * @param value the new value of the attribute.
     */
    public void setTaxonAttribute(int taxonIndex, String name, Object value) {
        Taxon taxon = getTaxon(taxonIndex);
        if (taxon != null)
            taxon.setAttribute(name, value);
        else
            getExternalNode(taxonIndex).setAttribute(name, value);

        fireTreeChanged();
        fireTaxaChanged();
    }

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

    private Attributable.AttributeHelper attributes = null;

    /**
     * Sets an named attribute for this object.
     * @param name the name of the attribute.
     * @param value the new value of the attribute.
     */
    public void setAttribute(String name, Object value) {
        if (attributes == null)
            attributes = new Attributable.AttributeHelper();
        attributes.setAttribute(name, value);
    }

    /**
     * @return an object representing the named attributed for this object.
     * @param name the name of the attribute of interest.
     */
    public Object getAttribute(String name) {
        if (attributes == null)
            return null;
        else
            return attributes.getAttribute(name);
    }

    /**
     * @return an iterator of the attributes that this object has.
     */
    public Iterator<String> getAttributeNames() {
        if (attributes == null)
            return null;
        else
            return attributes.getAttributeNames();
    }

    public void addMutableTreeListener(MutableTreeListener listener) {
        mutableTreeListeners.add(listener);
    }

    private void fireTreeChanged() {
        for (MutableTreeListener mutableTreeListener : mutableTreeListeners) {
            mutableTreeListener.treeChanged(this);
        }
    }

    private final ArrayList<MutableTreeListener> mutableTreeListeners = new ArrayList<MutableTreeListener>();

    public void addMutableTaxonListListener(MutableTaxonListListener listener) {
        mutableTaxonListListeners.add(listener);
    }

    private void fireTaxaChanged() {
        for (MutableTaxonListListener mutableTaxonListListener : mutableTaxonListListeners) {
            mutableTaxonListListener.taxaChanged(this);
        }
    }

    private final ArrayList<MutableTaxonListListener> mutableTaxonListListeners = new ArrayList<MutableTaxonListListener>();

    /**
     * @return a string containing a newick representation of the tree
     */
    public String toString() {
        return TreeUtils.newick(this);
    }

    /**
     * @return whether two trees have the same topology
     */
    public boolean equals(Object obj) {
        if (!(obj instanceof Tree)) {
            throw new IllegalArgumentException("SimpleTree.equals can only compare instances of Tree");
        }
        return TreeUtils.equal(this, (Tree)obj);
    }

    // **************************************************************
    // Private members
    // **************************************************************

    /** root node */
    SimpleNode root;

    /** list of internal nodes (including root) */
    SimpleNode[] nodes = null;

    /** number of nodes (including root and tips) */
    int nodeCount;

    /** number of external nodes */
    int externalNodeCount;

    /** number of internal nodes (including root) */
    int internalNodeCount;

    /** holds the units of the trees branches. */
    private Type units = Type.SUBSTITUTIONS;

    boolean inEdit = false;
}
