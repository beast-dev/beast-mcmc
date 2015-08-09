/*
 * FlexibleNode.java
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

import dr.evolution.util.Taxon;
import dr.util.Attributable;

import java.util.Iterator;

/**
 * A flexible implementation of the Node interface. Can have both node heights
 * and node lengths. Note that the height and length are not garenteed to both
 * be correct. Which is currently in use and the conversion between them is handled
 * by FlexibleTree.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: FlexibleNode.java,v 1.9 2005/06/20 16:06:02 rambaut Exp $
 */
public class FlexibleNode implements NodeRef, Attributable {

	/**
	 * parent node
	 */
	private FlexibleNode parent;

	/**
	 * number of node
	 */
	private int nodeNumber;

	/**
	 * height of this node
	 */
	private double height;

	/**
	 * length of the branch to the parent
	 */
	private double length;

	/**
	 * instantaneous rate at this node
	 */
	private double rate;

	/**
	 * the taxon if node is external
	 */
	private Taxon taxon = null;

	//
	// Private stuff
	//

	private FlexibleNode[] child;

	/**
	 * constructor default node
	 */
	public FlexibleNode() {
		parent = null;
		child = null;
		height = 0.0;
		rate = 1.0;

		nodeNumber = 0;
	}

	/**
	 * constructor default node
	 */
	public FlexibleNode(Taxon taxon) {
		parent = null;
		child = null;
		height = 0.0;
		rate = 1.0;

		nodeNumber = 0;
		this.taxon = taxon;
	}

	public FlexibleNode(FlexibleNode node) {
		parent = null;

		setHeight(node.getHeight());
		setLength(node.getLength());
		setRate(node.getRate());
		setId(node.getId());
		setNumber(node.getNumber());
		setTaxon(node.getTaxon());

		child = null;

		for (int i = 0; i < node.getChildCount(); i++) {
			addChild(new FlexibleNode(node.getChild(i)));
		}
	}

	/**
	 * constructor used to clone a node and all children
	 */
	public FlexibleNode(Tree tree, NodeRef node) {
        this(tree, node, false);
        // Unbelievable code duplication

//        parent = null;
//		setHeight(tree.getNodeHeight(node));
//		setLength(tree.getBranchLength(node));
//		setRate(tree.getNodeRate(node));
//		setId(tree.getTaxonId(node.getNumber()));
//		setNumber(node.getNumber());
//		setTaxon(tree.getNodeTaxon(node));
//
//		child = null;
//
//		for (int i = 0; i < tree.getChildCount(node); i++) {
//			addChild(new FlexibleNode(tree, tree.getChild(node, i)));
//		}
	}

	public FlexibleNode(Tree tree, NodeRef node, boolean copyAttributes) {
		parent = null;
		setHeight(tree.getNodeHeight(node));
		setLength(tree.getBranchLength(node));
		setRate(tree.getNodeRate(node));
		setId(tree.getTaxonId(node.getNumber()));
		setNumber(node.getNumber());
		setTaxon(tree.getNodeTaxon(node));

		child = null;

		for (int i = 0; i < tree.getChildCount(node); i++) {
			addChild(new FlexibleNode(tree, tree.getChild(node, i), copyAttributes));
		}
		if (copyAttributes) {
			Iterator iter = tree.getNodeAttributeNames(node);
			if (iter != null) {
				while (iter.hasNext()) {
					String name = (String) iter.next();
					this.setAttribute(name, tree.getNodeAttribute(node, name));
				}
			}
		}
	}

	public FlexibleNode getDeepCopy() {
        return new FlexibleNode(this);
	}

	public FlexibleNode getShallowCopy() {
		FlexibleNode copy = new FlexibleNode();

		copy.setHeight(getHeight());
		copy.setLength(getLength());
		copy.setRate(getRate());
		copy.setId(getId());
		copy.setNumber(getNumber());
		copy.setTaxon(getTaxon());

		return copy;
	}

	/**
	 * Returns the parent node of this node.
	 */
	public final FlexibleNode getParent() {
		return parent;
	}

	/**
	 * Set the parent node of this node.
	 */
	public void setParent(FlexibleNode node) {
		parent = node;
	}

	/**
	 * Get the height of this node.
	 */
	public final double getHeight() {
		return height;
	}

	/**
	 * Set the height of this node.
	 */
	public final void setHeight(double value) {
		height = value;
	}

	/**
	 * Get the length of the branch to the node's parent.
	 */
	public final double getLength() {
		return length;
	}

	/**
	 * Set the length of the branch to the node's parent.
	 */
	public final void setLength(double value) {
		length = value;
	}

	/**
	 * Get the rate at this node.
	 */
	public final double getRate() {
		return rate;
	}

	/**
	 * Set the rate at this node.
	 */
	public final void setRate(double value) {
		rate = value;
	}

	public void setNumber(int n) {
		nodeNumber = n;
	}

	public int getNumber() {
		return nodeNumber;
	}

	public void setTaxon(Taxon taxon) {
		this.taxon = taxon;
	}

	public Taxon getTaxon() {
		return taxon;
	}

	/**
	 * get child node
	 *
	 * @param n number of child
	 * @return child node
	 */
	public FlexibleNode getChild(int n) {
		return child[n];
	}

	public boolean hasChild(FlexibleNode node) {

		for (int i = 0, n = getChildCount(); i < n; i++) {
			if (node == child[i]) return true;
		}
		return false;
	}

	/**
	 * add new child node
	 *
	 * @param n new child node
	 */
	public void addChild(FlexibleNode n) {
		insertChild(n, getChildCount());
	}

	/**
	 * add new child node (insertion at a specific position)
	 *
	 * @param n new child node
	 *          + @param pos position
	 */
	public void insertChild(FlexibleNode n, int pos) {
		int numChildren = getChildCount();

		FlexibleNode[] newChild = new FlexibleNode[numChildren + 1];

// AR 28/05/08
// This doesn't work because pos can be zero and (numChildren - pos) can be zero in which case an NullPointerException is thrown.
// We could check for these special cases but we are generally only moving one or two nodes so it is unlikely to be an efficiency
// gain for using arraycopy.
//        System.arraycopy(child, 0, newChild, 0, pos);
//		newChild[pos] = n;
//        System.arraycopy(child, pos, newChild, pos + 1, numChildren - pos);

        for (int i = 0; i < pos; i++) {
            newChild[i] = child[i];
        }
        newChild[pos] = n;
        for (int i = pos; i < numChildren; i++) {
            newChild[i + 1] = child[i];
        }

		child = newChild;

		n.setParent(this);
	}

	/**
	 * remove child
	 *
	 * @param n child to be removed
	 */
	public FlexibleNode removeChild(FlexibleNode n) {
		int numChildren = getChildCount();
		FlexibleNode[] newChild = new FlexibleNode[numChildren - 1];

		int j = 0;
		boolean found = false;
		for (int i = 0; i < numChildren; i++) {
			if (child[i] != n) {
				newChild[j] = child[i];
				j++;
			} else
				found = true;
		}

		if (!found)
			throw new IllegalArgumentException("Nonexistent child");

		//remove parent link from removed child!
		n.setParent(null);

		child = newChild;

		return n;
	}

	/**
	 * remove child
	 *
	 * @param n number of child to be removed
	 */
	public FlexibleNode removeChild(int n) {
		int numChildren = getChildCount();

		if (n >= numChildren) {
			throw new IllegalArgumentException("Nonexistent child");
		}

		return removeChild(child[n]);
	}


	/**
	 * check whether this node has any children
	 *
	 * @return result (true or false)
	 */
	public boolean hasChildren() {
		return (getChildCount() != 0);
	}

	/**
	 * check whether this node is an external node
	 *
	 * @return result (true or false)
	 */
	public boolean isExternal() {
		return !hasChildren();
	}

	/**
	 * check whether this node is a root node
	 *
	 * @return result (true or false)
	 */
	public boolean isRoot() {
		return (getParent() == null);
	}


	/**
	 * Returns the number of children this node has.
	 */
	public final int getChildCount() {
		if (child == null) return 0;
		return child.length;
	}

	// **************************************************************
	// Identifiable IMPLEMENTATION
	// **************************************************************

	private String id = null;

	/**
	 * @return the id as a string.
	 */
	public String getId() {
		return id;
	}

	/**
	 * set the id as a string.
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
	 *
	 * @param name  the name of the attribute.
	 * @param value the new value of the attribute.
	 */
	public void setAttribute(String name, Object value) {
		if (attributes == null)
			attributes = new Attributable.AttributeHelper();
		attributes.setAttribute(name, value);
	}

	/**
	 * @param name the name of the attribute of interest.
	 * @return an object representing the named attributed for this object.
	 */
	public Object getAttribute(String name) {

		if (attributes == null) {
			return null;
		} else {
			return attributes.getAttribute(name);
		}
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
}
