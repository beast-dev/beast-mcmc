/*
 * Tree.java
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
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.util.Attributable;
import dr.util.Identifiable;
import jebl.evolution.graphs.Node;
import jebl.evolution.trees.SimpleRootedTree;

import java.text.NumberFormat;
import java.util.*;

/**
 * Interface for a phylogenetic or genealogical tree.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: Tree.java,v 1.59 2006/09/08 17:34:23 rambaut Exp $
 */
public interface Tree extends TaxonList, Units, Identifiable, Attributable {

    /**
     * @return root node of this tree.
     */
    NodeRef getRoot();

    /**
     * @return a count of the number of nodes (internal + external) in this
     *         tree, currently connected from the root node.
     */
    int getNodeCount();

    /**
     * @param i node index, terminal nodes are first
     * @return the ith node.
     */
    NodeRef getNode(int i);

    /**
     * @param i index of an internal node
     * @return the ith internal node.
     */
    NodeRef getInternalNode(int i);

    /**
     * @param i the index of an external node
     * @return the ith external node.
     */
    NodeRef getExternalNode(int i);

    /**
     * @return a count of the number of external nodes (tips) in this
     *         tree, currently connected from the root node.
     */
    int getExternalNodeCount();

    /**
     * @return a count of the number of internal nodes in this
     *         tree, currently connected from the root node.
     */
    int getInternalNodeCount();

    /**
     * @param node the node to retrieve the taxon of
     * @return the taxon of this node.
     */
    Taxon getNodeTaxon(NodeRef node);

    /**
     * @return whether this tree has known node heights.
     */
    boolean hasNodeHeights();

    /**
     * @param node the node to retrieve height of
     * @return the height of node in the tree.
     */
    double getNodeHeight(NodeRef node);

    /**
     * @return whether this tree has known branch lengths.
     */
    boolean hasBranchLengths();

    /**
     * @param node the node to retrieve the length of branch to parent
     * @return the length of the branch from node to its parent.
     */
    double getBranchLength(NodeRef node);

    /**
     * @param node the node to retrieve the rate of
     * @return the rate of node in the tree.
     */
    double getNodeRate(NodeRef node);

    /**
     * @param node the node whose attribute is being fetched.
     * @param name the name of the attribute of interest.
     * @return an object representing the named attributed for the given node.
     */
    Object getNodeAttribute(NodeRef node, String name);

    /**
     * @param node the node whose attribute is being fetched.
     * @return an interator of attribute names available for this node.
     */
    Iterator getNodeAttributeNames(NodeRef node);

    /**
     * @param node the node to test if external
     * @return whether the node is external.
     */
    boolean isExternal(NodeRef node);

    /**
     * @param node the node to test if root
     * @return whether the node is the root.
     */
    boolean isRoot(NodeRef node);

    /**
     * @param node the node to get child count of
     * @return the number of children of node.
     */
    int getChildCount(NodeRef node);

    /**
     * @param node the node to get jth child of
     * @param j    the index of child to retrieve
     * @return the jth child of node
     */
    NodeRef getChild(NodeRef node, int j);

    NodeRef getParent(NodeRef node);

    /**
     * @return a clone of this tree
     */
    public Tree getCopy();

}