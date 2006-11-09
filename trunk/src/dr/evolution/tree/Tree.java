/*
 * Tree.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

import dr.evolution.colouring.BranchColouring;
import dr.evolution.colouring.TreeColouring;
import dr.evolution.util.*;
import dr.util.Attributable;
import dr.util.Identifiable;

import java.util.*;

/**
 * Interface for a phylogenetic or genealogical tree.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: Tree.java,v 1.59 2006/09/08 17:34:23 rambaut Exp $
 */
public interface Tree extends TaxonList, Units, Identifiable, Attributable {

    /**
     * @return the index of root node of this tree.
     */
    NodeRef getRoot();

    /**
     * @return a count of the number of nodes (internal + external) in this
     * tree, currently connected from the root node.
     */
    int getNodeCount();

    /**
     * @return the ith node.
     */
    NodeRef getNode(int i);

    /**
     * @return the ith internal node.
     */
    NodeRef getInternalNode(int i);

    /**
     * @return the ith internal node.
     */
    NodeRef getExternalNode(int i);

    /**
     * @return a count of the number of external nodes (tips) in this
     * tree, currently connected from the root node.
     */
    int getExternalNodeCount();

    /**
     * @return a count of the number of internal nodes in this
     * tree, currently connected from the root node.
     */
    int getInternalNodeCount();

    /**
     * @return the taxon of this node.
     */
    Taxon getNodeTaxon(NodeRef node);

    /**
     * @return whether this tree has known node heights.
     */
    boolean hasNodeHeights();

    /**
     * @return the height of the ith node in the tree (where the first n are internal).
     */
    double getNodeHeight(NodeRef node);

    /**
     * @return whether this tree has known branch lengths.
     */
    boolean hasBranchLengths();

    /**
     * @return the length of the branch from the ith node in the tree to its parent.
     */
    double getBranchLength(NodeRef node);

    /**
     * @return the rate of the ith node in the tree (where the first n are internal).
     */
    double getNodeRate(NodeRef node);

    /**
     * @return an object representing the named attributed for the given node.
     * @param node the node whose attribute is being fetched.
     * @param name the name of the attribute of interest.
     */
    Object getNodeAttribute(NodeRef node, String name);

    /**
     * @return an interator of attribute names available for this node.
     * @param node the node whose attribute is being fetched.
     */
    Iterator getNodeAttributeNames(NodeRef node);

    /**
     * @return whether the node is external.
     */
    boolean isExternal(NodeRef node);

    /**
     * @return whether the node is the root.
     */
    boolean isRoot(NodeRef node);

    /**
     * @return the number of children of the ith node.
     */
    int getChildCount(NodeRef node);

    /**
     * @return the jth child of the ith node
     */
    NodeRef getChild(NodeRef node, int j);

    NodeRef getParent(NodeRef node);

    /**
     * @return a clone of this tree
     */
    public Tree getCopy();

    public class MissingTaxonException extends Exception {
        /**
		 * 
		 */
		private static final long serialVersionUID = 8468656622238269963L;

		public MissingTaxonException(Taxon taxon) { super(taxon.getId()); }
    }

    /**
     * Static utility functions for trees.
     *
     */
    public class Utils {

        public static final int NO_BRANCH_LENGTHS = 0;
        public static final int LENGTHS_AS_TIME = 1;
        public static final int LENGTHS_AS_SUBSTITUTIONS = 2;

        /**
         * Return the number of leaves under this node.
         * @param tree
         * @param node
         * @return the number of leaves under this node.
         */
        public static final int getLeafCount(Tree tree, NodeRef node) {

            int childCount = tree.getChildCount(node);
            if (childCount == 0) return 1;

            int leafCount = 0;
            for (int i = 0; i < childCount; i++) {
                leafCount += getLeafCount(tree, tree.getChild(node, i));
            }
            return leafCount;
        }

        public static double getMinNodeHeight(Tree tree, NodeRef node) {

            int childCount = tree.getChildCount(node);
            if (childCount == 0) return tree.getNodeHeight(node);

            double minNodeHeight = Double.MAX_VALUE;
            for (int i = 0; i < childCount; i++) {
                double height = getMinNodeHeight(tree, tree.getChild(node, i));
                if (height < minNodeHeight) {
                    minNodeHeight = height;
                }
            }
            return minNodeHeight;
        }

        /**
         * @return true only if all tips have height 0.0
         */
        public static boolean isUltrametric(Tree tree) {
            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                if (tree.getNodeHeight(tree.getExternalNode(i)) != 0.0) return false;
            }
            return true;
        }

        /**
         * @return true only if internal nodes have 2 children
         */
        public static boolean isBinary(Tree tree) {
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                if (tree.getChildCount(tree.getInternalNode(i)) > 2) return false;
            }
            return true;
        }

        /**
         * @return a set of strings which are the taxa of the tree.
         */
        public static Set getLeafSet(Tree tree) {

            HashSet leafSet = new HashSet();
            int m = tree.getTaxonCount();

            for (int i = 0; i < m; i++) {

                Taxon taxon = tree.getTaxon(i);
                leafSet.add(taxon.getId());
            }

            return leafSet;
        }

        /**
         * Gets an array of booleans, true for the leaf nodes that are associated with
         * the taxa in taxa.
         */
        public static Set getLeavesForTaxa(Tree tree, TaxonList taxa) throws MissingTaxonException {

            HashSet leafNodes = new HashSet();
            int m = taxa.getTaxonCount();
            int n = tree.getExternalNodeCount();

            for (int i = 0; i < m; i++) {

                Taxon taxon = taxa.getTaxon(i);
                boolean found = false;
                for (int j = 0; j < n; j++) {

                    NodeRef node = tree.getExternalNode(j);
                    if (tree.getNodeTaxon(node).getId().equals(taxon.getId())) {

                        found = true;
                        break;
                    }
                }

                if (!found) {
                    throw new MissingTaxonException(taxon);
                }

                leafNodes.add(taxon.getId());
            }

            return leafNodes;
        }

        /**
         * Gets a set of taxa names (as strings) of the leaf nodes descended from the given node.
         */
        public static Set getDescendantLeaves(Tree tree, NodeRef node) {

            HashSet set = new HashSet();
            getDescendantLeaves(tree, node, set);
            return set;
        }

        /**
         * Gets a set of taxa names (as strings) of the leaf nodes descended from the given node.
         */
        private static void getDescendantLeaves(Tree tree, NodeRef node, Set set) {

            if (tree.isExternal(node)) {
                set.add(tree.getTaxonId(node.getNumber()));
            } else {

                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    getDescendantLeaves(tree, node1, set);
                }

            }
        }

        /**
         * Gets the most recent common ancestor (MRCA) node of a set of leaf nodes.
         * @param tree the Tree
         * @param leafNodes a set of names
         * @return the NodeRef of the MRCA
         */
        public static NodeRef getCommonAncestorNode(Tree tree, Set leafNodes) {

            int cardinality = leafNodes.size();

            if (cardinality == 0) {
                throw new IllegalArgumentException("No leaf nodes selected");
            }

            NodeRef[] mrca = new NodeRef[] { null };
            getCommonAncestorNode(tree, tree.getRoot(), leafNodes, cardinality, mrca);

            return mrca[0];
        }

        /**
         * Private recursive function used by getCommonAncestorNode.
         */
        private static int getCommonAncestorNode(Tree tree, NodeRef node,
                                                 Set leafNodes, int cardinality,
                                                 NodeRef[] mrca) {

            if (tree.isExternal(node)) {

                if (leafNodes.contains(tree.getTaxonId(node.getNumber()))) {
                    if (cardinality == 1) {
                        mrca[0] = node;
                    }
                    return 1;
                } else {
                    return 0;
                }
            }

            int matches = 0;

            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                matches += getCommonAncestorNode(tree, node1, leafNodes, cardinality, mrca);

                if (mrca[0] != null) { break; }
            }

            if (mrca[0] == null) {
                // If we haven't already found the MRCA, test this node
                if (matches == cardinality) {
                    mrca[0] = node;
                }
            }

            return matches;
        }

        /**
         * Performs the a monophyly test on a set of leaf nodes. The nodes are monophyletic
         * if there is a node in the tree which subtends all the taxa in the set (and
         * only those taxa).
         * @param tree a tree object to perform test on
         * @param leafNodes a array with one boolean for each leaf node.
         * @return boolean is monophyletic?
         */
        public static boolean isMonophyletic(Tree tree, Set leafNodes) {

            int cardinality = leafNodes.size();

            if (cardinality == 1) {
                // A single selected leaf is always monophyletic
                return true;
            }

            if (cardinality == tree.getExternalNodeCount()) {
                // All leaf nodes are selected
                return true;
            }

            if (cardinality == 0) {
                throw new IllegalArgumentException("No leaf nodes selected");
            }

            int[] matchCount = new int[] { 0 };
            int[] leafCount = new int[] { 0 };
            boolean[] isMono = new boolean[] { false };
            isMonophyletic(tree, tree.getRoot(), leafNodes, cardinality, matchCount, leafCount, isMono);

            return isMono[0];
        }

        /**
         * Private recursive function used by isMonophyletic.
         */
        private static boolean isMonophyletic(Tree tree, NodeRef node,
                                              Set leafNodes, int cardinality,
                                              int[] matchCount, int[] leafCount,
                                              boolean[] isMono) {

            if (tree.isExternal(node)) {

                if (leafNodes.contains(tree.getNodeTaxon(node).getId())) {
                    matchCount[0] = 1;
                } else {
                    matchCount[0] = 0;
                }
                leafCount[0] = 1;
                return false;
            }

            int mc = 0;
            int lc = 0;

            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                boolean done = isMonophyletic(tree, node1, leafNodes, cardinality, matchCount, leafCount, isMono);
                mc += matchCount[0];
                lc += leafCount[0];

                if (done) { return true; }
            }

            matchCount[0] = mc;
            leafCount[0] = lc;

            // If we haven't already found the MRCA, test this node
            if (mc == lc && lc == cardinality) {
                isMono[0] = true;
                return true;
            }

            return false;
        }

        /**
         * @return the size of the largest clade with tips in the given range of times.
         */
        public static int largestClade(Tree tree, double range) {

            return largestClade(tree, tree.getRoot(), range, new double[] {0.0, 0.0});

        }

        /**
         * @return the size of the largest clade with tips in the given range of times.
         */
        private static int largestClade(Tree tree, NodeRef node, double range, double[] currentBounds) {

            if (tree.isExternal(node)) {
                currentBounds[0] = tree.getNodeHeight(node);
                currentBounds[1] = tree.getNodeHeight(node);
                return 1;
            } else {
                // get the bounds and max clade size of the left clade
                int cladeSize1 = largestClade(tree, tree.getChild(node, 0), range, currentBounds);
                double min = currentBounds[0];
                double max = currentBounds[1];

                // get the bounds and max clade size of the right clade
                int cladeSize2 = largestClade(tree, tree.getChild(node, 1), range, currentBounds);
                min = Math.min(min, currentBounds[0]);
                max = Math.max(max, currentBounds[1]);

                // update the joint bounds
                currentBounds[0] = min;
                currentBounds[1] = max;

                // if the joint clade is valid return the joint size
                if (max - min < range) {
                    return cladeSize1 + cladeSize2;
                }
                // if the joint clade is not valid return the max of the two
                return Math.max(cladeSize1, cladeSize2);
            }
        }

        /**
         * Calculates the minimum number of steps for the parsimony reconstruction of a
         * binary character defined by leafStates.
         * @param tree a tree object to perform test on
         * @param leafStates a set of booleans, one for each leaf node
         * @return number of parsimony steps
         */
        public static int getParsimonySteps(Tree tree, Set leafStates) {

            int[] score = new int[] { 0 };
            getParsimonySteps(tree, tree.getRoot(), leafStates, score);
            return score[0];
        }

        private static int getParsimonySteps(Tree tree, NodeRef node, Set leafStates, int[] score) {

            if (tree.isExternal(node)) {
                return (leafStates.contains(tree.getTaxonId(node.getNumber())) ? 1 : 2);

            } else {

                int uState = getParsimonySteps(tree, tree.getChild(node, 0), leafStates, score);
                int iState = uState;

                for (int i = 1; i < tree.getChildCount(node); i++) {

                    int state = getParsimonySteps(tree, tree.getChild(node, i), leafStates, score);
                    uState = state | uState;

                    iState = state & iState;

                }

                if (iState == 0) {
                    score[0] += 1;
                }

                return uState;
            }

        }

        /**
         * Calculates the parsimony reconstruction of a binary character defined
         * by leafStates at a given node.
         * @param tree a tree object to perform test on
         * @param node a NodeRef object from tree
         * @param leafStates a set of booleans, one for each leaf node
         * @return number of parsimony steps
         */
        public static double getParsimonyState(Tree tree, NodeRef node, Set leafStates) {

            int state = getParsimonyStateAtNode(tree, node, leafStates);
            switch (state) {
                case 1: return 0.0;
                case 2: return 1.0;
                default: return 0.5;
            }
        }

        private static int getParsimonyStateAtNode(Tree tree, NodeRef node, Set leafStates) {

            if (tree.isExternal(node)) {
                return (leafStates.contains(tree.getTaxonId(node.getNumber())) ? 1 : 2);

            } else {

                int uState = getParsimonyStateAtNode(tree, tree.getChild(node, 0), leafStates);
                int iState = uState;

                for (int i = 1; i < tree.getChildCount(node); i++) {

                    int state = getParsimonyStateAtNode(tree, tree.getChild(node, i), leafStates);
                    uState = state | uState;

                    iState = state & iState;

                }

                return uState;
            }

        }

        /**
         * determine preorder successor of this node
         *
         * @return next node
         */
        public static NodeRef preorderSuccessor(Tree tree, NodeRef node) {

            NodeRef next = null;

            if (tree.isExternal(node)) {
                NodeRef cn = node, ln = null; // Current and last node

                // Go up
                do {
                    if (tree.isRoot(cn)) {
                        next = cn;
                        break;
                    }
                    ln = cn;
                    cn = tree.getParent(cn);
                }
                while (tree.getChild(cn, tree.getChildCount(cn)-1) == ln);

                // Determine next node
                if (next == null) {
                    // Go down one node
                    for (int i = 0; i < tree.getChildCount(cn)-1; i++) {

                        if (tree.getChild(cn, i) == ln) {
                            next = tree.getChild(cn, i+1);
                            break;
                        }
                    }
                }
            } else {
                next = tree.getChild(node, 0);
            }

            return next;
        }

        /**
         * determine postorder successor of a node
         *
         * @return next node
         */
        public static NodeRef postorderSuccessor(Tree tree, NodeRef node) {

            NodeRef cn = null;
            NodeRef parent = tree.getParent(node);

            if (tree.getRoot() == node) {
                cn = node;
            } else {

                // Go up one node
                if (tree.getChild(parent, tree.getChildCount(parent)-1) == node) {
                    return parent;
                }

                // Go down one node
                for (int i = 0; i < tree.getChildCount(parent)-1; i++) {
                    if (tree.getChild(parent, i) == node) {
                        cn = tree.getChild(parent, i+1);
                        break;
                    }
                }
            }

            // Go down until leaf
            while (tree.getChildCount(cn) > 0) {
                cn = tree.getChild(cn, 0);
            }

            return cn;
        }

        /**
         * Gets finds the most ancestral node with attribute set.
         */
        public static NodeRef findNodeWithAttribute(Tree tree, String attribute) {

            NodeRef root = tree.getRoot();
            NodeRef node = root;

            do {

                if (tree.getNodeAttribute(node, attribute) != null) {
                    return node;
                }

                node = Tree.Utils.preorderSuccessor(tree, node);

            } while (node != root);

            return null;
        }

        /**
         * Gets finds the most recent date amongst the external nodes.
         */
        public static dr.evolution.util.Date findMostRecentDate(Tree tree) {

            dr.evolution.util.Date mostRecent = null;

            for (int i =0; i < tree.getExternalNodeCount(); i++) {
                Taxon taxon = tree.getNodeTaxon(tree.getExternalNode(i));

                dr.evolution.util.Date date = (dr.evolution.util.Date)taxon.getAttribute(dr.evolution.util.Date.DATE);
                if ((date != null) && (mostRecent == null || date.after(mostRecent))) {
                    mostRecent = date;
                }
            }

            return mostRecent;
        }

        /**
         * Recursive function for constructing a newick tree representation in the given buffer.
         */
        public static String newick(Tree tree) {
            StringBuffer buffer = new StringBuffer();
            newick(tree, tree.getRoot(), true, LENGTHS_AS_TIME, null, null, null, null, buffer);
            buffer.append(";");
            return buffer.toString();
        }

        public static String newick(Tree tree, BranchRateController branchRateController) {
            StringBuffer buffer = new StringBuffer();
            newick(tree, tree.getRoot(), true, LENGTHS_AS_SUBSTITUTIONS, branchRateController, null, null, null, buffer);
            buffer.append(";");
            return buffer.toString();
        }

        public static String newick(Tree tree,
                                    BranchRateController branchRateController, String rateLabel) {
            StringBuffer buffer = new StringBuffer();
            newick(tree, tree.getRoot(), true, LENGTHS_AS_TIME,
                    branchRateController, rateLabel,
                    null, null,
                    buffer);
            buffer.append(";");
            return buffer.toString();
        }

        public static String newick(Tree tree,
                                    TreeColouring colouring, String colouringLabel) {
            StringBuffer buffer = new StringBuffer();
            newick(tree, tree.getRoot(), true, LENGTHS_AS_TIME,
                    null, null,
                    colouring, colouringLabel,
                    buffer);
            buffer.append(";");
            return buffer.toString();
        }

        public static String newick(Tree tree,
                                    BranchRateController branchRateController, String rateLabel,
                                    TreeColouring colouring, String colouringLabel) {
            StringBuffer buffer = new StringBuffer();
            newick(tree, tree.getRoot(), true, LENGTHS_AS_TIME,
                    branchRateController, rateLabel,
                    colouring, colouringLabel,
                    buffer);
            buffer.append(";");
            return buffer.toString();
        }

        /**
         * Recursive function for constructing a newick tree representation in the given buffer.
         */
        public static String newickNoLengths(Tree tree) {
            StringBuffer buffer = new StringBuffer();
            newick(tree, tree.getRoot(), true, NO_BRANCH_LENGTHS, null, null, null, null, buffer);
            buffer.append(";");
            return buffer.toString();
        }

        /**
         * Recursive function for constructing a newick tree representation in the given buffer.
         * @param tree The tree
         * @param node The node [tree.getRoot()]
         * @param labels whether labels or numbers should be used
         * @param lengths What type of branch lengths: NO_BRANCH_LENGTHS, LENGTHS_AS_TIME, LENGTHS_AS_SUBSTITUTIONS
         * @param branchRateController An optional BranchRateController (or null)
         * @param rateLabel The label to use for rates in the tree (or null)
         * @param colouring An option TreeColouring object (or null)
         * @param colouringLabel The label for colourings in the tree (or null)
         * @param buffer The StringBuffer
         */
        public static void newick(Tree tree, NodeRef node, boolean labels, int lengths,
                                  BranchRateController branchRateController, String rateLabel,
                                  TreeColouring colouring, String colouringLabel,
                                  StringBuffer buffer) {

            NodeRef parent = tree.getParent(node);

            if (tree.isExternal(node)) {
                if (!labels) {
                    int k = node.getNumber() + 1;
                    buffer.append(k);
                } else {
                    buffer.append(tree.getTaxonId(node.getNumber()));
                }
            } else {
                buffer.append("(");
                newick(tree, tree.getChild(node, 0), labels, lengths,
                        branchRateController, rateLabel,
                        colouring, colouringLabel,
                        buffer);
                for (int i = 1; i < tree.getChildCount(node); i++) {
                    buffer.append(",");
                    newick(tree, tree.getChild(node, i), labels, lengths,
                            branchRateController, rateLabel,
                            colouring, colouringLabel,
                            buffer);
                }
                buffer.append(")");
            }

            boolean hasAttribute = false;
            if (colouring != null && colouringLabel != null) {
                BranchColouring bcol = colouring.getBranchColouring(node);
                if (bcol != null) {
                    buffer.append("[&").append(colouringLabel).append("={");
                    buffer.append(bcol.getChildColour());
                    for (int i = 1; i <= bcol.getNumEvents(); i++) {
                        buffer.append(",");
                        buffer.append(bcol.getBackwardTime(i));
                        buffer.append(",");
                        buffer.append(bcol.getBackwardColourAbove(i));
                    }
                    buffer.append("}");
                    hasAttribute = true;
                }
            }

            if (branchRateController != null && rateLabel != null && !tree.isRoot(node)) {
                double rate = branchRateController.getBranchRate(tree, node);
                if (!hasAttribute) {
                    buffer.append("[&");
                } else {
                    buffer.append(",");
                }
                buffer.append(rateLabel).append("=");
                buffer.append(String.valueOf(rate));
                hasAttribute = true;
            }
            if (hasAttribute) {
                buffer.append("]");
            }

            if (parent != null && lengths != NO_BRANCH_LENGTHS) {
                buffer.append(":");
                if (lengths == LENGTHS_AS_TIME) {
                    double length = tree.getNodeHeight(parent) - tree.getNodeHeight(node);
                    buffer.append(String.valueOf(length));
                } else if (lengths == LENGTHS_AS_SUBSTITUTIONS) {
                    if (branchRateController == null) {
                        throw new IllegalArgumentException("No BranchRateController provided");
                    }
                    double length = tree.getNodeHeight(parent) - tree.getNodeHeight(node);
                    double rate = branchRateController.getBranchRate(tree, node);
                    buffer.append(String.valueOf(length * rate));
                }
            }
        }

        /**
         * Recursive function for constructing a newick tree representation in the given buffer.
         */
        public static String uniqueNewick(Tree tree, NodeRef node) {
            if (tree.isExternal(node)) {
                //buffer.append(tree.getNodeTaxon(node).getId());
                return tree.getNodeTaxon(node).getId();
            } else {
                StringBuffer buffer = new StringBuffer("(");

                ArrayList subtrees = new ArrayList();
                for (int i = 0; i < tree.getChildCount(node); i++) {
                    NodeRef child = tree.getChild(node, i);
                    subtrees.add(uniqueNewick(tree, child));
                }
                Collections.sort(subtrees);
                for (int i = 0; i < subtrees.size(); i++) {
                    buffer.append(subtrees.get(i));
                    if (i < subtrees.size() - 1) {
                        buffer.append(",");
                    }
                }
                buffer.append(")");

                return buffer.toString();
            }
        }

        /**
         * Recursive function for constructing a newick tree representation in the given buffer.
         */
        public static Tree rotateByName(Tree tree) {

            return new SimpleTree(rotateNodeByName(tree, tree.getRoot()));
        }

        /**
         * Recursive function for constructing a newick tree representation in the given buffer.
         */
        private static SimpleNode rotateNodeByName(Tree tree, NodeRef node) {

            if (tree.isExternal(node)) {
                return  new SimpleNode(tree, node);
            } else {

                SimpleNode parent = new SimpleNode(tree, node);

                NodeRef child1 = tree.getChild(node, 0);
                NodeRef child2 = tree.getChild(node, 1);

                String subtree1 = uniqueNewick(tree, child1);
                String subtree2 = uniqueNewick(tree, child2);

                if (subtree1.compareTo(subtree2) > 0) {
                    parent.addChild(rotateNodeByName(tree, child2));
                    parent.addChild(rotateNodeByName(tree, child1));
                } else {
                    parent.addChild(rotateNodeByName(tree, child1));
                    parent.addChild(rotateNodeByName(tree, child2));
                }
                return parent;
            }
        }

        public static MutableTree rotateTreeByComparator(Tree tree, Comparator comparator) {

            return new SimpleTree(rotateTreeByComparator(tree, tree.getRoot(), comparator));
        }

        /**
         * Recursive function for constructing a newick tree representation in the given buffer.
         */
        private static SimpleNode rotateTreeByComparator(Tree tree, NodeRef node, Comparator comparator) {

            SimpleNode newNode = new SimpleNode();
            newNode.setHeight(tree.getNodeHeight(node));
            newNode.setRate(tree.getNodeRate(node));
            newNode.setId(tree.getTaxonId(node.getNumber()));
            newNode.setNumber(node.getNumber());
            newNode.setTaxon(tree.getNodeTaxon(node));

            if (!tree.isExternal(node)) {

                NodeRef child1 = tree.getChild(node, 0);
                NodeRef child2 = tree.getChild(node, 1);

                if (comparator.compare(child1, child2) > 0) {
                    newNode.addChild(rotateTreeByComparator(tree, child2, comparator));
                    newNode.addChild(rotateTreeByComparator(tree, child1, comparator));
                } else {
                    newNode.addChild(rotateTreeByComparator(tree, child1, comparator));
                    newNode.addChild(rotateTreeByComparator(tree, child2, comparator));
                }
            }

            return newNode;
        }

        public static Comparator createNodeDensityComparator(final Tree tree) {

            return new Comparator() {

                public int compare(Object o1, Object o2) {
                    NodeRef node1 = (NodeRef)o1;
                    NodeRef node2 = (NodeRef)o2;

                    return getLeafCount(tree, node2) - getLeafCount(tree, node1);
                }

                public boolean equals(Object o1, Object o2) {
                    NodeRef node1 = (NodeRef)o1;
                    NodeRef node2 = (NodeRef)o2;

                    return getLeafCount(tree, node1) == getLeafCount(tree, node2);
                }
            };
        }

        public static Comparator createNodeDensityMinNodeHeightComparator(final Tree tree) {

            return new Comparator() {

                public int compare(Object o1, Object o2) {
                    NodeRef node1 = (NodeRef)o1;
                    NodeRef node2 = (NodeRef)o2;

                    int larger = getLeafCount(tree, node1) - getLeafCount(tree, node2);

                    if (larger != 0) return larger;

                    double tipRecent = getMinNodeHeight(tree, node2) - getMinNodeHeight(tree, node1);
                    if (tipRecent > 0.0) return 1;
                    if (tipRecent < 0.0) return -1;
                    return 0;
                }
            };
        }

        /**
         * Compares 2 trees and returns true if they have the same topology (same taxon
         * order is assumed).
         */
        public static boolean equal(Tree tree1, Tree tree2) {

            return uniqueNewick(tree1, tree1.getRoot()).equals(uniqueNewick(tree2, tree2.getRoot()));
        }



    }

}