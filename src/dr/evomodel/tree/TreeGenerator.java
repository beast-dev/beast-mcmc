/*
 * CoalescentSimulator.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.tree;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.tree.*;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.DemographicModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.math.UnivariateFunction;

import java.util.*;

/**
 * A utility class for generating trees compatible with various constraints.
 *
 * @author Andrew Rambaut
 * @version $Id: CoalescentSimulator.java,v 1.43 2005/10/27 10:40:48 rambaut Exp $
 */
public class TreeGenerator {

    public TreeGenerator() {
    }

    /**
     * Generates a tree from a set of subtrees.
     *
     * @param guideTree               a tree topology to act as a guide. Subtrees must be compatible with this.
     * @param subtrees                an array of tree to be used as subtrees
     * @param rootHeight              an optional root height with which to scale the whole tree
     * @return a generated tree
     */
    public SimpleTree simulateTree(Tree guideTree, Tree[] subtrees, double rootHeight) throws GenerationFailedException {
        SimpleTree tree;

        Map<String, SimpleNode> leafNodes = new HashMap<String, SimpleNode>();

        for (int i = 0; i < subtrees.length; i++) {
            Set<String> leafSet = Tree.Utils.getLeafSet(subtrees[i]);

            NodeRef mrca = Tree.Utils.getCommonAncestorNode(guideTree, leafSet);

            if (mrca == null) {
                throw new GenerationFailedException("Subtree " + (i+1) + " is not compatible with the guide tree");
            }

            for (int j = 0; j < subtrees[i].getExternalNodeCount(); j++) {
                NodeRef leaf = subtrees[i].getExternalNode(j);
                leafNodes.put(subtrees[i].getNodeTaxon(leaf).getId(), new SimpleNode(subtrees[i], leaf));
            }
        }

        tree = new SimpleTree(constructTree(guideTree, guideTree.getRoot(), leafNodes));

        if (tree.getNodeHeight(tree.getRoot()) >= rootHeight) {

        }

        if (!Double.isNaN(rootHeight) && rootHeight > 0.0) {
            scaleTree(tree, rootHeight);
        }

        return tree;
    }


    /**
     * @return the root node of a tree corresponding to the given subtree but with heights rescaled
     */
    private SimpleNode constructTree(Tree guideTree, NodeRef node, Map<String, SimpleNode> leafNodes) throws GenerationFailedException {
        SimpleNode newNode;

        if (guideTree.isExternal(node)) {
            SimpleNode leaf = leafNodes.get(guideTree.getNodeTaxon(node).getId());

            if (leaf == null) {
                throw new GenerationFailedException("Leaf taxon, " + guideTree.getNodeTaxon(node).getId() + ", not found in subtree.");
            }

            newNode = leaf;
        } else {
            newNode = new SimpleNode();

            double maxHeight = Double.NEGATIVE_INFINITY;

            for (int i = 0; i < guideTree.getInternalNodeCount(); i++) {
                NodeRef child = guideTree.getInternalNode(i);

                SimpleNode newChild = constructTree(guideTree, child, leafNodes);

                if (newChild.getHeight() > maxHeight) {
                    maxHeight = newChild.getHeight();
                }

                newNode.addChild(newChild);
            }

            newNode.setHeight(maxHeight);
        }

        return newNode;
    }

    private void scaleTree(MutableTree tree, double rootHeight) {
        // avoid empty tree
        if (tree.getRoot() == null) return;

        double scale = rootHeight / tree.getNodeHeight(tree.getRoot());
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef n = tree.getInternalNode(i);
            tree.setNodeHeight(n, tree.getNodeHeight(n) * scale);
        }
        MutableTree.Utils.correctHeightsForTips(tree);
    }

    private void scaleTree(MutableTree tree, NodeRef node, double rootHeight) {
        if (!tree.isExternal(node)) {

//            double maxHeight = Double.NEGATIVE_INFINITY;
//
//            for (int i = 0; i < guideTree.getInternalNodeCount(); i++) {
//                NodeRef child = guideTree.getInternalNode(i);
//
//                SimpleNode newChild = scaleTree(guideTree, child, leafNodes);
//
//                if (newChild.getHeight() > maxHeight) {
//                    maxHeight = newChild.getHeight();
//                }
//
//                newNode.addChild(newChild);
//            }
//
//            newNode.setHeight(maxHeight);
        }
    }


    public class GenerationFailedException extends Exception {
        public GenerationFailedException(String message) {
            super(message);
        }
    }
}
