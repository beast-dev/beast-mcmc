/*
 * FixedHeightSubtreePruneRegraftOperator.java
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

package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.SubtreeJumpOperatorParser;
import dr.inference.operators.*;
import dr.math.distributions.NormalDistribution;
import dr.math.MathUtils;
import java.util.ArrayList;
import java.util.List;

/**
 * Implements the fixed height subtree prune regraft move. Described by Sebastian Hoehna et al
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class FixedHeightSubtreePruneRegraftOperator extends AbstractTreeOperator {

    private final TreeModel tree;


    /**
     * Constructor
     * @param tree
     * @param weight
     */
    public FixedHeightSubtreePruneRegraftOperator(TreeModel tree, double weight) {
        this.tree = tree;
        setWeight(weight);

    }
    /**
     * Do a subtree jump move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {

        final NodeRef root = tree.getRoot();

//        double  maxHeight = tree.getNodeHeight(root);

        NodeRef i;
        NodeRef iP = null;
        NodeRef CiP = null;
        NodeRef PiP = null;
        List<NodeRef> destinations = null;

        do {
            // 1. choose a random node avoiding root or child of root
            i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));

        } while (root == i || tree.getParent(i) == root);

        iP = tree.getParent(i);
        CiP = getOtherChild(tree, iP, i);
        PiP = tree.getParent(iP);

        // get the height of the parent
        double parentHeight = tree.getNodeHeight(iP);

        // get a list of all edges that intersect this height
        destinations = getIntersectingEdges(tree, parentHeight);

        if (destinations.size() == 0) {
            // if there are no destinations available then reject the move
            return Double.NEGATIVE_INFINITY;
        }

        int r = MathUtils.nextInt(destinations.size());

        // remove the target node and its sibling (shouldn't be there because their parent's height is exactly equal to the target height).
        destinations.remove(i);
        destinations.remove(CiP);

        final NodeRef j = destinations.get(r);
        final NodeRef jP = tree.getParent(j);

        tree.beginTreeEdit();

        // remove the parent of i by connecting its sibling to its grandparent.
        tree.removeChild(iP, CiP);
        tree.removeChild(PiP, iP);
        tree.addChild(PiP, CiP);

        // remove destination edge j from its parent
        tree.removeChild(jP, j);

        // add destination edge to the parent of i
        tree.addChild(iP, j);

        // and add the parent of i as a child of the former parent of j
        tree.addChild(jP, iP);

        tree.endTreeEdit();

        return 0.0;
    }

    /**
     * Gets a list of edges that subtend the given height
     * @param tree
     * @param height
     * @return
     */
    private List<NodeRef> getIntersectingEdges(Tree tree, double height) {

        List<NodeRef> intersectingEdges = new ArrayList<NodeRef>();

        for (int i = 0; i < tree.getNodeCount(); i++) {
            final NodeRef node = tree.getNode(i);
            final NodeRef parent = tree.getParent(node);

            // The original node and its sibling will not be included because their height is exactly equal to the target height
            if (parent != null && tree.getNodeHeight(node) < height && tree.getNodeHeight(parent) > height) {
                intersectingEdges.add(node);
            }
        }
        return intersectingEdges;
    }

    public String getPerformanceSuggestion() {
        return null;
    }


    public String getOperatorName() {
        return SubtreeJumpOperatorParser.SUBTREE_JUMP + "(" + tree.getId() + ")";
    }
}
