/*
 * SubtreeSlideOperator.java
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

package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.SubtreeJumpOperatorParser;
import dr.evomodelxml.operators.SubtreeSlideOperatorParser;
import dr.inference.operators.*;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * Implements the Subtree Jump move.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SubtreeJumpOperator extends AbstractTreeOperator implements CoercableMCMCOperator {

    private TreeModel tree = null;
    private double size = 1;
    private CoercionMode mode = CoercionMode.DEFAULT;

    /**
     * Constructor
     * @param tree
     * @param weight
     * @param size a non-negative value used as a power coeficient for the relatedness weights
     * @param mode
     */
    public SubtreeJumpOperator(TreeModel tree, double weight, double size, CoercionMode mode) {
        this.tree = tree;
        setWeight(weight);
        this.size = size;
        this.mode = mode;
    }

    /**
     * Do a subtree jump move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() throws OperatorFailedException {

        double logq;

        final NodeRef root = tree.getRoot();
        final double oldTreeHeight = tree.getNodeHeight(root);

        NodeRef i;

        // 1. choose a random node avoiding root or child of root
        do {
            i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
        } while (root == i || tree.getParent(i) == root);


        final NodeRef iP = tree.getParent(i);
        final NodeRef CiP = getOtherChild(tree, iP, i);
        final NodeRef PiP = tree.getParent(iP);

        // get the height of the parent
        final double height = tree.getNodeHeight(iP);

        // get a list of all edges that intersect this height
        final List<NodeRef> destinations = getIntersectingEdges(tree, height);

        if (destinations.size() < 1) {
            throw new OperatorFailedException("no destinations to jump to");
        }

        double[] pdf = getDestinationProbabilities(tree, i, height, destinations);

        // remove the target node and its sibling (shouldn't be there because their parent's height is exactly equal to the target height).
        destinations.remove(i);
        destinations.remove(CiP);

        // pick uniformly from this list
        int r = MathUtils.randomChoicePDF(pdf);

        double forwardProbability = pdf[r];
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

        double reverseProbability = getReverseProbability(tree, j, i, height, destinations);

        logq = Math.log(forwardProbability) - Math.log(reverseProbability);

        return logq;
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
            if (parent != null && tree.getNodeHeight(node) < height && tree.getNodeHeight(parent) > height) {
                intersectingEdges.add(node);
            }
        }
        return intersectingEdges;
    }

    private double[] getDestinationProbabilities(Tree tree, NodeRef node0, double height, List<NodeRef> intersectingEdges) {
        double[] weights = new double[intersectingEdges.size()];
        double sum = 0.0;
        int i = 0;
        for (NodeRef node1 : intersectingEdges) {
            if (node1 != node0) {
                double age = tree.getNodeHeight(Tree.Utils.getCommonAncestor(tree, node0, node1)) - height;
                weights[i] = 1.0 / Math.pow(age, size);
                sum += weights[i];
            } else {
                weights[i] = 0.0;
            }
            i++;
        }
        for (int j = 0; j < weights.length; j++) {
            weights[j] /= sum;
        }

        return weights;
    }

    private double getReverseProbability(Tree tree, NodeRef originalNode, NodeRef targetNode, double height, List<NodeRef> intersectingEdges) {
        double[] weights = new double[intersectingEdges.size()];
        double sum = 0.0;
        int i = 0;
        int originalIndex = -1;
        for (NodeRef node1 : intersectingEdges) {
            if (node1 != targetNode) {
                double age = tree.getNodeHeight(Tree.Utils.getCommonAncestor(tree, targetNode, node1)) - height;
                weights[i] = 1.0 / Math.pow(age, size);
                sum += weights[i];

                if (node1 == originalNode) {
                    originalIndex = i;
                }
            } else {
                weights[i] = 0.0;
            }
            i++;
        }

        return weights[originalIndex] /= sum;
    }


    public double getSize() {
        return size;
    }

    public void setSize(double size) {
        this.size = size;
    }

    public double getCoercableParameter() {
        return Math.log(getSize());
    }

    public void setCoercableParameter(double value) {
        setSize(Math.exp(value));
    }

    public double getRawParameter() {
        return getSize();
    }

    public CoercionMode getMode() {
        return mode;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }


    public String getPerformanceSuggestion() {
        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();

        double ws = OperatorUtils.optimizeWindowSize(getSize(), Double.MAX_VALUE, prob, targetProb);

        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try decreasing size to about " + ws;
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try increasing size to about " + ws;
        } else return "";
    }


    public String getOperatorName() {
        return SubtreeJumpOperatorParser.SUBTREE_JUMP + "(" + tree.getId() + ")";
    }

}
