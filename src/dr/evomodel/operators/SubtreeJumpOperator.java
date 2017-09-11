/*
 * SubtreeJumpOperator.java
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
 * Implements the Subtree Jump move.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SubtreeJumpOperator extends AbstractTreeOperator implements CoercableMCMCOperator {

    private double size = 1.0;
    private double accP = 0.234;
    
    private final TreeModel tree;
    private final CoercionMode mode;
    private final boolean uniform;

    /**
     * Constructor
     * @param tree
     * @param weight
     * @param size: the variance of a half normal used to compute distance weights (as a rule, larger size, bolder moves)
     * @param mode
     */
    public SubtreeJumpOperator(TreeModel tree, double weight, double size, double accP, CoercionMode mode) {
        this.tree = tree;
        setWeight(weight);
        this.size = size;
        this.accP = accP;
        this.mode = mode;

        uniform = (size == Double.POSITIVE_INFINITY);
    }
    /**
     * Do a subtree jump move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {
        double logq;

        final NodeRef root = tree.getRoot();

        double  maxHeight = tree.getNodeHeight(root);

        NodeRef i;
        NodeRef iP = null;
        NodeRef CiP = null;
        NodeRef PiP = null;
        double height = Double.NaN;
        List<NodeRef> destinations = null;

        do {
            // 1. choose a random node avoiding root or child of root
            i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));

        } while (root == i || tree.getParent(i) == root);

        iP = tree.getParent(i);
        CiP = getOtherChild(tree, iP, i);
        PiP = tree.getParent(iP);

        // get the height of the parent
        height = tree.getNodeHeight(iP);

        // get a list of all edges that intersect this height
        destinations = getIntersectingEdges(tree, height);

        if (destinations.size() == 0) {
            // if there are no destinations available then reject the move
            return Double.NEGATIVE_INFINITY;
        }

        double forwardProbability;
        int r;

        if (uniform) {
            r = MathUtils.nextInt(destinations.size());
            forwardProbability = 1.0;
        } else {
            double[] pdf = getDestinationProbabilities(tree, i, height, maxHeight, destinations, size);

            // pick uniformly from this list
            r = MathUtils.randomChoicePDF(pdf);
            forwardProbability = pdf[r];
        }

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

        if (uniform) {
            logq = 0.0;
        } else {

            final List<NodeRef> reverseDestinations = getIntersectingEdges(tree, height);
            double reverseProbability = getReverseProbability(tree, CiP, j, height, maxHeight, reverseDestinations, size);

            // hastings ratio = reverse Prob / forward Prob
            logq = Math.log(reverseProbability) - Math.log(forwardProbability);
        }

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

            // The original node and its sibling will not be included because their height is exactly equal to the target height
            if (parent != null && tree.getNodeHeight(node) < height && tree.getNodeHeight(parent) > height) {
                intersectingEdges.add(node);
            }
        }
        return intersectingEdges;
    }

    private double[] getDestinationProbabilities(Tree tree, NodeRef node0, double height, double maxAge, List<NodeRef> intersectingEdges, double size) {
        double[] weights = new double[intersectingEdges.size()];

        double sum = 0.0;
        int i = 0;
        for (NodeRef node1 : intersectingEdges) {
            assert (node1 != node0);

            double age = tree.getNodeHeight(TreeUtils.getCommonAncestor(tree, node0, node1)) - height;
            age = age / maxAge;
            weights[i] = getJumpWeight(age, size, uniform);

            sum += weights[i];
            i++;
        }
        for (int j = 0; j < weights.length; j++) {
            weights[j] /= sum;
        }

        return weights;
    }

    private double getReverseProbability(Tree tree, NodeRef originalNode, NodeRef targetNode, double height, double maxAge, List<NodeRef> intersectingEdges, double size) {

        double[] weights = new double[intersectingEdges.size()];
        double sum = 0.0;

        int i = 0;
        int originalIndex = -1;
        for (NodeRef node1 : intersectingEdges) {
            assert(node1 != targetNode);

            double age = tree.getNodeHeight(TreeUtils.getCommonAncestor(tree, targetNode, node1)) - height;
            age = age/maxAge;
            weights[i] = getJumpWeight(age, size, uniform);
            sum += weights[i];

            if (node1 == originalNode) {
                originalIndex = i;
            }
            i++;
        }

        return weights[originalIndex] /= sum;

    }

    private double getJumpWeight(double age, double size, boolean uniform) {
    	if(uniform){
    		return 1;
    	}else{
    		return dr.math.distributions.NormalDistribution.pdf(age, 0, size) + Double.MIN_VALUE;
    	}
    }

    public double getCoercableParameter() {
        return size;
    }

    public void setCoercableParameter(double value) {
        size = value;
    }

    public double getRawParameter() {
        return size;
    }

    public CoercionMode getMode() {
        return mode;
    }

    public double getTargetAcceptanceProbability() {
        return accP;
    }
    public String getPerformanceSuggestion() {
        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();

        if (size <=0) {
            return "";
        }

        double ws = OperatorUtils.optimizeWindowSize(size, Double.MAX_VALUE, prob, targetProb);

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
