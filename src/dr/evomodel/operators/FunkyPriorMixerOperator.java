/*
 * FunkyPriorMixerOperator.java
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
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.FunkyPriorMixerOperatorParser;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.inferencexml.operators.RandomWalkOperatorParser;
import dr.math.MathUtils;

/**
 * @author Marc A. Suchard
 * @author John J. Welch
 */
// Cleaning out untouched stuff. Can be resurrected if needed
@Deprecated
public class FunkyPriorMixerOperator extends
        SimpleMCMCOperator {
//        AbstractCoercableOperator {
//        implements GibbsOperator {

    public FunkyPriorMixerOperator(TreeModel treeModel, Parameter parameter, double windowSize, RandomWalkOperator.BoundaryCondition bc,
                                   double weight, CoercionMode mode) {
      
//        super(mode);
        this.treeModel = treeModel;
//        this.parameter = treeModel.getRootHeightParameter();
//        this.parameter = null;
        this.parameter = parameter;
        this.windowSize = windowSize;
        this.condition = bc;

        setWeight(weight);
    }

    public final double doOperation() {

//        // a random dimension to perturb
//        int index;
//        index = MathUtils.nextInt(parameter.getDimension());
//
//        double newValue = parameter.getParameterValue(index) + ((2.0 * MathUtils.nextDouble() - 1.0) * windowSize);
//
////        treeModel.setNodeHeight(node, 0.0);
////        treeModel.getNodeHeightLower(node);
//
//        double lower = parameter.getBounds().getLowerLimit(index);
//        double upper = parameter.getBounds().getUpperLimit(index);
//
//        while (newValue < lower || newValue > upper) {
//            if (newValue < lower) {
//                if (condition == RandomWalkOperator.BoundaryCondition.reflecting) {
//                    newValue = lower + (lower - newValue);
//                } else {
//                    throw new OperatorFailedException("proposed value outside boundaries");
//                }
//
//            }
//            if (newValue > upper) {
//                if (condition == RandomWalkOperator.BoundaryCondition.reflecting) {
//                    newValue = upper - (newValue - upper);
//                } else {
//                    throw new OperatorFailedException("proposed value outside boundaries");
//                }
//
//            }
//        }
//
//        parameter.setParameterValue(index, newValue);

//        double[] minNodeHeights = new double[treeModel.getNodeCount()];
//        recursivelyFindNodeMinHeights(treeModel, treeModel.getRoot(), minNodeHeights);
//
//        logForwardDensity = new Double(0.0);
//        logBackwardDensity = new Double(0.0);
//
//        try {
//
//        recursivelyDrawNodeHeights(treeModel, treeModel.getRoot(), 0.0, 0.0, minNodeHeights); //,
//                //logForwardDensity, logBackwardDensity);
//
//        } catch (Exception e) {
//            System.err.println("Got exception: " + e.getMessage());
//        }


        int iterations = 100;

        for (int i = 0; i < iterations; i++) {

//            NodeRef node = treeModel.getNode(MathUtils.nextInt(treeModel.getNodeCount()));
//            if (node != treeModel.getRoot()) {
//                treeModel.getN
//            }

            final int index = MathUtils.nextInt(parameter.getDimension());
            final Bounds<Double> bounds = parameter.getBounds();
            final double lower = bounds.getLowerLimit(index);
            final double upper = bounds.getUpperLimit(index);
            final double newValue = (MathUtils.nextDouble() * (upper - lower)) + lower;

//            parameter.setParameterValueQuietly(index, newValue);
            parameter.setParameterValue(index, newValue);
        }
//        ((Parameter.Default)parameter).fireParameterChangedEvent(-1, Parameter.ChangeType.VALUE_CHANGED);

//        System.err.println("logFD = " + logForwardDensity);
//        System.err.println("logBD = " + logBackwardDensity);

        return
                0.0;
//                -1 *
//        (logBackwardDensity - logForwardDensity);
    }


    private double recursivelyFindNodeMinHeights(Tree tree, NodeRef node, double[] minNodeHeights) {

        // Post-order traversal
        
        double minHeight;

        if (tree.isExternal(node))
            minHeight = tree.getNodeHeight(node);
        else {
            double minHeightChild0 = recursivelyFindNodeMinHeights(tree, tree.getChild(node, 0), minNodeHeights);
            double minHeightChild1 = recursivelyFindNodeMinHeights(tree, tree.getChild(node, 1), minNodeHeights);
            minHeight = (minHeightChild0 > minHeightChild1) ? minHeightChild0 : minHeightChild1;
        }

        minNodeHeights[node.getNumber()] = minHeight;
        return minHeight;
    }

    private void recursivelyDrawNodeHeights(TreeModel tree, NodeRef node,
                                            double oldParentHeight, double newParentHeight, double[] minNodeHeights) {//,
//                                            Double logForwardDensity, Double logBackwardDensity) {
        // Pre-order traversal

        if (tree.isExternal(node))
            return;

        final double oldNodeHeight = tree.getNodeHeight(node);
        double newNodeHeight = oldNodeHeight;

//        System.err.println("old: " + oldNodeHeight);
//
        if (!tree.isRoot(node)) {

            double minHeight = minNodeHeights[node.getNumber()];

            double oldDiff = oldParentHeight - minHeight;
            double newDiff = newParentHeight - minHeight;

            newNodeHeight = MathUtils.nextDouble() * newDiff + minHeight; // Currently uniform
            logForwardDensity -= Math.log(newDiff);
            logBackwardDensity -= Math.log(oldDiff);

//            System.err.println("inner logFD = " + logForwardDensity);

            tree.setNodeHeight(node, newNodeHeight);            
        }

//        System.err.println("new: " + newNodeHeight + "\n");
//
        recursivelyDrawNodeHeights(tree, tree.getChild(node, 0), oldNodeHeight, newNodeHeight, minNodeHeights); //,
//                logForwardDensity, logBackwardDensity);
        recursivelyDrawNodeHeights(tree, tree.getChild(node, 1), oldNodeHeight, newNodeHeight, minNodeHeights); //,
//                logForwardDensity, logBackwardDensity);

    }

    //MCMCOperator INTERFACE

    public final String getOperatorName() {
        return FunkyPriorMixerOperatorParser.FUNKY_OPERATOR;
    }

    public double getCoercableParameter() {
        return Math.log(windowSize);
    }

    public void setCoercableParameter(double value) {
        windowSize = Math.exp(value);
    }

    public double getRawParameter() {
        return windowSize;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public final String getPerformanceSuggestion() {

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();

//        double ws = OperatorUtils.optimizeWindowSize(windowSize, parameter.getParameterValue(0) * 2.0, prob, targetProb);
        double ws = 2;

        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try decreasing windowSize to about " + ws;
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try increasing windowSize to about " + ws;
        } else return "";
    }

    public String toString() {
        return RandomWalkOperatorParser.RANDOM_WALK_OPERATOR + "(" + parameter.getParameterName() + ", " + windowSize + ", " + getWeight() + ")";
    }

    public int getStepCount() {
        return 0;
    }

    private final TreeModel treeModel;
    private final Parameter parameter;
    private double windowSize;
    private final RandomWalkOperator.BoundaryCondition condition;

    private Double logForwardDensity;
    private Double logBackwardDensity;

    /**
     * @return the number of steps the operator performs in one go.
     */
}
