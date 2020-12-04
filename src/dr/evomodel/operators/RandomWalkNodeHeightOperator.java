/*
 * SubtreeLeapOperator.java
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
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.RandomWalkOperator;
import dr.math.MathUtils;
import dr.util.Transform;

/**
 * Implements moves that changes node heights - Random Walk version
 *
 * Created to allow these operations on TreeModels that can't expose node heights as parameters for efficiency
 * reasons such as BigFastTreeModel.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class RandomWalkNodeHeightOperator extends AbstractAdaptableTreeOperator {

    private double tunableParameter;

    private final TreeModel tree;

    /**
     * Constructor
     *
     * @param tree   the tree
     * @param weight the weight
     * @param tunableParameter   size of move for types that need it
     * @param targetAcceptance the desired acceptance probability
     * @param mode   coercion mode
     */
    public RandomWalkNodeHeightOperator(TreeModel tree, double weight, double tunableParameter, AdaptationMode mode, double targetAcceptance) {
        super(mode, targetAcceptance);

        this.tree = tree;
        setWeight(weight);
        this.tunableParameter = tunableParameter;
    }

    /**
     * Do a subtree leap move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {
        final NodeRef root = tree.getRoot();

        NodeRef node;
        // Pick a node (but not the root)
        do {
            // choose a internal node avoiding root
            node = tree.getInternalNode(MathUtils.nextInt(tree.getInternalNodeCount()));
        } while (node == root);

        final double upperHeight = tree.getNodeHeight(tree.getParent(node));
        final double lowerHeight = Math.max(
                tree.getNodeHeight(tree.getChild(node, 0)),
                tree.getNodeHeight(tree.getChild(node, 1)));

        final double oldHeight = tree.getNodeHeight(node);

        RandomWalkOperator.BoundaryCondition boundaryCondition = RandomWalkOperator.BoundaryCondition.rejecting;
        return doRandomWalk(node, oldHeight, upperHeight, lowerHeight, tunableParameter, boundaryCondition);
    }

    private double doRandomWalk(NodeRef node, double oldValue, double upper, double lower,
                                double windowSize, RandomWalkOperator.BoundaryCondition boundaryCondition) {

        // a random point around old value within windowSize * 2
        double draw = (2.0 * MathUtils.nextDouble() - 1.0) * windowSize;

        double newValue;
        double logHR;

        if (boundaryCondition == RandomWalkOperator.BoundaryCondition.logit) {
            // scale oldValue to [0,1]
            double x1 = (oldValue - lower) / (upper - lower);
            // logit transform it, add the draw, inverse transform it
            double x2 = Transform.LOGIT.inverse(Transform.LOGIT.transform(x1) + draw);

            // parameter takes new value scaled back into interval [lower, upper]
            newValue =  (x2 * (upper - lower)) + lower;

            // HR is the ratio of Jacobians for the before and after values in interval [0,1]
            logHR = Transform.LOGIT.getLogJacobian(x1) - Transform.LOGIT.getLogJacobian(x2);

        } else if (boundaryCondition == RandomWalkOperator.BoundaryCondition.log) {
            // offset oldValue to [0,+Inf]
            double x1 = oldValue - lower;
            // logit transform it, add the draw, inverse transform it
            double x2 = Transform.LOG.inverse(Transform.LOG.transform(x1) + draw);

            // parameter takes new value tranlated back into interval [lower, +Inf]
            newValue = x2 + lower;

            // HR is the ratio of Jacobians for the before and after values
            logHR = Transform.LOG.getLogJacobian(x1) - Transform.LOG.getLogJacobian(x2);

        } else {
            newValue = oldValue + draw;
            logHR = 0.0;
            if (boundaryCondition == RandomWalkOperator.BoundaryCondition.reflecting) {
                newValue = reflectValue(newValue, lower, upper);
            } else if (boundaryCondition == RandomWalkOperator.BoundaryCondition.rejecting && (newValue < lower || newValue > upper)) {
                logHR = Double.NEGATIVE_INFINITY;
            }
        }

        tree.setNodeHeight(node, newValue);

        return logHR;
    }

    public static double reflectValue(double value, double lower, double upper) {

        double newValue = value;

        if (upper == lower) {
            newValue = upper;
        } else if (value < lower) {
            if (Double.isInfinite(upper)) {
                // we are only going to reflect once as the upper bound is at infinity...
                newValue = lower + (lower - value);
            } else {
//                double remainder = lower - value;
//
//                double widths = Math.floor(remainder / (upper - lower));
//                remainder -= (upper - lower) * widths;

                final double ratio = (lower - value) / (upper - lower);
                final double widths = Math.floor(ratio);
                final double remainder = (ratio - widths) * (upper - lower);

                // even reflections
                if (widths % 2 == 0) {
                    newValue = lower + remainder;
                    // odd reflections
                } else {
                    newValue = upper - remainder;
                }
            }
        } else if (value > upper) {
            if (Double.isInfinite(lower)) {
                // we are only going to reflect once as the lower bound is at -infinity...
                newValue = upper - (newValue - upper);
            } else {

//                double remainder = value - upper;
//
//                double widths = Math.floor(remainder / (upper - lower));
//                remainder -= (upper - lower) * widths;

                final double ratio = (value - upper) / (upper - lower);
                final double widths = Math.floor(ratio);
                final double remainder = (ratio - widths) * (upper - lower);

                // even reflections
                if (widths % 2 == 0) {
                    newValue = upper - remainder;
                    // odd reflections
                } else {
                    newValue = lower + remainder;
                }
            }
        }

        return newValue;
    }

    @Override
    protected void setAdaptableParameterValue(double value) {
        tunableParameter = Math.exp(value);
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(tunableParameter);
    }

    @Override
    public double getRawParameter() {
        return tunableParameter;
    }

    public String getAdaptableParameterName() {
        return "windowSize";
    }

    public String getOperatorName() {
        return "randomWalk(" + tree.getId() + " internal nodes)";
    }
}
