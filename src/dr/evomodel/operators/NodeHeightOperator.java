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
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.SubtreeLeapOperatorParser;
import dr.evomodelxml.operators.TipLeapOperatorParser;
import dr.inference.distribution.CauchyDistribution;
import dr.inference.model.Bounds;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.RandomWalkOperator;
import dr.math.MathUtils;
import dr.math.distributions.Distribution;
import dr.util.Transform;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Implements moves that changes node heights.
 *
 * Created to allow these operations on TreeModels that can't expose node heights as parameters for efficiency
 * reasons such as BigFastTreeModel.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class NodeHeightOperator extends AbstractAdaptableTreeOperator {

    public enum OperatorType {
        UNIFORM("uniform"),
        RANDOM_WALK("randomWalk");

        OperatorType(String name) {
            this.name = name;
        }

        @Override
        public String toString() {
            return name;
        }

        String name;
    }

    private double size;

    private final TreeModel tree;
    private final OperatorType operatorType;

    private final List<NodeRef> tips;

    /**
     * Constructor
     *
     * @param tree   the tree
     * @param weight the weight
     * @param size   size of move for types that need it
     * @param targetAcceptance the desired acceptance probability
     * @param operatorType the type of move to make
     * @param mode   coercion mode
     */
    public NodeHeightOperator(TreeModel tree, double weight, double size, OperatorType operatorType, AdaptationMode mode, double targetAcceptance) {
        super(mode, targetAcceptance);

        this.tree = tree;
        setWeight(weight);
        this.size = size;
        this.operatorType = operatorType;
        this.tips = null;
    }

    /**
     * Constructor that takes a taxon set to pick from for the move.
     *
     * @param tree   the tree
     * @param taxa   some taxa
     * @param weight the weight
     * @param size   size of move for types that need it
     * @param mode   coercion mode
     */
    public NodeHeightOperator(TreeModel tree, TaxonList taxa, double weight, double size, OperatorType operatorType, AdaptationMode mode, double targetAcceptance) {
        super(mode, targetAcceptance);

        this.tree = tree;
        setWeight(weight);
        this.size = size;
        this.operatorType = operatorType;
        this.tips = new ArrayList<NodeRef>();

        for (Taxon taxon : taxa) {
            boolean found = false;
            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                NodeRef tip = tree.getExternalNode(i);
                if (tree.getNodeTaxon(tip).equals(taxon)) {
                    tips.add(tip);
                    found = true;
                    break;
                }
            }

            if (!found) {
                throw new IllegalArgumentException("Taxon, " + taxon.getId() + ", not found in tree with id " + tree.getId());
            }
        }
    }


    /**
     * Do a subtree leap move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {
        final NodeRef root = tree.getRoot();

        NodeRef node;

        if (tips == null) {
            // Pick a node (but not the root)
            do {
                // choose a random node avoiding root
                node = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));

            } while (node == root);
        } else {
            throw new UnsupportedOperationException("tip height sampling not supported yet");
            // Pick a tip from the specified set of tips.
//            node = tips.get(MathUtils.nextInt(tips.size()));
        }

        NodeRef parent = tree.getParent(node);

        final double upperHeight = tree.getNodeHeight(tree.getParent(parent));
        final double lowerHeight = Math.max(
                tree.getNodeHeight(tree.getChild(parent, 0)),
                tree.getNodeHeight(tree.getChild(parent, 1)));

        final double oldHeight = tree.getNodeHeight(node);

        double logq;

        switch (operatorType) {
            case UNIFORM:
                logq = doUniform(node, oldHeight, upperHeight, lowerHeight);
                break;
            case RANDOM_WALK:
                RandomWalkOperator.BoundaryCondition boundaryCondition = RandomWalkOperator.BoundaryCondition.rejecting;
                logq = doRandomWalk(node, oldHeight, upperHeight, lowerHeight, getSize(), boundaryCondition);
                break;
            default:
                throw new IllegalArgumentException("Unknown operatorType");
        }
        return logq;
    }

    private double doUniform(NodeRef node, double oldValue, double upper, double lower) {
        tree.setNodeHeight(node, (MathUtils.nextDouble() * (upper - lower)) + lower);
        return 0.0;
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

    public double reflectValueLoop(double value, double lower, double upper) {
        double newValue = value;

        while (newValue < lower || newValue > upper) {
            if (newValue < lower) {
                newValue = lower + (lower - newValue);
            }
            if (newValue > upper) {
                newValue = upper - (newValue - upper);

            }
        }

        return newValue;
    }

    public double getSize() {
        return size;
    }

    public void setSize(double size) {
        this.size = size;
    }

    @Override
    protected void setAdaptableParameterValue(double value) {
        setSize(Math.exp(value));
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(getSize());
    }

    @Override
    public double getRawParameter() {
        return getSize();
    }

    public String getAdaptableParameterName() {
        return "size";
    }

    public String getOperatorName() {
        if (tips == null) {
            return SubtreeLeapOperatorParser.SUBTREE_LEAP + "(" + tree.getId() + ")";
        } else {
            return TipLeapOperatorParser.TIP_LEAP + "(" + tree.getId() + ")";
        }
    }
}
