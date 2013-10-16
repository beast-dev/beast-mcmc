/*
 * RateVarianceScaleOperator.java
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

package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.RateVarianceScaleOperatorParser;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.OperatorUtils;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * A special operator for scaling the variance of the autocorrelated clock model
 * and subsequnetly the rates in a tree.
 *
 * @author Michael Defoin Platel
 */
public class RateVarianceScaleOperator extends AbstractCoercableOperator {

    private TreeModel tree;
    private Parameter variance;

    public RateVarianceScaleOperator(TreeModel tree, Parameter variance, double scale, CoercionMode mode) {
        super(mode);

        this.scaleFactor = scale;
        this.tree = tree;
        this.variance = variance;
    }


    /**
     * scale the rates of a subtree and return the hastings ratio.
     */
    public final double doOperation() throws OperatorFailedException {

        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

        //Scale the variance
        double oldValue = variance.getParameterValue(0);
        double newValue = scale * oldValue;
        double logq = -Math.log(scale);

        final Bounds<Double> bounds = variance.getBounds();
        if (newValue < bounds.getLowerLimit(0) || newValue > bounds.getUpperLimit(0)) {
            throw new OperatorFailedException("proposed value outside boundaries");
        }
        variance.setParameterValue(0, newValue);

        //Scale the rates of the tree accordingly
        NodeRef root = tree.getRoot();
        final int index = root.getNumber();

        List<NodeRef> listNode = new ArrayList<NodeRef>();
        getSubtree(listNode, tree.getNode(index));

        final double rateScale = Math.sqrt(scale);

        for (NodeRef node : listNode) {

            oldValue = tree.getNodeRate(node);
            newValue = oldValue * rateScale;

            tree.setNodeRate(node, newValue);
        }

        //  According to the hastings ratio in the scale Operator
        logq += (listNode.size() - 2) * Math.log(rateScale);

        return logq;
    }

    void getSubtree(List<NodeRef> listNode, NodeRef parent) {

        listNode.add(parent);
        int nbChildren = tree.getChildCount(parent);
        for (int c = 0; c < nbChildren; c++) {
            getSubtree(listNode, tree.getChild(parent, c));
        }
    }

    /**
     * This method should be overridden by operators that need to do something just before the return of doOperation.
     *
     * @param newValue the proposed parameter value
     * @param oldValue the old parameter value
     */
    void cleanupOperation(double newValue, double oldValue) {
        // DO NOTHING
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "rateVarianceScale";
    }

    public double getCoercableParameter() {
        return Math.log(1.0 / scaleFactor - 1.0);
    }

    public void setCoercableParameter(double value) {
        scaleFactor = 1.0 / (Math.exp(value) + 1.0);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getScaleFactor() {
        return scaleFactor;
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

        double prob = Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = OperatorUtils.optimizeScaleFactor(scaleFactor, prob, targetProb);
        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

    public String toString() {
        return RateVarianceScaleOperatorParser.SCALE_OPERATOR + "(" + " [" + scaleFactor + ", " + (1.0 / scaleFactor) + "]";
    }

    //PRIVATE STUFF

    private double scaleFactor = 0.5;
}
