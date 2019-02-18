/*
 * RateVarianceScaleOperator.java
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
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.RateVarianceScaleOperatorParser;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractAdaptableOperator;
import dr.inference.operators.AdaptationMode;
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
public class RateVarianceScaleOperator extends AbstractAdaptableOperator {

    private TreeModel tree;
    private Parameter variance;

    public RateVarianceScaleOperator(TreeModel tree, Parameter variance, double scale, AdaptationMode mode) {
        super(mode);

        this.scaleFactor = scale;
        this.tree = tree;
        this.variance = variance;
    }


    /**
     * scale the rates of a subtree and return the hastings ratio.
     */
    public final double doOperation() {

        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

        //Scale the variance
        double oldValue = variance.getParameterValue(0);
        double newValue = scale * oldValue;
        double logq = -Math.log(scale);

        final Bounds<Double> bounds = variance.getBounds();
        if (newValue < bounds.getLowerLimit(0) || newValue > bounds.getUpperLimit(0)) {
//            throw new OperatorFailedException("proposed value outside boundaries");
            return Double.NEGATIVE_INFINITY;
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

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(1.0 / scaleFactor - 1.0);
    }

    public void setAdaptableParameterValue(double value) {
        scaleFactor = 1.0 / (Math.exp(value) + 1.0);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getScaleFactor() {
        return scaleFactor;
    }

    public String getAdaptableParameterName() {
        return "scaleFactor";
    }

    public String toString() {
        return RateVarianceScaleOperatorParser.SCALE_OPERATOR + "(" + " [" + scaleFactor + ", " + (1.0 / scaleFactor) + "]";
    }

    //PRIVATE STUFF

    private double scaleFactor = 0.5;
}
