/*
 * UniformTipFromPrecisionOperator.java
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
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.operators.SimpleMCMCOperator;
import dr.evomodelxml.operators.UniformTipFromPrecisionOperatorParser;
import dr.math.MathUtils;

/**
 * Like UniformOperator but for tip times being sampled from precision (and are thus subject to tree-based age constraints)
 */
public class UniformTipFromPrecisionOperator extends SimpleMCMCOperator {

    public UniformTipFromPrecisionOperator(Parameter parameter, double weight, Taxon taxon, TreeModel tree) {
        this(parameter, weight, taxon, tree, null, null);
    }

    public UniformTipFromPrecisionOperator(Parameter parameter, double weight, Taxon taxon, TreeModel tree, Double lowerBound, Double upperBound) {
        this.parameter = parameter;
        this.tree = tree;
        this.taxon = taxon;
        setWeight(weight);

        this.lowerBound = lowerBound;
        this.upperBound = upperBound;

        if (parameter.getDimension() != 1) {
            throw new RuntimeException("UniformTipFromPrecision operator only valid for single tip.");
        }

        if (!parameter.getParameterName().contains(taxon.toString())) {
            throw new RuntimeException("Mismatch between parameter and taxon in UniformTipFromPrecision.");
        }
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return parameter;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() {
        final int index = 0;
        final Bounds<Double> bounds = parameter.getBounds();

        NodeRef node = tree.getExternalNode(tree.getTaxonIndex(taxon));

        final double parentHeight = tree.getNodeHeight(tree.getParent(node));

        final double lower = (lowerBound == null ? bounds.getLowerLimit(index) : Math.max(bounds.getLowerLimit(index), lowerBound));
        final double upper = Math.min(parentHeight, (upperBound == null ? bounds.getUpperLimit(index) : Math.min(bounds.getUpperLimit(index), upperBound)));
        final double newValue = (MathUtils.nextDouble() * (upper - lower)) + lower;

        parameter.setParameterValue(index, newValue);

//        System.out.println(newValue + "[" + lower + "," + upper + "]");
        return 0.0;
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "uniformTipFromPrecision(" + parameter.getParameterName() + ")";
    }

    public final void optimize(double targetProb) {

        throw new RuntimeException("This operator cannot be optimized!");
    }

    public boolean isOptimizing() {
        return false;
    }

    public void setOptimizing(boolean opt) {
        throw new RuntimeException("This operator cannot be optimized!");
    }

    public String getPerformanceSuggestion() {
        return "";
    }

    public String toString() {
        return UniformTipFromPrecisionOperatorParser.UTFP + "(" + parameter.getParameterName() + ")";
    }

    //PRIVATE STUFF

    private Parameter parameter = null;
    private final Double lowerBound;
    private final Double upperBound;
    private final TreeModel tree;
    private final Taxon taxon;
}
