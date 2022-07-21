/*
 * UniformTipFromPrecisionOperator.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evomodelxml.operators.UniformTipFromPrecisionOperatorParser;
import dr.inference.operators.UniformOperator;
import dr.math.MathUtils;

/**
 * Uniform operator on tip dates being sampled from precision but that are subject to tree-based age constraints (i.e. parent node height constraint)
 *
 * @author Andy Magee
 * @author Guy Baele
 */
public class UniformTipFromPrecisionOperator extends UniformOperator {

    public UniformTipFromPrecisionOperator(Parameter parameter, double weight, Taxon taxon, TreeModel tree) {
        this(parameter, weight, taxon, tree, null, null);
    }

    public UniformTipFromPrecisionOperator(Parameter parameter, double weight, Taxon taxon, TreeModel tree, Double lowerBound, Double upperBound) {
        super(parameter, weight, lowerBound, upperBound);
        this.tree = tree;
        this.taxon = taxon;

        if (parameter.getDimension() != 1) {
            throw new RuntimeException("UniformTipFromPrecision operator only valid for a single tip.");
        }

        if (!parameter.getParameterName().contains(taxon.toString())) {
            throw new RuntimeException("Mismatch between parameter and taxon names in UniformTipFromPrecision.");
        }
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

        return 0.0;
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "uniformTipFromPrecision(" + parameter.getParameterName() + ")";
    }

    public String toString() {
        return UniformTipFromPrecisionOperatorParser.UTFP + "(" + parameter.getParameterName() + ")";
    }

    //PRIVATE STUFF
    private final TreeModel tree;
    private final Taxon taxon;
}
