/*
 * CTMCScalePriorNodeHeightGradient.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.tree;

import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * Gradient of the CTMC rate reference prior with respect to the internal node heights, which enter through the
 * tree length.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class CTMCScalePriorNodeHeightGradient implements GradientWrtParameterProvider, Reportable {

    private final CTMCScalePrior prior;
    private final TreeModel tree;
    private final Parameter parameter;
    private final double tolerance;

    public CTMCScalePriorNodeHeightGradient(CTMCScalePrior prior, TreeModel tree, double tolerance) {
        this.prior = prior;
        this.tree = tree;
        this.parameter = new NodeHeightProxyParameter("nodeHeightProxyParameter", tree, true);
        this.tolerance = tolerance;
    }

    @Override
    public Likelihood getLikelihood() {
        return prior;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return prior.getGradientLogDensityWrtNodeHeights(tree);
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, tolerance);
    }
}
