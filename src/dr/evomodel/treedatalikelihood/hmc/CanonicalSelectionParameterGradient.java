/*
 * CanonicalSelectionParameterGradient.java
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

package dr.evomodel.treedatalikelihood.hmc;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.Parameter;

/**
 * Canonical tree-level selection-strength gradient for native orthogonal-block parameters.
 */
public final class CanonicalSelectionParameterGradient extends AbstractDiffusionGradient {

    private final ContinuousDataLikelihoodDelegate continuousData;
    private final Parameter parameter;
    private final Parameter rawParameter;
    private final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter;

    public CanonicalSelectionParameterGradient(final TreeDataLikelihood likelihood,
                                               final ContinuousDataLikelihoodDelegate continuousData,
                                               final Parameter parameter,
                                               final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter) {
        super(likelihood, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);
        this.continuousData = continuousData;
        this.parameter = parameter;
        this.rawParameter = nativeBlockParameter.getParameter();
        this.nativeBlockParameter = nativeBlockParameter;
    }

    @Override
    public double[] getGradientLogDensity() {
        return continuousData.getCanonicalSelectionGradient(parameter, nativeBlockParameter);
    }

    @Override
    public double[] getGradientLogDensity(final double[] gradient) {
        if (gradient == null) {
            return getGradientLogDensity();
        }
        return extractSourceGradient(gradient, getDimension());
    }

    @Override
    public Parameter getRawParameter() {
        return rawParameter;
    }

    @Override
    public ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter getDerivationParameter() {
        throw new UnsupportedOperationException(
                "Canonical native selection gradients do not expose legacy derivation metadata.");
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
    public String getReport() {
        return "Gradient." + rawParameter.getParameterName() + "\n" + super.getReport();
    }
}
