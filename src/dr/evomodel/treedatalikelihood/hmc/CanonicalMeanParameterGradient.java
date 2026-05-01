/*
 * CanonicalMeanParameterGradient.java
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
import dr.evomodel.treedatalikelihood.continuous.canonical.adapter.CanonicalOUGradientAdapter;
import dr.inference.model.Parameter;

/**
 * Canonical tree-level mean gradient for OU stationary/root mean parameters.
 */
public final class CanonicalMeanParameterGradient extends AbstractDiffusionGradient {

    private final CanonicalOUGradientAdapter gradientAdapter;
    private final Parameter parameter;
    private final boolean includeStationaryMean;
    private final boolean includeRootMean;

    public CanonicalMeanParameterGradient(final TreeDataLikelihood likelihood,
                                          final ContinuousDataLikelihoodDelegate continuousData,
                                          final Parameter parameter,
                                          final boolean includeStationaryMean,
                                          final boolean includeRootMean) {
        super(likelihood, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);
        this.gradientAdapter = continuousData.getCanonicalOUGradientAdapter();
        this.parameter = parameter;
        this.includeStationaryMean = includeStationaryMean;
        this.includeRootMean = includeRootMean;
    }

    @Override
    public double[] getGradientLogDensity() {
        return gradientAdapter.getMeanGradient(
                parameter,
                includeStationaryMean,
                includeRootMean);
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
        return parameter;
    }

    @Override
    public ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter getDerivationParameter() {
        if (includeStationaryMean && includeRootMean) {
            return ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_CONSTANT_DRIFT_AND_ROOT_MEAN;
        }
        if (includeStationaryMean) {
            return ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_CONSTANT_DRIFT;
        }
        return ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_ROOT_MEAN;
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
        return "Gradient." + parameter.getParameterName() + "\n" + super.getReport();
    }
}
