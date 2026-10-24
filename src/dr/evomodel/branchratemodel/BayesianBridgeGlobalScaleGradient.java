/*
 * BayesianBridgeGlobalScaleGradient.java
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

package dr.evomodel.branchratemodel;

import dr.inference.distribution.shrinkage.BayesianBridgeDistributionModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * Gradient of an auto-correlated rates prior with a Bayesian-bridge distribution on its increments (the shrinkage
 * local clock) with respect to the bridge's global scale. Only the bridge depends on the global scale.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class BayesianBridgeGlobalScaleGradient implements GradientWrtParameterProvider, Reportable {

    /**
     * Shrinkage-clock gradients are switched off until the clock is supported for gradient-based inference: BEAUti's
     * default exponent (0.25) makes the bridge non-differentiable at zero increments, and node-height gradients do
     * not yet account for the clock's tree-time normalization.
     */
    public static final String SHRINKAGE_CLOCK_NOT_IMPLEMENTED =
            "Gradients for the shrinkage local clock are not yet implemented";

    private static final boolean ENABLED = false;

    private final AutoCorrelatedBranchRatesDistribution ratesPrior;
    private final BayesianBridgeDistributionModel bridge;
    private final double tolerance;

    public BayesianBridgeGlobalScaleGradient(AutoCorrelatedBranchRatesDistribution ratesPrior, double tolerance) {
        if (!ENABLED) {
            throw new UnsupportedOperationException(SHRINKAGE_CLOCK_NOT_IMPLEMENTED);
        }
        if (!(ratesPrior.getPrior() instanceof BayesianBridgeDistributionModel)) {
            throw new IllegalArgumentException("The prior on the increments of " + ratesPrior.getId()
                    + " is not a Bayesian bridge");
        }
        this.ratesPrior = ratesPrior;
        this.bridge = (BayesianBridgeDistributionModel) ratesPrior.getPrior();
        this.tolerance = tolerance;
    }

    @Override
    public Likelihood getLikelihood() {
        return ratesPrior;
    }

    @Override
    public Parameter getParameter() {
        return bridge.getGlobalScale();
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] increments = new double[ratesPrior.getDimension()];
        for (int i = 0; i < increments.length; i++) {
            increments[i] = ratesPrior.getIncrement(i);
        }
        return new double[]{bridge.getGradientLogDensityWrtGlobalScale(increments)};
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, tolerance);
    }
}
