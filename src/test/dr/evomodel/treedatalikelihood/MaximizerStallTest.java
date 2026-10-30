/*
 * MaximizerStallTest.java
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

package test.dr.evomodel.treedatalikelihood;

import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.GradientProvider;
import dr.inference.model.Parameter;
import dr.math.distributions.Distribution;
import dr.math.distributions.NormalDistribution;
import junit.framework.TestCase;

/**
 * A failed line search that no restart improves on counts as convergence to numerical precision only where the
 * gradient is small: a smooth optimum converges; a kink, where the gradient does not shrink, stalls.
 */
public class MaximizerStallTest extends TestCase {

    public MaximizerStallTest(String name) {
        super(name);
    }

    public void testSmoothOptimumConverges() {
        MaximizerWrtParameter maximizer = maximizer(new SmoothNormal(0.3), 2.0);
        maximizer.maximize();
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, maximizer.getReason());
        assertEquals(0.3, maximizer.getMinimumPoint(true)[0], 1E-4);
    }

    public void testKinkStalls() {
        MaximizerWrtParameter maximizer = maximizer(new Kink(0.3), 5.0);
        maximizer.maximize();
        assertTrue(maximizer.getReport(), maximizer.getRestarts() > 0);
        assertEquals(maximizer.getReport(), MaximizerWrtParameter.Reason.LINE_SEARCH_FAILURE, maximizer.getReason());
        assertFalse(maximizer.getReport().contains("to numerical precision"));
        assertTrue(maximizer.getReport(), maximizer.getReport().contains("stalled"));
        // the best point found, near the kink (the line search fails there before reaching it)
        assertEquals(0.3, maximizer.getMinimumPoint(true)[0], 0.1);
    }

    private static MaximizerWrtParameter maximizer(Distribution distribution, double start) {
        Parameter x = new Parameter.Default(new double[]{start});
        DistributionLikelihood likelihood = new DistributionLikelihood(distribution);
        likelihood.addData(x);
        GradientWrtParameterProvider gradient =
                new GradientWrtParameterProvider.ParameterWrapper((GradientProvider) distribution, x, likelihood);
        MaximizerWrtParameter.Settings settings = new MaximizerWrtParameter.Settings(1000, 1E-6, true, false, false)
                .setConvergence(null, 10, 1E-8).setLineSearchRestarts(5);
        return new MaximizerWrtParameter(likelihood, x, gradient, null, settings);
    }

    // log density -(x - m)^2 / 2, up to a constant
    private static final class SmoothNormal extends NormalDistribution implements GradientProvider {
        private final double m;

        SmoothNormal(double m) {
            super(m, 1.0);
            this.m = m;
        }

        public int getDimension() {
            return 1;
        }

        public double[] getGradientLogDensity(Object x) {
            return new double[]{m - ((double[]) x)[0]};
        }
    }

    // log density -(x - m)^2 / 2 - |x - m|: smooth far from the optimum, a kink at it, where the gradient is +-1
    private static final class Kink extends NormalDistribution implements GradientProvider {
        private final double m;

        Kink(double m) {
            super(m, 1.0);
            this.m = m;
        }

        @Override
        public double logPdf(double x) {
            return -0.5 * (x - m) * (x - m) - Math.abs(x - m);
        }

        public int getDimension() {
            return 1;
        }

        public double[] getGradientLogDensity(Object x) {
            final double d = ((double[]) x)[0] - m;
            return new double[]{-d - Math.signum(d)};
        }
    }
}
