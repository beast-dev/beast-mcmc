/*
 * LaplaceEvidenceTest.java
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

import dr.evomodel.treedatalikelihood.ApproximateTreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.LaplaceApproximation;
import dr.evomodel.treedatalikelihood.discrete.LaplaceEvidence;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.GradientProvider;
import dr.inference.model.Parameter;
import dr.math.distributions.BetaDistribution;
import dr.math.distributions.Distribution;
import dr.math.distributions.LogNormalDistribution;
import dr.math.distributions.NormalDistribution;
import dr.util.Transform;
import junit.framework.TestCase;

/**
 * The Laplace evidence (plan 8 X3) on densities whose normalizing constants are known, through both classes that use
 * it: a Gaussian is exact (log Z = 0 for a normalized density, so (1/2) log 2 pi + log s for exp(-x^2/2s^2)); a
 * log-normal is exactly Gaussian in the log coordinates, once the Jacobian is included; a Beta in logit coordinates
 * is not, and gives the closed-form Laplace value. The earlier code added +(1/2) log det of the precision, which gives
 * (1/2) log 2 pi - log s, and took the Hessian of the untransformed parameter.
 */
public class LaplaceEvidenceTest extends TestCase {

    public LaplaceEvidenceTest(String name) {
        super(name);
    }

    // a normalized N(3, 2^2): log Z = 0, and the curvature is 1/s^2
    public void testGaussian() {
        double sd = 2.0;
        MaximizerWrtParameter maximizer = maximizer(new NormalDistribution(3.0, sd), new double[]{0.5}, null, false);
        maximizer.maximize();
        LaplaceEvidence laplace = new LaplaceEvidence(maximizer);
        assertEquals(1.0 / (sd * sd), laplace.getPrecision()[0][0], 1E-6);
        assertEquals(0.0, laplace.getLogEvidence(), 1E-6);
        // the unnormalized exp(-x^2 / 2 s^2): log Z = (1/2) log 2 pi + log s
        double normalizer = -0.5 * Math.log(2 * Math.PI) - Math.log(sd);
        assertEquals(0.5 * Math.log(2 * Math.PI) + Math.log(sd), laplace.getLogEvidence() - normalizer, 1E-6);

        // both classes that use it
        LaplaceApproximation approximation = new LaplaceApproximation(maximizer, true, false, true);
        approximation.approximate();
        assertEquals(0.0, approximation.getLogMarginalLikelihood(), 1E-6);
        assertTrue(approximation.getReport().contains("sd(y) = "));
        assertEquals(0.0, new ApproximateTreeDataLikelihood(
                maximizer(new NormalDistribution(3.0, sd), new double[]{0.5}, null, false)).getLogLikelihood(), 1E-6);
    }

    // two log-normals maximized in log coordinates with the Jacobian: Gaussian in y, so exact
    public void testLogTransformedGaussian() {
        Distribution logNormal = new LogNormalDistribution(0.3, 0.6);
        MaximizerWrtParameter maximizer = maximizer(logNormal, new double[]{2.0, 0.4},
                new Transform.Array(Transform.LOG, 2, null), true);
        maximizer.maximize();
        LaplaceEvidence laplace = new LaplaceEvidence(maximizer);
        // independent coordinates, each of precision 1 / 0.6^2
        assertEquals(1.0 / 0.36, laplace.getPrecision()[0][0], 1E-5);
        assertEquals(0.0, laplace.getPrecision()[0][1], 1E-6);
        assertEquals(0.0, laplace.getLogEvidence(), 1E-6);
        assertEquals(0.0, new ApproximateTreeDataLikelihood(maximizer(logNormal, new double[]{2.0, 0.4},
                new Transform.Array(Transform.LOG, 2, null), true)).getLogLikelihood(), 1E-6);

        // without the Jacobian the optimum is not the mode of the density in y: no evidence
        MaximizerWrtParameter unweighted = maximizer(logNormal, new double[]{2.0, 0.4},
                new Transform.Array(Transform.LOG, 2, null), false);
        unweighted.maximize();
        try {
            new LaplaceEvidence(unweighted).getLogEvidence();
            fail("the evidence needs includeJacobian");
        } catch (IllegalStateException expected) {
        }
    }

    // Beta(a, b) in logit coordinates: l(y) = a log x + b log(1 - x) - log B(a, b), maximized at x = a / (a + b), with
    // curvature (a + b) x (1 - x); the Laplace value is close to log Z = 0 and equal to its closed form
    public void testLogitBeta() {
        final double a = 8.0;
        final double b = 5.0;
        MaximizerWrtParameter maximizer = maximizer(new BetaWithGradient(a, b), new double[]{0.2},
                new Transform.Array(Transform.LOGIT, 1, null), true);
        maximizer.maximize();
        double x = a / (a + b);
        double curvature = (a + b) * x * (1 - x);
        double expected = new BetaDistribution(a, b).logPdf(x) + Math.log(x * (1 - x)) +
                0.5 * Math.log(2 * Math.PI) - 0.5 * Math.log(curvature);
        LaplaceEvidence laplace = new LaplaceEvidence(maximizer);
        assertEquals(curvature, laplace.getPrecision()[0][0], 1E-5);
        assertEquals(expected, laplace.getLogEvidence(), 1E-6);
        assertTrue(Math.abs(laplace.getLogEvidence()) < 0.05);
    }

    // a negative Hessian that is not positive definite (not a maximum, or flat) has no log determinant
    public void testNotPositiveDefinite() {
        assertEquals(2.0 * Math.log(2.0), logDet(new double[][]{{4.0, 0.0}, {0.0, 1.0}}), 1E-12);
        try {
            LaplaceEvidence.cholesky(new double[][]{{1.0, 2.0}, {2.0, 1.0}});
            fail("not positive definite");
        } catch (IllegalStateException expected) {
        }
    }

    private static double logDet(double[][] matrix) {
        double[][] lower = LaplaceEvidence.cholesky(matrix);
        double sum = 0.0;
        for (int i = 0; i < lower.length; i++) {
            sum += 2.0 * Math.log(lower[i][i]);
        }
        return sum;
    }

    // -------------------------------------------------------------------------------------------------------------

    private static MaximizerWrtParameter maximizer(Distribution distribution, double[] start, Transform transform,
                                                   boolean includeJacobian) {
        Parameter x = new Parameter.Default(start);
        DistributionLikelihood likelihood = new DistributionLikelihood(distribution);
        likelihood.addData(x);
        GradientWrtParameterProvider gradient =
                new GradientWrtParameterProvider.ParameterWrapper((GradientProvider) distribution, x, likelihood);
        MaximizerWrtParameter.Settings settings =
                new MaximizerWrtParameter.Settings(0, 1E-6, true, false, includeJacobian)
                        .setConvergence(1E-10, null, null);
        return new MaximizerWrtParameter(likelihood, x, gradient, transform, settings);
    }

    // BetaDistribution has no gradient
    private static final class BetaWithGradient extends BetaDistribution implements GradientProvider {
        private final double alpha;
        private final double beta;

        BetaWithGradient(double alpha, double beta) {
            super(alpha, beta);
            this.alpha = alpha;
            this.beta = beta;
        }

        public int getDimension() {
            return 1;
        }

        public double[] getGradientLogDensity(Object x) {
            double[] values = (double[]) x;
            double[] gradient = new double[values.length];
            for (int i = 0; i < values.length; i++) {
                gradient[i] = (alpha - 1) / values[i] - (beta - 1) / (1 - values[i]);
            }
            return gradient;
        }
    }
}
