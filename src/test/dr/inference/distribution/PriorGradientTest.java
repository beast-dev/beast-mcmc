/*
 * PriorGradientTest.java
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

package test.dr.inference.distribution;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.math.distributions.Distribution;
import dr.math.distributions.ExponentialDistribution;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.LaplaceDistribution;
import dr.math.distributions.LogNormalDistribution;
import dr.math.distributions.NormalDistribution;
import dr.math.distributions.UniformDistribution;
import dr.inference.model.GradientProvider;
import dr.util.Transform;
import junit.framework.TestCase;

/**
 * Gradients of univariate priors through {@code <gradient>}, including the offset and scale of the prior
 * (plan 3, B5, B6, B-exp, B-uniform), against central finite differences of the prior's log density.
 */
public class PriorGradientTest extends TestCase {

    private static final double H = 1E-6;
    private static final double TOLERANCE = 1E-7;

    public PriorGradientTest(String name) {
        super(name);
    }

    public void testNormal() {
        check(new NormalDistribution(0.5, 2.0), 0.0, 1.0, -1.3, 0.2, 3.1);
        check(new NormalDistribution(0.5, 2.0), 1.5, 3.0, -1.3, 0.2, 3.1); // offset and scale
    }

    public void testLogNormalWithOffset() {
        check(new LogNormalDistribution(-0.5, 0.8), 0.0, 1.0, 0.3, 1.7);
        check(new LogNormalDistribution(-0.5, 0.8), 2.0, 1.0, 2.3, 3.7);
    }

    public void testGammaWithOffset() {
        check(new GammaDistribution(2.5, 0.4), 0.0, 1.0, 0.3, 1.7);
        check(new GammaDistribution(2.5, 0.4), 2.0, 1.0, 2.3, 3.7);
    }

    public void testExponentialWithOffset() {
        check(new ExponentialDistribution(1.0 / 0.7), 0.0, 1.0, 0.3, 1.7);
        check(new ExponentialDistribution(1.0 / 0.7), 2.0, 1.0, 2.3, 3.7);
    }

    public void testLaplace() {
        // on both sides of the mean: the derivative of -|x - m| / b is -sign(x - m) / b
        check(new LaplaceDistribution(0.5, 0.3), 0.0, 1.0, -1.0, 0.2, 0.8, 2.0);
        double[] gradient = ((GradientProvider) new LaplaceDistribution(0.5, 0.3))
                .getGradientLogDensity(new double[]{0.2, 0.8});
        assertEquals(1.0 / 0.3, gradient[0], 1E-15);
        assertEquals(-1.0 / 0.3, gradient[1], 1E-15);
    }

    public void testUniform() {
        check(new UniformDistribution(-1.0, 4.0), 0.0, 1.0, -0.5, 1.0, 3.5);
    }

    // y = log((x - l) / (u - x)); the transformed density is f(x(y)) + log dx/dy
    public void testScaledLogitTransform() {
        final double lower = -1.0;
        final double upper = 4.0;
        Transform.ScaledLogitTransform transform = new Transform.ScaledLogitTransform(upper, lower);
        NormalDistribution f = new NormalDistribution(0.5, 2.0);

        for (double x : new double[]{-0.9, 0.3, 2.0, 3.95}) {
            final double y = transform.transform(x);
            assertEquals(x, transform.inverse(y), 1E-12);

            final double numericJacobian = (transform.inverse(y + H) - transform.inverse(y - H)) / (2 * H);
            assertEquals("dx/dy at x=" + x, numericJacobian, transform.gradientInverse(y), 1E-8);
            assertEquals("d log J/dy at x=" + x,
                    (logJacobian(transform, y + H) - logJacobian(transform, y - H)) / (2 * H),
                    transform.gradientLogJacobianInverse(y), 1E-6);

            final double gradient = NormalDistribution.gradLogPdf(x, 0.5, 2.0);
            final double hessian = -1.0 / (2.0 * 2.0);
            final double analytic = transform.updateGradientLogDensity(gradient, x);
            final double numeric = (transformedLogDensity(transform, f, y + H)
                    - transformedLogDensity(transform, f, y - H)) / (2 * H);
            assertEquals("gradient at x=" + x, numeric, analytic, 1E-6 * Math.max(1.0, Math.abs(numeric)));

            final double h = 1E-4;
            final double numericHessian = (transformedLogDensity(transform, f, y + h)
                    - 2 * transformedLogDensity(transform, f, y)
                    + transformedLogDensity(transform, f, y - h)) / (h * h);
            assertEquals("Hessian at x=" + x, numericHessian,
                    transform.updateDiagonalHessianLogDensity(hessian, gradient, x),
                    1E-5 * Math.max(1.0, Math.abs(numericHessian)));
        }
    }

    private static double logJacobian(Transform.ScaledLogitTransform transform, double y) {
        return Math.log(transform.gradientInverse(y));
    }

    private static double transformedLogDensity(Transform.ScaledLogitTransform transform, Distribution f, double y) {
        return f.logPdf(transform.inverse(y)) + logJacobian(transform, y);
    }

    private static void check(Distribution distribution, double offset, double scale, double... values) {
        Parameter parameter = new Parameter.Default(values);
        DistributionLikelihood likelihood = new DistributionLikelihood(distribution, offset, false, scale);
        likelihood.addData(parameter);
        GradientWrtParameterProvider gradient = new GradientWrtParameterProvider.ParameterWrapper(
                (GradientProvider) distribution, parameter, likelihood);

        double[] analytic = gradient.getGradientLogDensity();
        for (int i = 0; i < values.length; i++) {
            parameter.setParameterValue(i, values[i] + H);
            final double plus = likelihood.calculateLogLikelihood();
            parameter.setParameterValue(i, values[i] - H);
            final double minus = likelihood.calculateLogLikelihood();
            parameter.setParameterValue(i, values[i]);
            final double numeric = (plus - minus) / (2 * H);
            assertEquals(distribution.getClass().getSimpleName() + " offset=" + offset + " scale=" + scale
                    + " x=" + values[i], numeric, analytic[i], TOLERANCE * Math.max(1.0, Math.abs(numeric)));
        }
    }
}
