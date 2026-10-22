/*
 * SimplexTangentCheckTest.java
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

package test.dr.inference.hmc;

import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.GradientWrtParameterProvider.CheckGradientNumerically;
import dr.inference.hmc.GradientWrtParameterProvider.MismatchException;
import dr.inference.hmc.JointGradient;
import dr.inference.hmc.SimplexTangent;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.DirichletDistribution;
import dr.util.UnitSimplexToRealsTransform;
import junit.framework.TestCase;

import java.util.Arrays;

/**
 * Gradient checks for parameters on a simplex compare only the tangent part g - (x.g) 1 / S (plan 5, FQ4 a-c).
 */
public class SimplexTangentCheckTest extends TestCase {

    private static final double[] X = {0.1, 0.2, 0.3, 0.4};
    private static final double[] COUNTS = {3, 2, 5, 4};

    public SimplexTangentCheckTest(String name) {
        super(name);
    }

    public void testToyExtensionFailsRawPassesTangent() {
        // f = sum_k c_k log x_k + c0 log sum x; the analytic side drops the c0 term, which only adds c0 / S * 1
        final double[] c = {1.0, 2.0, 3.0, 4.0};
        final double c0 = 5.0;

        Parameter raw = new Parameter.Default(X.clone());
        try {
            GradientWrtParameterProvider.getReportAndCheckForError(toyProvider(raw, c, c0, Double.NaN),
                    0.0, Double.POSITIVE_INFINITY, 1E-6);
            fail("the coordinate-wise check must see the extension-dependent c0 term");
        } catch (RuntimeException expected) {
            // expected
        }

        Parameter onSimplex = new Parameter.Default(X.clone());
        String report = GradientWrtParameterProvider.getReportAndCheckForError(toyProvider(onSimplex, c, c0, 1.0),
                0.0, Double.POSITIVE_INFINITY, 1E-6);
        assertTrue(report.contains("analytic (tangent)"));
        assertTrue(report.contains("numeric (raw, extension-dependent, unchecked)"));
    }

    public void testDirichletTangent() throws MismatchException {
        Parameter x = new Parameter.Default(X.clone());
        GradientWrtParameterProvider dirichlet = dirichletGradient(x);
        assertEquals(1.0, dirichlet.getSimplexSum());

        double[] expected = {10.0, -5.0, 10.0 / 3.0, -2.5};
        double[] analytic = SimplexTangent.project(
                dirichlet.getGradientLogDensity(), x.getParameterValues(), 1.0);
        double[] numeric = checker(dirichlet, 1E-6).getTangentNumerical(1.0);
        for (int k = 0; k < expected.length; k++) {
            assertEquals(expected[k], analytic[k], 1E-12);
            assertEquals(expected[k], numeric[k], 1E-6);
        }

        String report = GradientWrtParameterProvider.getReportAndCheckForError(dirichlet,
                0.0, 1.0, 1E-6);
        assertTrue(report.contains("numeric (raw): skipped (Dirichlet off-simplex)"));
    }

    public void testRestoresParameterAndLikelihood() throws MismatchException {
        Parameter x = new Parameter.Default(X.clone());
        GradientWrtParameterProvider dirichlet = dirichletGradient(x);
        double[] before = x.getParameterValues();
        double logLikelihood = dirichlet.getLikelihood().getLogLikelihood();

        checker(dirichlet, 1E-6).getTangentNumerical(1.0);

        double[] after = x.getParameterValues();
        for (int k = 0; k < before.length; k++) {
            assertEquals(Double.doubleToLongBits(before[k]), Double.doubleToLongBits(after[k]));
        }
        assertEquals(Double.doubleToLongBits(logLikelihood),
                Double.doubleToLongBits(dirichlet.getLikelihood().getLogLikelihood()));
    }

    public void testOffSimplexPrecondition() {
        double[] off = X.clone();
        off[3] += 1E-9;
        Parameter x = new Parameter.Default(off);
        try {
            checker(dirichletGradient(x), 1E-6).getTangentNumerical(1.0);
            fail("a point with sum 1 + 1e-9 must be refused");
        } catch (MismatchException expected) {
            assertTrue(expected.getMessage().startsWith("off-simplex"));
        }
    }

    public void testNonFiniteAndLengthMismatchFailCheckedReports() throws MismatchException {
        double[] finite = {1.0, 2.0};
        try {
            GradientWrtParameterProvider.makeReport("", finite, new double[]{1.0, Double.NaN}, true, 0.1, 0.0);
            fail("a NaN in a checked line must fail");
        } catch (MismatchException expected) {
            // expected
        }
        try {
            GradientWrtParameterProvider.makeReport("", finite, new double[]{1.0}, true, 0.1, 0.0);
            fail("a length mismatch in a checked line must fail");
        } catch (MismatchException expected) {
            // expected
        }
        // unchecked reports only print
        GradientWrtParameterProvider.makeReport("", finite, new double[]{1.0, Double.NaN}, false, 0.1, 0.0);
    }

    public void testJointGradientSimplexSum() {
        Parameter x = new Parameter.Default(X.clone());
        GradientWrtParameterProvider first = dirichletGradient(x);
        GradientWrtParameterProvider second = dirichletGradient(x);
        assertEquals(1.0, new JointGradient(Arrays.asList(first, second)).getSimplexSum());

        GradientWrtParameterProvider unmarked = toyProvider(x, new double[]{1, 1, 1, 1}, 0.0, Double.NaN);
        assertTrue(Double.isNaN(new JointGradient(Arrays.asList(first, unmarked)).getSimplexSum()));
    }

    public void testStickBreakingChainIgnoresConstantShift() {
        // J^T g against finite differences of f(pi(y)), where f is the Dirichlet log density
        final int dim = X.length;
        UnitSimplexToRealsTransform transform = new UnitSimplexToRealsTransform(dim);
        DirichletDistribution dirichlet = new DirichletDistribution(COUNTS, 1.0);

        double[] y = transform.transform(X, 0, dim);
        double[] g = dirichlet.getGradientLogDensity(X);
        double[][] jacobian = transform.computeJacobianMatrixInverse(y); // jacobian[j][i] = d pi_i / d y_j

        final double h = 1E-6;
        for (int j = 0; j < dim - 1; j++) {
            double chain = 0.0;
            double shifted = 0.0;
            for (int i = 0; i < dim; i++) {
                chain += jacobian[j][i] * g[i];
                shifted += jacobian[j][i] * (g[i] + 7.0);
            }

            double[] plus = y.clone();
            plus[j] += h;
            double[] minus = y.clone();
            minus[j] -= h;
            double numeric = (dirichlet.logPdf(transform.inverse(plus, 0, dim))
                    - dirichlet.logPdf(transform.inverse(minus, 0, dim))) / (2 * h);

            assertEquals("y" + j, numeric, chain, 1E-6 * Math.max(1.0, Math.abs(numeric)));
            assertEquals("shift y" + j, chain, shifted, 1E-10);
        }
    }

    private static CheckGradientNumerically checker(GradientWrtParameterProvider provider, double tolerance) {
        return new CheckGradientNumerically(provider, 0.0, Double.POSITIVE_INFINITY, tolerance, null);
    }

    private static GradientWrtParameterProvider dirichletGradient(Parameter x) {
        DirichletDistribution distribution = new DirichletDistribution(COUNTS, 1.0);
        MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(distribution);
        likelihood.addData(x);
        return new GradientWrtParameterProvider.ParameterWrapper(distribution, x, likelihood);
    }

    private static GradientWrtParameterProvider toyProvider(final Parameter x, final double[] c, final double c0,
                                                            final double simplexSum) {

        final Likelihood likelihood = new Likelihood.Abstract(null) {
            @Override
            protected boolean getLikelihoodKnown() {
                return false;
            }

            @Override
            protected double calculateLogLikelihood() {
                double logL = 0.0;
                double sum = 0.0;
                for (int k = 0; k < c.length; k++) {
                    logL += c[k] * Math.log(x.getParameterValue(k));
                    sum += x.getParameterValue(k);
                }
                return logL + c0 * Math.log(sum);
            }
        };

        return new GradientWrtParameterProvider() {
            @Override
            public Likelihood getLikelihood() {
                return likelihood;
            }

            @Override
            public Parameter getParameter() {
                return x;
            }

            @Override
            public int getDimension() {
                return x.getDimension();
            }

            @Override
            public double[] getGradientLogDensity() {
                double[] gradient = new double[c.length];
                for (int k = 0; k < c.length; k++) {
                    gradient[k] = c[k] / x.getParameterValue(k);
                }
                return gradient;
            }

            @Override
            public double getSimplexSum() {
                return simplexSum;
            }
        };
    }
}
