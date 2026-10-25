/*
 * MaximizerStatusTest.java
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

import com.github.lbfgs4j.liblbfgs.LbfgsConstant.ReturnValue;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter.Reason;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.GammaDistribution;
import dr.util.Transform;
import junit.framework.TestCase;

/**
 * The maximizer's stop reasons, budgets, caching, non-finite guard, gradient check and report (plan 1, section 5).
 * The objective is a product of Gamma(3, 0.5) densities, whose mode is 1.0 in every coordinate.
 */
public class MaximizerStatusTest extends TestCase {

    private static final double SHAPE = 3.0;
    private static final double SCALE = 0.5;
    private static final double MODE = (SHAPE - 1) * SCALE;

    private int likelihoodCalls;
    private int gradientCalls;

    public MaximizerStatusTest(String name) {
        super(name);
    }

    public void testConvergedReport() {
        Parameter x = new Parameter.Default(new double[]{0.4, 2.2, 5.0});
        MaximizerWrtParameter maximizer = maximizer(x, true, settings(0));

        assertFalse(maximizer.wasExecuted());
        assertEquals("Not yet executed.", maximizer.getReport());

        maximizer.maximize();

        assertTrue(maximizer.wasExecuted());
        assertEquals(ReturnValue.LBFGS_SUCCESS, maximizer.getStatus());
        assertEquals(Reason.CONVERGED, maximizer.getReason());
        assertTrue(maximizer.getIterations() > 0);
        for (int i = 0; i < x.getDimension(); i++) {
            assertEquals(MODE, x.getParameterValue(i), 1E-5);
        }

        String report = maximizer.getReport();
        for (String key : new String[]{"Untransformed X: ", "X: ", "Gradient: ", "Gradient type: analytic", "Fx: ",
                "Time: ", "Count: ", "Status: LBFGS_SUCCESS", "Reason: converged", "Iterations: ", "ValueCalls: ",
                "GradientCalls: ", "LikelihoodRecomputes: ", "NonFiniteEvaluations: 0", "Seconds: ", "GradientNorm: "}) {
            assertTrue("report has " + key, report.contains(key));
        }
    }

    // L-BFGS asks for the value and the gradient at the same point: each is computed once, and the report computes
    // nothing
    public void testOneEvaluationPerPoint() {
        Parameter x = new Parameter.Default(new double[]{0.4, 2.2, 5.0});
        MaximizerWrtParameter maximizer = maximizer(x, true, settings(0));
        maximizer.maximize();

        assertEquals(maximizer.getValueCalls(), maximizer.getGradientCalls());
        assertEquals(maximizer.getValueCalls(), likelihoodCalls);
        assertEquals(maximizer.getGradientCalls(), gradientCalls);

        maximizer.getReport();
        maximizer.getReport();
        assertEquals(maximizer.getValueCalls(), likelihoodCalls);
        assertEquals(maximizer.getGradientCalls(), gradientCalls);
    }

    public void testMaxIterations() {
        Parameter x = new Parameter.Default(new double[]{0.01, 20.0});
        MaximizerWrtParameter maximizer = maximizer(x, true, settings(2));
        maximizer.maximize();

        assertEquals(ReturnValue.LBFGSERR_MAXIMUMITERATION, maximizer.getStatus());
        assertEquals(Reason.MAX_ITERATIONS, maximizer.getReason());
        assertEquals(2, maximizer.getIterations());
    }

    public void testMaxEvaluations() {
        Parameter x = new Parameter.Default(new double[]{0.01, 20.0});
        MaximizerWrtParameter maximizer = maximizer(x, true, settings(0).setBudgets(0, 3, 0, 0));
        maximizer.maximize();

        assertEquals(ReturnValue.LBFGS_STOP, maximizer.getStatus());
        assertEquals(Reason.MAX_EVALUATIONS, maximizer.getReason());
        assertTrue(maximizer.getValueCalls() >= 3);
        assertTrue(maximizer.getReport().contains("Reason: maxEvaluations"));
    }

    public void testMaxSeconds() {
        Parameter x = new Parameter.Default(new double[]{0.01, 20.0});
        MaximizerWrtParameter maximizer = maximizer(x, true, settings(0).setBudgets(1E-12, 0, 0, 0));
        maximizer.maximize();

        assertEquals(Reason.MAX_SECONDS, maximizer.getReason());
        assertEquals(1, maximizer.getIterations());
    }

    public void testStall() {
        Parameter x = new Parameter.Default(new double[]{0.01, 20.0});
        // any improvement smaller than 1E10 relative counts as a stall
        MaximizerWrtParameter maximizer = maximizer(x, true, settings(0).setBudgets(0, 0, 1E10, 1));
        maximizer.maximize();

        assertEquals(Reason.STALL, maximizer.getReason());
        assertEquals(2, maximizer.getIterations());
    }

    // without a transform the first (unit) step from 0.5 lands at -0.5, outside the domain of Gamma(3, 0.1); the
    // non-finite value makes the line search back off
    public void testNonFiniteGuard() {
        Parameter x = new Parameter.Default(new double[]{0.5});
        MaximizerWrtParameter maximizer = maximizer(x, false, 0.1, settings(0)
                .setConvergence(1E-8, null, null));
        maximizer.maximize();

        assertTrue(maximizer.getNonFiniteEvaluations() > 0);
        assertEquals(Reason.CONVERGED, maximizer.getReason());
        assertEquals((SHAPE - 1) * 0.1, x.getParameterValue(0), 1E-6);
        assertTrue(maximizer.getReport().contains("NonFiniteEvaluations: " + maximizer.getNonFiniteEvaluations()));
    }

    public void testGradientCheck() {
        Parameter x = new Parameter.Default(new double[]{0.4, 2.2, 5.0});
        MaximizerWrtParameter maximizer = maximizer(x, true, settings(0).setCheckGradientEvery(1));
        maximizer.maximize();

        assertEquals(Reason.CONVERGED, maximizer.getReason());
        for (int i = 0; i < x.getDimension(); i++) {
            assertEquals(MODE, x.getParameterValue(i), 1E-5);
        }
        String report = maximizer.getReport();
        String line = report.substring(report.indexOf("WorstGradientCheckError: "));
        double error = Double.parseDouble(line.substring(line.indexOf(' ') + 1, line.indexOf('\n')));
        assertTrue("worst gradient check error " + error, error < 1E-6);
    }

    public void testNumericalGradientLeavesTheParameter() {
        Parameter x = new Parameter.Default(new double[]{0.4, 2.2, 5.0});
        GammaDistribution gamma = new GammaDistribution(SHAPE, SCALE);
        DistributionLikelihood likelihood = new DistributionLikelihood(gamma);
        likelihood.addData(x);
        MaximizerWrtParameter maximizer = new MaximizerWrtParameter(likelihood, x, null,
                new Transform.Array(Transform.LOG, x.getDimension(), null), settings(0));

        double[] gradient = maximizer.getGradient().getGradientLogDensity();
        assertEquals(0.4, x.getParameterValue(0), 0.0);
        assertEquals(2.2, x.getParameterValue(1), 0.0);
        assertEquals(5.0, x.getParameterValue(2), 0.0);
        // d/dx log Gamma(x; 3, 0.5) = 2 / x - 2
        assertEquals(2 / 5.0 - 2, gradient[2], 1E-5);

        maximizer.maximize();
        for (int i = 0; i < x.getDimension(); i++) {
            assertEquals(MODE, x.getParameterValue(i), 1E-4);
        }
        assertTrue(maximizer.getReport().contains("Gradient type: numerical"));
    }

    public void testWarnings() {
        Parameter x = new Parameter.Default(new double[]{0.4});
        x.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        String[] warnings = maximizer(x, false,
                new MaximizerWrtParameter.Settings(0, 1E-6, false, false, false)).getWarnings();
        assertEquals(3, warnings.length);
        assertTrue(warnings[0].contains("nIterations"));
        assertTrue(warnings[1].contains("no transform"));
        assertTrue(warnings[2].contains("startAtCurrentState"));

        assertEquals(0, maximizer(x, true, settings(100)).getWarnings().length);
    }

    private static MaximizerWrtParameter.Settings settings(int iterations) {
        return new MaximizerWrtParameter.Settings(iterations, 1E-6, true, false, false);
    }

    private MaximizerWrtParameter maximizer(Parameter x, boolean logTransform, MaximizerWrtParameter.Settings settings) {
        return maximizer(x, logTransform, SCALE, settings);
    }

    private MaximizerWrtParameter maximizer(Parameter x, boolean logTransform, double scale,
                                            MaximizerWrtParameter.Settings settings) {
        GammaDistribution gamma = new GammaDistribution(SHAPE, scale);
        DistributionLikelihood likelihood = new DistributionLikelihood(gamma) {
            @Override
            public double calculateLogLikelihood() {
                ++likelihoodCalls;
                return super.calculateLogLikelihood();
            }
        };
        likelihood.addData(x);
        final GradientWrtParameterProvider wrapped = new GradientWrtParameterProvider.ParameterWrapper(gamma, x, likelihood);
        GradientWrtParameterProvider gradient = new GradientWrtParameterProvider() {
            public Likelihood getLikelihood() { return wrapped.getLikelihood(); }
            public Parameter getParameter() { return wrapped.getParameter(); }
            public int getDimension() { return wrapped.getDimension(); }
            public double[] getGradientLogDensity() {
                ++gradientCalls;
                return wrapped.getGradientLogDensity();
            }
        };
        Transform transform = logTransform ? new Transform.Array(Transform.LOG, x.getDimension(), null) : null;
        return new MaximizerWrtParameter(likelihood, x, gradient, transform, settings);
    }
}
