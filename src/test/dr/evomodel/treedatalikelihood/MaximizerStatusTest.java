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
import dr.evomodelxml.continuous.hmc.MaximizeWrtParameterParser;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.GammaDistribution;
import dr.util.Transform;
import dr.xml.XMLParseException;
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
                "GradientCalls: ", "Restarts: 0", "LikelihoodRecomputes: ", "NonFiniteEvaluations: 0", "Seconds: ", "GradientNorm: "}) {
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

    public void testMaxRunTime() {
        Parameter x = new Parameter.Default(new double[]{0.01, 20.0});
        MaximizerWrtParameter maximizer = maximizer(x, true, settings(0).setBudgets(1E-12, 0, 0, 0));
        maximizer.maximize();

        assertEquals(Reason.MAX_RUN_TIME, maximizer.getReason());
        assertEquals(1, maximizer.getIterations());
    }

    // maxRunTime takes a number with an optional unit; seconds by default
    public void testParseRunTime() throws XMLParseException {
        assertEquals(90.0, MaximizeWrtParameterParser.parseRunTime("90"), 0.0);
        assertEquals(90.0, MaximizeWrtParameterParser.parseRunTime("90s"), 0.0);
        assertEquals(0.25, MaximizeWrtParameterParser.parseRunTime("250ms"), 1E-15);
        assertEquals(90.0, MaximizeWrtParameterParser.parseRunTime("1.5min"), 0.0);
        assertEquals(120.0, MaximizeWrtParameterParser.parseRunTime(" 2 m "), 0.0);
        assertEquals(7200.0, MaximizeWrtParameterParser.parseRunTime("2 hours"), 0.0);
        assertEquals(3600.0, MaximizeWrtParameterParser.parseRunTime("1H"), 0.0);
        assertEquals(86400.0, MaximizeWrtParameterParser.parseRunTime("1d"), 0.0);
        assertEquals(1000.0, MaximizeWrtParameterParser.parseRunTime("1e3s"), 0.0);
        for (String bad : new String[]{"", "fast", "10 weeks", "-5s", "1.2.3s", "s10"}) {
            try {
                MaximizeWrtParameterParser.parseRunTime(bad);
                fail("\"" + bad + "\" should not parse");
            } catch (XMLParseException expected) {
                // expected
            }
        }
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

    // a holds at 1.5, above its mode 1.0, once a step takes it below; b, not held, still reaches the mode
    public void testHoldAtLowerBound() {
        Parameter a = new Parameter.Default(new double[]{3.0});
        Parameter b = new Parameter.Default(new double[]{0.4});
        dr.inference.model.CompoundParameter x = new dr.inference.model.CompoundParameter("x");
        x.addParameter(a);
        x.addParameter(b);
        MaximizerWrtParameter maximizer = maximizer(x, true, settings(0).addLowerBoundHold(a, 1.5)
                .setCheckGradientEvery(1));
        maximizer.maximize();

        assertEquals(Reason.CONVERGED, maximizer.getReason());
        assertEquals(1.5, a.getParameterValue(0), 0.0);
        assertEquals(MODE, b.getParameterValue(0), 1E-5);
        assertEquals(1, maximizer.getHeld().size());
        // d/dlog a of log Gamma(a; 3, 0.5) = 2 - 2a: negative, the optimum lies below the bound
        assertEquals(2 - 2 * 1.5, maximizer.getHeld().get(0).gradient, 1E-6);
        assertTrue(maximizer.getReport().contains("HeldAtLowerBound: "));
        for (double error : maximizer.getGradientCheckErrors()) {
            assertTrue(error < 1E-6);
        }

        // a bound the iterates never cross changes nothing
        a.setParameterValue(0, 3.0);
        maximizer = maximizer(x, true, settings(0).addLowerBoundHold(a, 0.5));
        maximizer.maximize();
        assertEquals(MODE, a.getParameterValue(0), 1E-5);
        assertTrue(maximizer.getHeld().isEmpty());

        try {
            maximizer(x, true, settings(0).addLowerBoundHold(new Parameter.Default(1.0), 0.5));
            fail("a parameter outside the maximized one");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().contains("holdAtLowerBound"));
        }
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
