/*
 * MaximizerBudgetTest.java
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
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter.Budget;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter.Reason;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.MultivariateNormalDistribution;
import junit.framework.TestCase;

/**
 * Per-call budgets and the call-global best point of the maximizer (plan 7 M0, V15). The objective is a narrow
 * Gaussian valley, F(a, b) = ((a - 1)^2 + 100 (b - a)^2) / 2, with a held at 1.5 once an iterate takes it below:
 * moving a up to the bound from an iterate on the valley floor near (1, 1) raises F far above the start at (3, 3), so a
 * budget that runs out right after the hold used to end the call there.
 */
public class MaximizerBudgetTest extends TestCase {

    private static final double BOUND = 1.5;

    private Parameter a;
    private Parameter b;
    private CompoundParameter x;

    public MaximizerBudgetTest(String name) {
        super(name);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
        a = new Parameter.Default("a", 3.0);
        b = new Parameter.Default("b", 3.0);
        x = new CompoundParameter("x");
        x.addParameter(a);
        x.addParameter(b);
    }

    // every iteration budget ends no higher than the start and on the best accepted point; for some, a hold comes
    // right before the budget runs out, and the call ends on the earlier best point with the late hold released
    public void testHoldThenBudgetEndsAtBest() {
        MaximizerWrtParameter maximizer = maximizer(false, holdSettings(0));
        int endedAtBest = 0;
        for (int iterations = 1; iterations <= 40; iterations++) {
            a.setParameterValue(0, 3.0);
            b.setParameterValue(0, 3.0);
            maximizer.maximize(new Budget(iterations).setQuiet(true));

            final double end = maximizer.getMinimumValue();
            final double tolerance = 1E-8 * Math.max(1.0, Math.abs(end));
            assertEquals(objective(3.0, 3.0), maximizer.getStartValue(), 1E-12);
            assertTrue(iterations + ": ended at " + end + " above the start " + maximizer.getStartValue(),
                    end <= maximizer.getStartValue() + tolerance);
            assertTrue(iterations + ": ended at " + end + " above the best " + maximizer.getBestValue(),
                    end <= maximizer.getBestValue() + tolerance);
            // the model is left at the point reported
            assertEquals(end, objective(a.getParameterValue(0), b.getParameterValue(0)), 1E-10);
            assertTrue(a.getParameterValue(0) >= BOUND);
            if (maximizer.endedAtBestPoint()) {
                endedAtBest++;
                assertTrue(maximizer.getHeld().isEmpty() || a.getParameterValue(0) == BOUND);
                assertTrue(maximizer.getReport().contains("EndedAtBestPoint: true"));
            }
        }
        assertTrue("no budget ran out right after a hold", endedAtBest > 0);

        // the parsed settings are untouched: without a budget the run converges, with a held at its bound
        a.setParameterValue(0, 3.0);
        b.setParameterValue(0, 3.0);
        maximizer.maximize();
        assertEquals(Reason.CONVERGED, maximizer.getReason());
        assertFalse(maximizer.endedAtBestPoint());
        assertEquals(BOUND, a.getParameterValue(0), 0.0);
        assertEquals(BOUND, b.getParameterValue(0), 1E-5);
        assertEquals(1, maximizer.getHeld().size());
    }

    // the budget replaces the parsed stopping rules for one call only
    public void testBudgetIsPerCall() {
        MaximizerWrtParameter maximizer = maximizer(false, holdSettings(500));
        maximizer.maximize(new Budget(2).setQuiet(true));
        assertEquals(Reason.MAX_ITERATIONS, maximizer.getReason());
        assertEquals(2, maximizer.getIterations());

        maximizer.maximize(new Budget(0).setLimits(0.0, 3, 0.0, 0).setQuiet(true));
        assertEquals(Reason.MAX_EVALUATIONS, maximizer.getReason());

        maximizer.maximize();
        assertEquals(Reason.CONVERGED, maximizer.getReason());
        assertTrue(maximizer.describeSettings().contains("nIterations = 500"));
    }

    // a gradient that points uphill makes the first line search fail: the call ends at its start, which is its
    // best accepted point (lbfgs4j returns the last accepted point after a failed line search)
    public void testLineSearchFailureEndsAtBest() {
        MaximizerWrtParameter maximizer = maximizer(true, new MaximizerWrtParameter.Settings(100, 1E-6,
                true, false, false));
        maximizer.maximize(new Budget(100).setQuiet(true));
        assertEquals(Reason.LINE_SEARCH_FAILURE, maximizer.getReason());
        assertEquals(maximizer.getStartValue(), maximizer.getBestValue(), 0.0);
        assertEquals(maximizer.getBestValue(), maximizer.getMinimumValue(), 1E-12);
        assertEquals(3.0, a.getParameterValue(0), 0.0);
        assertEquals(3.0, b.getParameterValue(0), 0.0);
        assertFalse(maximizer.endedAtBestPoint());
    }

    // -------------------------------------------------------------------------------------------------------------

    private static double objective(double a, double b) {
        return 0.5 * ((a - 1) * (a - 1) + 100 * (b - a) * (b - a));
    }

    private MaximizerWrtParameter.Settings holdSettings(int iterations) {
        return new MaximizerWrtParameter.Settings(iterations, 1E-6, true, false, false)
                .addLowerBoundHold(a, BOUND).setConvergence(1E-8, null, null);
    }

    // N((1, 1), P) with P = [[101, -100], [-100, 100]], normalized so that F is exactly the quadratic above
    private MaximizerWrtParameter maximizer(boolean uphill, MaximizerWrtParameter.Settings settings) {
        final double[][] precision = {{101.0, -100.0}, {-100.0, 100.0}};
        final MultivariateNormalDistribution normal = new MultivariateNormalDistribution(new double[]{1.0, 1.0},
                precision);
        final double logNormalizer = normal.logPdf(new double[]{1.0, 1.0});
        final MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(normal) {
            @Override
            public double calculateLogLikelihood() {
                return super.calculateLogLikelihood() - logNormalizer;
            }
        };
        likelihood.addData(x);
        final GradientWrtParameterProvider wrapped =
                new GradientWrtParameterProvider.ParameterWrapper(normal, x, likelihood);
        GradientWrtParameterProvider gradient = new GradientWrtParameterProvider() {
            public Likelihood getLikelihood() { return wrapped.getLikelihood(); }
            public Parameter getParameter() { return wrapped.getParameter(); }
            public int getDimension() { return wrapped.getDimension(); }
            public double[] getGradientLogDensity() {
                double[] g = wrapped.getGradientLogDensity().clone();
                if (uphill) {
                    for (int i = 0; i < g.length; i++) {
                        g[i] = -g[i];
                    }
                }
                return g;
            }
        };
        return new MaximizerWrtParameter(likelihood, x, gradient, null, settings);
    }
}
