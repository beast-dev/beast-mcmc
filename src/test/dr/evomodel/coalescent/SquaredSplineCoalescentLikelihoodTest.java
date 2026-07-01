/*
 * SquaredSplineCoalescentLikelihoodTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package test.dr.evomodel.coalescent;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.SquaredSplineCoalescentLikelihood;
import dr.evomodel.tree.DefaultTreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.SquaredCachedSplines;
import junit.framework.TestCase;

import java.util.Collections;

/**
 * Tests for SquaredSplineCoalescentLikelihood:
 *
 *  (1) Gradient finite-difference check on rateParameter = [coefficients, intercept]
 *  (2) Gradient finite-difference check on epsilon
 *  (3) Multilocus: two identical trees give log-likelihood exactly double a single tree
 */
public class SquaredSplineCoalescentLikelihoodTest extends TestCase {

    // A small isochronous 5-tip tree (internal times well-separated for a clean test)
    private static final String FIVE_TIP_TREE =
            "((t1:0.5,t2:0.5):1.0,(t3:0.8,(t4:0.3,t5:0.3):0.5):0.7);";

    private static final double[] INTERIOR_KNOTS = {1.0, 2.0};
    private static final double   LOWER          = 0.0;
    private static final double   UPPER          = 4.0;
    private static final int      DEGREE         = 3;

    // 5 spline coefficients (= INTERIOR_KNOTS.length + DEGREE)
    private static final double[] THETA     = {0.3, -0.1, 0.2, 0.05, -0.15};
    private static final double   INTERCEPT = 0.4;
    private static final double   EPS       = 1e-4;

    private static final double FD_H       = 1e-6;
    private static final double TOLERANCE  = 2e-6;

    // -----------------------------------------------------------------------
    // Tests
    // -----------------------------------------------------------------------

    public void testRateParameterGradientMatchesFiniteDifference() throws Exception {
        SquaredSplineCoalescentLikelihood likelihood = buildLikelihood(FIVE_TIP_TREE, 1);

        Parameter rateParam = likelihood.getRateParameter();
        double[] analytic   = likelihood.getGradientLogDensity();
        double[] numeric    = finiteDifferenceGradient(likelihood, rateParam);

        assertEquals("gradient dimension", numeric.length, analytic.length);
        for (int i = 0; i < analytic.length; i++) {
            assertEquals("gradient[" + i + "]", numeric[i], analytic[i], TOLERANCE);
        }
    }

    public void testEpsilonGradientMatchesFiniteDifference() throws Exception {
        SquaredSplineCoalescentLikelihood likelihood = buildLikelihood(FIVE_TIP_TREE, 1);

        Parameter eps    = likelihood.getEpsilon();
        double analytic  = epsilonDerivativeAnalytic(likelihood);
        double numeric   = finiteDifferenceScalar(likelihood, eps, 0);

        assertEquals("epsilon gradient", numeric, analytic, TOLERANCE);
    }

    public void testMultilocusLogLikelihoodIsDouble() throws Exception {
        SquaredSplineCoalescentLikelihood single = buildLikelihood(FIVE_TIP_TREE, 1);
        SquaredSplineCoalescentLikelihood multi  = buildLikelihood(FIVE_TIP_TREE, 2);

        double singleLogL = single.getLogLikelihood();
        double multiLogL  = multi.getLogLikelihood();

        assertEquals("multilocus logL = 2 * single logL",
                2.0 * singleLogL, multiLogL, Math.abs(singleLogL) * 1e-10);
    }

    public void testLikelihoodCaching() throws Exception {
        SquaredSplineCoalescentLikelihood likelihood = buildLikelihood(FIVE_TIP_TREE, 1);

        double first  = likelihood.getLogLikelihood();
        double second = likelihood.getLogLikelihood();
        assertEquals("cached logL equals recomputed", first, second, 0.0);
    }

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    private static SquaredSplineCoalescentLikelihood buildLikelihood(String newick,
                                                                      int nTrees)
            throws Exception {
        Tree rawTree = new NewickImporter(newick).importTree(null);
        DefaultTreeModel treeModel = new DefaultTreeModel(rawTree);
        BigFastTreeIntervals intervals = new BigFastTreeIntervals(treeModel);

        Parameter coefficients = new Parameter.Default(THETA);
        Parameter intercept    = new Parameter.Default(INTERCEPT);
        SquaredCachedSplines splines = new SquaredCachedSplines(
                coefficients, intercept, INTERIOR_KNOTS, LOWER, UPPER, DEGREE);

        // The rate parameter is created by the caller and shared with operators and
        // gradient providers — it must NOT be created anonymously inside the likelihood.
        CompoundParameter rateParameter = new CompoundParameter("squaredSplineRate");
        rateParameter.addParameter(coefficients);
        rateParameter.addParameter(intercept);

        Parameter epsilon = new Parameter.Default("epsilon", EPS, 0.0, Double.MAX_VALUE);

        Parameter ploidyFactors = new Parameter.Default(nTrees);
        for (int i = 0; i < nTrees; i++) ploidyFactors.setParameterValue(i, 1.0);

        java.util.List<BigFastTreeIntervals> intervalsList = new java.util.ArrayList<>();
        for (int i = 0; i < nTrees; i++) intervalsList.add(intervals);

        return new SquaredSplineCoalescentLikelihood(
                intervalsList, splines, epsilon, rateParameter, ploidyFactors);
    }

    /** Finite-difference gradient w.r.t. every component of param. */
    private static double[] finiteDifferenceGradient(SquaredSplineCoalescentLikelihood likelihood,
                                                      Parameter param) {
        int dim  = param.getDimension();
        double[] grad = new double[dim];
        for (int i = 0; i < dim; i++) {
            grad[i] = finiteDifferenceScalar(likelihood, param, i);
        }
        return grad;
    }

    private static double finiteDifferenceScalar(SquaredSplineCoalescentLikelihood likelihood,
                                                  Parameter param, int index) {
        double orig = param.getParameterValue(index);

        param.setParameterValue(index, orig + FD_H);
        double plus  = likelihood.getLogLikelihood();

        param.setParameterValue(index, orig - FD_H);
        double minus = likelihood.getLogLikelihood();

        param.setParameterValue(index, orig);
        return (plus - minus) / (2.0 * FD_H);
    }

    /**
     * Analytic derivative of log L w.r.t. epsilon.
     *
     *   d/dε log L = sum_intervals [ -C(n,2)*ρ*(end-start) ]
     *              + sum_coal_events [ 1 / (ε + f(t)²) ]
     */
    private static double epsilonDerivativeAnalytic(SquaredSplineCoalescentLikelihood likelihood) {
        BigFastTreeIntervals intervals = likelihood.getIntervalsList().get(0);
        double ploidy   = likelihood.getPloidyFactors().getParameterValue(0);
        double epsValue = likelihood.getEpsilon().getParameterValue(0);
        SquaredCachedSplines splines = likelihood.getSplines();

        double deriv = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            double start = intervals.getIntervalTime(i);
            double end   = start + intervals.getInterval(i);
            int    n     = intervals.getLineageCount(i);
            deriv += -0.5 * n * (n - 1) * ploidy * (end - start);
            if (intervals.getIntervalType(i) == dr.evolution.coalescent.IntervalType.COALESCENT) {
                double f    = splines.evaluateSpline(end);
                double rate = epsValue + f * f;
                deriv += 1.0 / rate;
            }
        }
        return deriv;
    }
}
