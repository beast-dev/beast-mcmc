/*
 * SquaredBasisCoalescentLikelihoodTest.java
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
import dr.evomodel.coalescent.SquaredBasisCoalescentLikelihood;
import dr.evomodel.coalescent.SquaredSplineCoalescentLikelihood;
import dr.evomodel.coalescent.basis.BSplineBasisExpansion;
import dr.evomodel.coalescent.timeline.BasisContext;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.evomodel.coalescent.basis.CovariateAugmentedBasisExpansion;
import dr.evomodel.coalescent.basis.CovariateMode;
import dr.evomodel.coalescent.timeline.NoCovariateSegmentProvider;
import dr.evomodel.coalescent.timeline.PiecewiseConstantCovariateSegmentProvider;
import dr.evomodel.tree.DefaultTreeModel;
import dr.inference.model.Parameter;
import dr.math.SquaredCachedSplines;
import junit.framework.TestCase;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Tests for SquaredBasisCoalescentLikelihood.
 *
 * Regression values are pinned against the first correct run.
 */
public class SquaredBasisCoalescentLikelihoodTest extends TestCase {

    // Small isochronous 5-tip tree
    private static final String FIVE_TIP_TREE =
            "((t1:0.5,t2:0.5):1.0,(t3:0.8,(t4:0.3,t5:0.3):0.5):0.7);";

    private static final double[] INTERIOR_KNOTS = {1.0, 2.0};
    private static final double   LOWER          = 0.0;
    private static final double   UPPER          = 4.0;
    private static final int      DEGREE         = 3;

    // 5 coefficients (= INTERIOR_KNOTS.length + DEGREE) — canonical layout: intercept first
    private static final double   INTERCEPT      = 0.4;
    private static final double[] THETA          = {0.3, -0.1, 0.2, 0.05, -0.15};
    private static final double   EPS            = 1e-4;

    private static final double FD_H             = 1e-6;
    private static final double GRADIENT_TOL     = 2e-6;
    private static final double LOGLIK_TOL       = 1e-12;

    // -----------------------------------------------------------------------
    // Tests: no-covariate spline with intercept
    // -----------------------------------------------------------------------

    public void testGradientMatchesFiniteDifference() throws Exception {
        SquaredBasisCoalescentLikelihood lk = buildLikelihood(FIVE_TIP_TREE, 1);

        Parameter param   = lk.getParameter();
        double[] analytic = lk.getGradientLogDensity();
        double[] numeric  = finiteDifferenceGradient(lk, param);

        assertEquals("gradient dimension", numeric.length, analytic.length);
        for (int i = 0; i < analytic.length; i++) {
            assertEquals("gradient[" + i + "]", numeric[i], analytic[i], GRADIENT_TOL);
        }
    }

    public void testEpsilonGradientMatchesFiniteDifference() throws Exception {
        SquaredBasisCoalescentLikelihood lk = buildLikelihood(FIVE_TIP_TREE, 1);
        Parameter eps = lk.getEpsilon();
        double numeric  = finiteDifferenceScalar(lk, eps, 0);
        double analytic = epsilonDerivativeAnalytic(lk);
        assertEquals("epsilon gradient", numeric, analytic, GRADIENT_TOL);
    }

    public void testMultilocusLogLikelihoodIsDouble() throws Exception {
        SquaredBasisCoalescentLikelihood single = buildLikelihood(FIVE_TIP_TREE, 1);
        SquaredBasisCoalescentLikelihood multi  = buildLikelihood(FIVE_TIP_TREE, 2);

        double singleL = single.getLogLikelihood();
        double multiL  = multi.getLogLikelihood();
        assertEquals("multilocus logL = 2 * single logL",
                2.0 * singleL, multiL, Math.abs(singleL) * 1e-10);
    }

    public void testLikelihoodCaching() throws Exception {
        SquaredBasisCoalescentLikelihood lk = buildLikelihood(FIVE_TIP_TREE, 1);
        double first  = lk.getLogLikelihood();
        double second = lk.getLogLikelihood();
        assertEquals("cached logL equals recomputed", first, second, 0.0);
    }

    public void testHessianDiagonalConsistentWithGradientFD() throws Exception {
        SquaredBasisCoalescentLikelihood lk = buildLikelihood(FIVE_TIP_TREE, 1);
        Parameter param    = lk.getParameter();
        double[] hessianFD = finiteDifferenceHessianDiagonal(lk, param);
        double[] hessianA  = lk.getDiagonalHessianLogDensity();
        for (int i = 0; i < hessianA.length; i++) {
            assertEquals("hessian diag[" + i + "]", hessianFD[i], hessianA[i], 1e-5);
        }
    }

    // -----------------------------------------------------------------------
    // Tests: covariate augmented (ADDITIVE mode)
    // -----------------------------------------------------------------------

    /**
     * Directly verifies that CovariateAugmentedBasisExpansion.evaluateScalar(t, ctx)
     * matches dot(gamma, evaluateBasis(t, ctx)) at several time points across different
     * covariate segments.  This is the scalar shortcut used by evaluateRateAt in the
     * logL hot path after the evaluateScalar optimization.
     */
    public void testCovariateEvaluateScalarMatchesDotProduct() throws Exception {
        // Build the augmented basis directly (reuse same setup as buildCovariateLikelihood)
        Parameter intercept    = new Parameter.Default(INTERCEPT);
        Parameter coefficients = new Parameter.Default(THETA);
        SquaredCachedSplines splines = new SquaredCachedSplines(
                coefficients, intercept, INTERIOR_KNOTS, LOWER, UPPER, DEGREE);
        Parameter rateParameter = buildCanonicalRateParameter(intercept, coefficients);

        double[] breaks  = {0.0, 0.5, 1.5, UPPER + 1.0};
        double[][] covs  = {{1.0}, {2.0}, {0.5}};
        Parameter covCoeff = new Parameter.Default(new double[]{0.1});

        BSplineBasisExpansion base = new BSplineBasisExpansion(splines, rateParameter, true);
        PiecewiseConstantCovariateSegmentProvider segProvider =
                new PiecewiseConstantCovariateSegmentProvider(breaks, covs);
        CovariateAugmentedBasisExpansion augBasis =
                new CovariateAugmentedBasisExpansion(base, covCoeff, 1, CovariateMode.ADDITIVE);

        int dim = augBasis.getDimension();
        // gamma = [intercept, theta_0..theta_4, covCoeff_0] from the compound parameter
        Parameter fullParam = augBasis.getParameter();
        double[] gamma = new double[dim];
        for (int i = 0; i < dim; i++) gamma[i] = fullParam.getParameterValue(i);

        double[] basisBuf = new double[dim];

        // Test at several time points spanning all three covariate segments
        double[] testTimes = {0.1, 0.4, 0.5, 0.8, 1.2, 1.5, 1.8, 2.5};
        for (double t : testTimes) {
            BasisContext ctx = segProvider.getEventContext(0, t);
            double scalar   = augBasis.evaluateScalar(t, ctx);
            augBasis.evaluateBasis(t, ctx, basisBuf);
            double dotProd  = dot(gamma, basisBuf);
            assertEquals("evaluateScalar vs dot at t=" + t, dotProd, scalar, 1e-14);
        }
    }

    /**
     * Replicates what HMC's checkGradient does for the covariate likelihood:
     * uses NumericalDerivative.gradient() (which sets the full parameter vector
     * at once, exactly as a leapfrog step does) and compares against the analytic
     * gradient at several random parameter configurations.
     *
     * This exercises the full event chain:
     *   setParameterValue → handleVariableChangedEvent → gammaBufferKnown=false
     *   → getLogLikelihood() → evaluateRateAt → evaluateScalar → correct lambda(t)
     */
    public void testCovariateHmcStyleGradientCheck() throws Exception {
        final SquaredBasisCoalescentLikelihood lk = buildCovariateLikelihood(FIVE_TIP_TREE);
        final Parameter param = lk.getParameter();
        final int dim = param.getDimension();

        MultivariateFunction logLFunc = new MultivariateFunction() {
            @Override
            public double evaluate(double[] x) {
                for (int i = 0; i < dim; i++) param.setParameterValue(i, x[i]);
                return lk.getLogLikelihood();
            }
            @Override public int getNumArguments() { return dim; }
            @Override public double getLowerBound(int n) { return Double.NEGATIVE_INFINITY; }
            @Override public double getUpperBound(int n) { return Double.POSITIVE_INFINITY; }
        };

        // Test at 4 configurations: initial + 3 perturbed (fixed seed for reproducibility)
        java.util.Random rng = new java.util.Random(12345);
        double[] base = new double[dim];
        for (int i = 0; i < dim; i++) base[i] = param.getParameterValue(i);

        for (int trial = 0; trial < 4; trial++) {
            double[] x = base.clone();
            if (trial > 0) {
                for (int i = 0; i < dim; i++) x[i] += rng.nextGaussian() * 0.3;
            }
            for (int i = 0; i < dim; i++) param.setParameterValue(i, x[i]);

            double[] analytic = lk.getGradientLogDensity();
            double[] numeric  = NumericalDerivative.gradient(logLFunc, x);

            // Restore
            for (int i = 0; i < dim; i++) param.setParameterValue(i, x[i]);

            for (int i = 0; i < dim; i++) {
                assertEquals("HMC-style covariate grad[" + i + "] trial=" + trial,
                        numeric[i], analytic[i], GRADIENT_TOL);
            }
        }

        // Restore original values
        for (int i = 0; i < dim; i++) param.setParameterValue(i, base[i]);
    }

    public void testCovariateGradientMatchesFiniteDifference() throws Exception {
        SquaredBasisCoalescentLikelihood lk = buildCovariateLikelihood(FIVE_TIP_TREE);

        Parameter param   = lk.getParameter();
        double[] analytic = lk.getGradientLogDensity();
        double[] numeric  = finiteDifferenceGradient(lk, param);

        assertEquals("covariate gradient dimension", numeric.length, analytic.length);
        for (int i = 0; i < analytic.length; i++) {
            assertEquals("covariate gradient[" + i + "]", numeric[i], analytic[i], GRADIENT_TOL);
        }
    }

    public void testCovariateEpsilonGradientMatchesFiniteDifference() throws Exception {
        SquaredBasisCoalescentLikelihood lk = buildCovariateLikelihood(FIVE_TIP_TREE);
        Parameter eps = lk.getEpsilon();
        double numeric  = finiteDifferenceScalar(lk, eps, 0);
        double analytic = epsilonDerivativeAnalytic(lk);
        assertEquals("covariate epsilon gradient", numeric, analytic, GRADIENT_TOL);
    }

    // -----------------------------------------------------------------------
    // Tests: numerical equivalence with legacy SquaredSplineCoalescentLikelihood
    // -----------------------------------------------------------------------

    /**
     * Builds both the old and new likelihood on the same tree and parameter values
     * and checks that log L, gradient, and diagonal Hessian agree up to floating-point
     * tolerance.
     *
     * Layout mapping (p = coefficients.getDimension() = 5):
     *   old rateParam = [θ_0,...,θ_{p-1}, intercept]
     *   new rateParam = [intercept, θ_0,...,θ_{p-1}]  (canonical)
     *
     *   old_grad[i]   = d/dθ_i         for i = 0..p-1
     *   old_grad[p]   = d/d_intercept
     *   new_grad[0]   = d/d_intercept
     *   new_grad[i+1] = d/dθ_i         for i = 0..p-1
     */
    public void testEquivalenceWithLegacyClass() throws Exception {
        Parameter intercept    = new Parameter.Default(INTERCEPT);
        Parameter coefficients = new Parameter.Default(THETA);
        SquaredCachedSplines splines = new SquaredCachedSplines(
                coefficients, intercept, INTERIOR_KNOTS, LOWER, UPPER, DEGREE);

        BigFastTreeIntervals intervals = buildIntervals(FIVE_TIP_TREE);
        int p = coefficients.getDimension();

        // Old layout: [coefficients, intercept]
        dr.inference.model.CompoundParameter oldRateParam =
                new dr.inference.model.CompoundParameter("oldRate");
        oldRateParam.addParameter(coefficients);
        oldRateParam.addParameter(intercept);
        Parameter epsilonOld = new Parameter.Default("epsilon", EPS, 0.0, Double.MAX_VALUE);
        Parameter ploidyOld  = buildPloidy(1);
        SquaredSplineCoalescentLikelihood oldLk = new SquaredSplineCoalescentLikelihood(
                Collections.singletonList(intervals), splines, epsilonOld, oldRateParam, ploidyOld);

        // New canonical layout: [intercept, coefficients]
        Parameter canonicalRate = buildCanonicalRateParameter(intercept, coefficients);
        BSplineBasisExpansion basis = new BSplineBasisExpansion(splines, canonicalRate, true);
        Parameter epsilonNew = new Parameter.Default("epsilon2", EPS, 0.0, Double.MAX_VALUE);
        Parameter ploidyNew  = buildPloidy(1);
        SquaredBasisCoalescentLikelihood newLk = new SquaredBasisCoalescentLikelihood(
                Collections.singletonList(intervals), basis,
                NoCovariateSegmentProvider.INSTANCE, epsilonNew, ploidyNew);

        // Log-likelihood must be identical
        double oldLogL = oldLk.getLogLikelihood();
        double newLogL = newLk.getLogLikelihood();
        assertEquals("logL equivalence", oldLogL, newLogL, LOGLIK_TOL);

        // Gradients must match after layout reordering
        double[] oldGrad = oldLk.getGradientLogDensity();
        double[] newGrad = newLk.getGradientLogDensity();
        assertEquals("gradient dimension", oldGrad.length, newGrad.length);
        assertEquals("grad[intercept]", oldGrad[p], newGrad[0], GRADIENT_TOL);
        for (int i = 0; i < p; i++) {
            assertEquals("grad[theta_" + i + "]", oldGrad[i], newGrad[i + 1], GRADIENT_TOL);
        }

        // Diagonal Hessian must match after layout reordering
        double[] oldDiag = oldLk.getDiagonalHessianLogDensity();
        double[] newDiag = newLk.getDiagonalHessianLogDensity();
        assertEquals("diag-hessian dimension", oldDiag.length, newDiag.length);
        assertEquals("diag-hessian[intercept]", oldDiag[p], newDiag[0], 1e-5);
        for (int i = 0; i < p; i++) {
            assertEquals("diag-hessian[theta_" + i + "]", oldDiag[i], newDiag[i + 1], 1e-5);
        }
    }

    // -----------------------------------------------------------------------
    // Tests: no-covariate spline without intercept
    // -----------------------------------------------------------------------

    public void testNoInterceptGradientMatchesFiniteDifference() throws Exception {
        SquaredBasisCoalescentLikelihood lk = buildLikelihoodNoIntercept(FIVE_TIP_TREE);

        Parameter param   = lk.getParameter();
        double[] analytic = lk.getGradientLogDensity();
        double[] numeric  = finiteDifferenceGradient(lk, param);

        assertEquals("no-intercept gradient dimension", numeric.length, analytic.length);
        for (int i = 0; i < analytic.length; i++) {
            assertEquals("no-intercept gradient[" + i + "]", numeric[i], analytic[i], GRADIENT_TOL);
        }
    }

    // -----------------------------------------------------------------------
    // Builders
    // -----------------------------------------------------------------------

    private static SquaredBasisCoalescentLikelihood buildLikelihood(String newick, int nTrees)
            throws Exception {
        BigFastTreeIntervals intervals = buildIntervals(newick);

        // Canonical layout: intercept first
        Parameter intercept    = new Parameter.Default(INTERCEPT);
        Parameter coefficients = new Parameter.Default(THETA);
        SquaredCachedSplines splines = new SquaredCachedSplines(
                coefficients, intercept, INTERIOR_KNOTS, LOWER, UPPER, DEGREE);

        // rateParameter in canonical order [intercept, theta...]
        Parameter rateParameter = buildCanonicalRateParameter(intercept, coefficients);
        Parameter epsilon       = new Parameter.Default("epsilon", EPS, 0.0, Double.MAX_VALUE);
        Parameter ploidy        = buildPloidy(nTrees);

        List<BigFastTreeIntervals> list = Collections.nCopies(nTrees, intervals);

        BSplineBasisExpansion basis = new BSplineBasisExpansion(splines, rateParameter, true);

        return new SquaredBasisCoalescentLikelihood(
                list, basis, NoCovariateSegmentProvider.INSTANCE, epsilon, ploidy);
    }

    private static SquaredBasisCoalescentLikelihood buildLikelihoodNoIntercept(String newick)
            throws Exception {
        BigFastTreeIntervals intervals = buildIntervals(newick);

        Parameter coefficients = new Parameter.Default(THETA);
        Parameter intercept    = new Parameter.Default(INTERCEPT);
        SquaredCachedSplines splines = new SquaredCachedSplines(
                coefficients, intercept, INTERIOR_KNOTS, LOWER, UPPER, DEGREE);

        Parameter epsilon = new Parameter.Default("epsilon", EPS, 0.0, Double.MAX_VALUE);
        Parameter ploidy  = buildPloidy(1);

        BSplineBasisExpansion basis = new BSplineBasisExpansion(splines, coefficients, false);

        return new SquaredBasisCoalescentLikelihood(
                Collections.singletonList(intervals),
                basis, NoCovariateSegmentProvider.INSTANCE, epsilon, ploidy);
    }

    /**
     * Builds a likelihood with additive covariates.
     * Breakpoints at t=0.5 and t=1.5; covariate = [1.0] on [0,0.5), [0.5,1.5), [1.5,4.0).
     */
    private static SquaredBasisCoalescentLikelihood buildCovariateLikelihood(String newick)
            throws Exception {
        BigFastTreeIntervals intervals = buildIntervals(newick);

        Parameter intercept    = new Parameter.Default(INTERCEPT);
        Parameter coefficients = new Parameter.Default(THETA);
        SquaredCachedSplines splines = new SquaredCachedSplines(
                coefficients, intercept, INTERIOR_KNOTS, LOWER, UPPER, DEGREE);

        Parameter rateParameter = buildCanonicalRateParameter(intercept, coefficients);
        Parameter epsilon       = new Parameter.Default("epsilon", EPS, 0.0, Double.MAX_VALUE);
        Parameter ploidy        = buildPloidy(1);

        BSplineBasisExpansion base = new BSplineBasisExpansion(splines, rateParameter, true);

        // Single covariate, 3 segments
        double[] breaks     = {0.0, 0.5, 1.5, UPPER + 1.0};
        double[][] covs     = {{1.0}, {2.0}, {0.5}};
        Parameter covCoeff  = new Parameter.Default(new double[]{0.1}); // covariate coefficient

        PiecewiseConstantCovariateSegmentProvider segProvider =
                new PiecewiseConstantCovariateSegmentProvider(breaks, covs);

        CovariateAugmentedBasisExpansion augBasis =
                new CovariateAugmentedBasisExpansion(base, covCoeff, 1, CovariateMode.ADDITIVE);

        return new SquaredBasisCoalescentLikelihood(
                Collections.singletonList(intervals), augBasis, segProvider, epsilon, ploidy);
    }

    // -----------------------------------------------------------------------
    // Helper utilities
    // -----------------------------------------------------------------------

    private static BigFastTreeIntervals buildIntervals(String newick) throws Exception {
        Tree rawTree = new NewickImporter(newick).importTree(null);
        return new BigFastTreeIntervals(new DefaultTreeModel(rawTree));
    }

    /** Builds a CompoundParameter with [intercept, coefficients...] in canonical order. */
    private static Parameter buildCanonicalRateParameter(Parameter intercept,
                                                          Parameter coefficients) {
        dr.inference.model.CompoundParameter cp =
                new dr.inference.model.CompoundParameter("canonicalRate");
        cp.addParameter(intercept);
        cp.addParameter(coefficients);
        return cp;
    }

    private static Parameter buildPloidy(int nTrees) {
        Parameter ploidy = new Parameter.Default(nTrees);
        for (int i = 0; i < nTrees; i++) ploidy.setParameterValue(i, 1.0);
        return ploidy;
    }

    private static double[] finiteDifferenceGradient(SquaredBasisCoalescentLikelihood lk,
                                                      Parameter param) {
        int dim = param.getDimension();
        double[] grad = new double[dim];
        for (int i = 0; i < dim; i++) {
            grad[i] = finiteDifferenceScalar(lk, param, i);
        }
        return grad;
    }

    private static double finiteDifferenceScalar(SquaredBasisCoalescentLikelihood lk,
                                                  Parameter param, int index) {
        double orig = param.getParameterValue(index);
        param.setParameterValue(index, orig + FD_H);
        double plus  = lk.getLogLikelihood();
        param.setParameterValue(index, orig - FD_H);
        double minus = lk.getLogLikelihood();
        param.setParameterValue(index, orig);
        return (plus - minus) / (2.0 * FD_H);
    }

    private static double[] finiteDifferenceHessianDiagonal(SquaredBasisCoalescentLikelihood lk,
                                                              Parameter param) {
        // Use larger h for second-order FD to avoid floating-point cancellation.
        double h2 = 1e-4;
        int dim = param.getDimension();
        double[] diag = new double[dim];
        double orig0 = lk.getLogLikelihood();
        for (int i = 0; i < dim; i++) {
            double orig = param.getParameterValue(i);
            param.setParameterValue(i, orig + h2);
            double plus  = lk.getLogLikelihood();
            param.setParameterValue(i, orig - h2);
            double minus = lk.getLogLikelihood();
            param.setParameterValue(i, orig);
            diag[i] = (plus - 2.0 * orig0 + minus) / (h2 * h2);
        }
        return diag;
    }

    /**
     * Analytic d/dε log L = sum_intervals [-C(n,2)*rho*(end-start)]
     *                       + sum_coal_events [1 / (ε + eta^2)]
     */
    private static double epsilonDerivativeAnalytic(SquaredBasisCoalescentLikelihood lk) {
        BigFastTreeIntervals intervals = lk.getIntervalsList().get(0);
        double ploidy   = lk.getPloidyFactors().getParameterValue(0);
        double eps      = lk.getEpsilon().getParameterValue(0);

        Parameter param  = lk.getParameter();
        double[] gamma   = new double[param.getDimension()];
        for (int i = 0; i < gamma.length; i++) gamma[i] = param.getParameterValue(i);

        double[] buffer = new double[param.getDimension()];
        double deriv = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            double start = intervals.getIntervalTime(i);
            double end   = start + intervals.getInterval(i);
            int    n     = intervals.getLineageCount(i);
            deriv += -0.5 * n * (n - 1) * ploidy * (end - start);
            if (intervals.getIntervalType(i) == dr.evolution.coalescent.IntervalType.COALESCENT) {
                dr.evomodel.coalescent.timeline.BasisContext ctx =
                        lk.getSegmentProvider().getEventContext(0, end);
                lk.getBasis().evaluateBasis(end, ctx, buffer);
                double eta = dot(gamma, buffer);
                deriv += 1.0 / (eps + eta * eta);
            }
        }
        return deriv;
    }

    private static double dot(double[] a, double[] b) {
        double s = 0.0;
        for (int i = 0; i < a.length; i++) s += a[i] * b[i];
        return s;
    }
}
