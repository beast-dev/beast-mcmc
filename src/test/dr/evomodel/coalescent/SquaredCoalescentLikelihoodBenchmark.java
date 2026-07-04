/*
 * SquaredCoalescentLikelihoodBenchmark.java
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
import dr.evomodel.coalescent.timeline.NoCovariateSegmentProvider;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evolution.tree.NodeRef;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.SquaredCachedSplines;
import junit.framework.TestCase;

import java.util.Collections;
import java.util.Random;

/**
 * Speed comparison between SquaredSplineCoalescentLikelihood (old, full-matrix storage)
 * and SquaredBasisCoalescentLikelihood (new, lower-triangle accumulator).
 *
 * Two scenarios are timed:
 *   COLD-G  : tree node height changes → G must be rebuilt on each call
 *   WARM-*  : only coefficients change → G is cached; measures logL / gradient / diag-Hessian
 *
 * Results are printed to stdout; the test always passes.
 */
public class SquaredCoalescentLikelihoodBenchmark extends TestCase {

    // 30-tip random coalescent tree (generated once; heights in [0, ~1.28])
    private static final String TREE_30 =
            "((t20:0.1441,t21:0.1441):1.1315,(t10:0.5954," +
            "((((t26:0.0052,t16:0.0052):0.2502,t29:0.2555):0.0283," +
            "((t30:0.0361,((t9:0.0020,t1:0.0020):0.0037,t22:0.0057):0.0304):0.1670," +
            "(((t17:0.0642,(t5:0.0160,(t2:0.0060,t27:0.0060):0.0100):0.0482):0.0240," +
            "(t28:0.0809,(t25:0.0235,t12:0.0235):0.0575):0.0072):0.0148," +
            "(t19:0.0114,(t15:0.0086,(t14:0.0075,t4:0.0075):0.0011):0.0028):0.0915):" +
            "0.1002):0.0806):0.2291,(((t6:0.0111,t18:0.0111):0.0274,t3:0.0385):0.1760," +
            "((t8:0.0424,(t24:0.0021,t7:0.0021):0.0403):0.1497," +
            "((t23:0.0171,t11:0.0171):0.0367,t13:0.0539):0.1381):0.0225):0.2983):" +
            "0.0825):0.6803);";

    // Larger basis for a more discriminating benchmark: 10 coefficients → 11×11 Gram matrix
    private static final double   LOWER   = 0.0;
    private static final double   UPPER   = 2.0;
    private static final int      DEGREE  = 3;
    private static final double[] KNOTS   = {0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75};
    // dim = KNOTS.length + DEGREE = 10

    private static final double   EPS     = 1e-4;
    private static final int      NWARMUP = 500;
    private static final int      NITER   = 3000;

    // -----------------------------------------------------------------------
    // Benchmark: all scenarios
    // -----------------------------------------------------------------------

    public void testBenchmark() throws Exception {
        System.out.println("\n=== SquaredCoalescentLikelihood Speed Benchmark ===");
        System.out.printf("Tree: 30 tips, %d coalescent events%n", 29);
        System.out.printf("Basis: degree=%d, %d interior knots → dim=%d (Gram: %dx%d)%n",
                DEGREE, KNOTS.length, KNOTS.length + DEGREE,
                KNOTS.length + DEGREE + 1, KNOTS.length + DEGREE + 1);
        System.out.printf("Iterations: %d warmup + %d timed%n%n", NWARMUP, NITER);

        // ---- build tree ----
        Tree rawTree = new NewickImporter(TREE_30).importTree(null);
        DefaultTreeModel treeModel = new DefaultTreeModel(rawTree);
        BigFastTreeIntervals intervals = new BigFastTreeIntervals(treeModel);

        // ---- shared parameters ----
        double[] initCoeffs = randomCoeffs(KNOTS.length + DEGREE, 42);
        Parameter intercept    = new Parameter.Default(0.4);
        Parameter coefficients = new Parameter.Default(initCoeffs);
        SquaredCachedSplines splines = new SquaredCachedSplines(
                coefficients, intercept, KNOTS, LOWER, UPPER, DEGREE);

        Parameter epsilon  = new Parameter.Default("epsilon",  EPS, 0.0, Double.MAX_VALUE);
        Parameter epsilon2 = new Parameter.Default("epsilon2", EPS, 0.0, Double.MAX_VALUE);
        Parameter ploidyOld = buildPloidy(1);
        Parameter ploidyNew = buildPloidy(1);

        // ---- OLD likelihood ----
        CompoundParameter oldRate = new CompoundParameter("oldRate");
        oldRate.addParameter(coefficients);
        oldRate.addParameter(intercept);
        SquaredSplineCoalescentLikelihood oldLk = new SquaredSplineCoalescentLikelihood(
                Collections.singletonList(intervals), splines, epsilon, oldRate, ploidyOld);

        // ---- NEW likelihood ----
        CompoundParameter newRate = new CompoundParameter("newRate");
        newRate.addParameter(intercept);
        newRate.addParameter(coefficients);
        BSplineBasisExpansion basis = new BSplineBasisExpansion(splines, newRate, true);
        SquaredBasisCoalescentLikelihood newLk = new SquaredBasisCoalescentLikelihood(
                Collections.singletonList(intervals), basis,
                NoCovariateSegmentProvider.INSTANCE, epsilon2, ploidyNew);

        // Sanity: values agree
        double oldL = oldLk.getLogLikelihood();
        double newL = newLk.getLogLikelihood();
        assertEquals("logL equivalence", oldL, newL, 1e-10);

        // Pick an internal node for tree-change benchmark
        NodeRef internalNode = treeModel.getInternalNode(0);
        double origHeight    = treeModel.getNodeHeight(internalNode);

        // ---- SCENARIO 1: COLD-G (tree changes → G rebuilt each call) ----
        System.out.println("--- COLD-G: log L after tree node-height change (G rebuilt) ---");
        benchmarkColdG("OLD logL", NWARMUP, NITER, oldLk, treeModel, internalNode, origHeight);
        benchmarkColdG("NEW logL", NWARMUP, NITER, newLk, treeModel, internalNode, origHeight);

        System.out.println();

        // ---- SCENARIO 2: WARM-G logL (only coefficients change) ----
        System.out.println("--- WARM-G: log L with cached G (coefficient-only change) ---");
        benchmarkWarmLogL("OLD logL", NWARMUP, NITER, oldLk, coefficients);
        benchmarkWarmLogL("NEW logL", NWARMUP, NITER, newLk, coefficients);

        System.out.println();

        // ---- SCENARIO 3: WARM-G gradient ----
        System.out.println("--- WARM-G: gradient with cached G (coefficient-only change) ---");
        benchmarkWarmGradient("OLD grad", NWARMUP, NITER, oldLk, coefficients);
        benchmarkWarmGradient("NEW grad", NWARMUP, NITER, newLk, coefficients);

        System.out.println();

        // ---- SCENARIO 4: WARM-G diagonal Hessian ----
        System.out.println("--- WARM-G: diagonal Hessian with cached G ---");
        benchmarkWarmDiagHessian("OLD diag-H", NWARMUP, NITER, oldLk, coefficients);
        benchmarkWarmDiagHessian("NEW diag-H", NWARMUP, NITER, newLk, coefficients);

        System.out.println("\n=== end benchmark ===\n");
    }

    // -----------------------------------------------------------------------
    // Scenario runners
    // -----------------------------------------------------------------------

    private void benchmarkColdG(String label, int nWarm, int nIter,
                                  Object lk,
                                  DefaultTreeModel tree, NodeRef node, double origH)
            throws Exception {
        double delta = 1e-8;
        for (int i = 0; i < nWarm; i++) {
            tree.setNodeHeight(node, origH + delta);
            getLogLikelihood(lk);
            tree.setNodeHeight(node, origH);
            getLogLikelihood(lk);
        }
        long t0 = System.nanoTime();
        for (int i = 0; i < nIter; i++) {
            tree.setNodeHeight(node, origH + delta);
            getLogLikelihood(lk);
            tree.setNodeHeight(node, origH);
            getLogLikelihood(lk);
        }
        long elapsed = System.nanoTime() - t0;
        double nsPerCall = (double) elapsed / (2.0 * nIter);
        System.out.printf("  %-12s : %7.1f ns/call%n", label, nsPerCall);
    }

    private void benchmarkWarmLogL(String label, int nWarm, int nIter,
                                    Object lk, Parameter coefficients) {
        double orig = coefficients.getParameterValue(0);
        double delta = 1e-8;
        // Prime G cache
        getLogLikelihood(lk);
        for (int i = 0; i < nWarm; i++) {
            coefficients.setParameterValue(0, orig + delta);
            getLogLikelihood(lk);
            coefficients.setParameterValue(0, orig);
            getLogLikelihood(lk);
        }
        long t0 = System.nanoTime();
        for (int i = 0; i < nIter; i++) {
            coefficients.setParameterValue(0, orig + delta);
            getLogLikelihood(lk);
            coefficients.setParameterValue(0, orig);
            getLogLikelihood(lk);
        }
        long elapsed = System.nanoTime() - t0;
        double nsPerCall = (double) elapsed / (2.0 * nIter);
        System.out.printf("  %-12s : %7.1f ns/call%n", label, nsPerCall);
    }

    private void benchmarkWarmGradient(String label, int nWarm, int nIter,
                                        Object lk, Parameter coefficients) {
        double orig = coefficients.getParameterValue(0);
        double delta = 1e-8;
        getLogLikelihood(lk);
        for (int i = 0; i < nWarm; i++) {
            coefficients.setParameterValue(0, orig + delta);
            getGradient(lk);
            coefficients.setParameterValue(0, orig);
            getGradient(lk);
        }
        long t0 = System.nanoTime();
        for (int i = 0; i < nIter; i++) {
            coefficients.setParameterValue(0, orig + delta);
            getGradient(lk);
            coefficients.setParameterValue(0, orig);
            getGradient(lk);
        }
        long elapsed = System.nanoTime() - t0;
        double nsPerCall = (double) elapsed / (2.0 * nIter);
        System.out.printf("  %-12s : %7.1f ns/call%n", label, nsPerCall);
    }

    private void benchmarkWarmDiagHessian(String label, int nWarm, int nIter,
                                           Object lk, Parameter coefficients) {
        double orig = coefficients.getParameterValue(0);
        double delta = 1e-8;
        getLogLikelihood(lk);
        for (int i = 0; i < nWarm; i++) {
            coefficients.setParameterValue(0, orig + delta);
            getDiagHessian(lk);
            coefficients.setParameterValue(0, orig);
            getDiagHessian(lk);
        }
        long t0 = System.nanoTime();
        for (int i = 0; i < nIter; i++) {
            coefficients.setParameterValue(0, orig + delta);
            getDiagHessian(lk);
            coefficients.setParameterValue(0, orig);
            getDiagHessian(lk);
        }
        long elapsed = System.nanoTime() - t0;
        double nsPerCall = (double) elapsed / (2.0 * nIter);
        System.out.printf("  %-12s : %7.1f ns/call%n", label, nsPerCall);
    }

    // -----------------------------------------------------------------------
    // Dispatch helpers (both old and new implement the same methods)
    // -----------------------------------------------------------------------

    private static double getLogLikelihood(Object lk) {
        if (lk instanceof SquaredSplineCoalescentLikelihood)
            return ((SquaredSplineCoalescentLikelihood) lk).getLogLikelihood();
        return ((SquaredBasisCoalescentLikelihood) lk).getLogLikelihood();
    }

    private static double[] getGradient(Object lk) {
        if (lk instanceof SquaredSplineCoalescentLikelihood)
            return ((SquaredSplineCoalescentLikelihood) lk).getGradientLogDensity();
        return ((SquaredBasisCoalescentLikelihood) lk).getGradientLogDensity();
    }

    private static double[] getDiagHessian(Object lk) {
        if (lk instanceof SquaredSplineCoalescentLikelihood)
            return ((SquaredSplineCoalescentLikelihood) lk).getDiagonalHessianLogDensity();
        return ((SquaredBasisCoalescentLikelihood) lk).getDiagonalHessianLogDensity();
    }

    // -----------------------------------------------------------------------
    // Utilities
    // -----------------------------------------------------------------------

    private static double[] randomCoeffs(int dim, long seed) {
        Random rng = new Random(seed);
        double[] v = new double[dim];
        for (int i = 0; i < dim; i++) v[i] = rng.nextGaussian() * 0.2;
        return v;
    }

    private static Parameter buildPloidy(int n) {
        Parameter p = new Parameter.Default(n);
        for (int i = 0; i < n; i++) p.setParameterValue(i, 1.0);
        return p;
    }
}
