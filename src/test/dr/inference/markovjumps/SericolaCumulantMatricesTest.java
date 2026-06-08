/*
 * SericolaCumulantMatricesTest.java
 *
 * Copyright (c) 2002-2024 the BEAST Development Team
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.markovjumps;

import junit.framework.TestCase;

public class SericolaCumulantMatricesTest extends TestCase {

    private static final int DIM = 3;
    private static final int DIM2 = DIM * DIM;
    private static final int PHI = DIM - 1;
    private static final double TOL = 1.0e-14;

    private static final double[] SORTED_ALPHA = {
            0.0, 0.35, 1.0
    };

    private static final double[] TRANSITION_MATRIX = {
            0.82, 0.11, 0.07,
            0.08, 0.83, 0.09,
            0.04, 0.19, 0.77
    };
    private static final double[] INV_ALPHA_DIFF = inverseAlphaDiff();
    private static final int[] OUT_ROW_BASE_BY_SORTED = {0, DIM, 2 * DIM};
    private static final int[] OUT_COL_BY_SORTED = {0, 1, 2};
    private static final double[] UPSTREAM_DENSITY_ADJOINT = {
            0.70, -0.20, 0.30,
            0.15, 0.90, -0.45,
            -0.10, 0.25, 0.55
    };

    public void testGrowthRecomputesConsistentCumulants() {
        SericolaCumulantMatrices cache = new SericolaCumulantMatrices(DIM);
        cache.ensureForTime(1.0, 2, TRANSITION_MATRIX, SORTED_ALPHA);

        double[] beforeGrowth = snapshot(cache, 2);
        double[] oldValues = cache.values();

        cache.ensureForTime(2.0, 5, TRANSITION_MATRIX, SORTED_ALPHA);

        assertTrue(cache.allocatedN() >= 5);
        assertNotSame(oldValues, cache.values());
        assertSnapshotMatches(cache, 2, beforeGrowth);

        SericolaCumulantMatrices fresh = new SericolaCumulantMatrices(DIM);
        fresh.ensureForTime(2.0, 5, TRANSITION_MATRIX, SORTED_ALPHA);
        assertSameCumulants(fresh, cache, 5);
    }

    public void testSmallerRequirementReusesExistingStorage() {
        SericolaCumulantMatrices cache = new SericolaCumulantMatrices(DIM);
        cache.ensureForTime(2.0, 5, TRANSITION_MATRIX, SORTED_ALPHA);

        double[] values = cache.values();
        int allocatedN = cache.allocatedN();
        int computedN = cache.computedN();
        double maxTime = cache.maxTime();

        cache.ensureForTime(1.0, 3, TRANSITION_MATRIX, SORTED_ALPHA);

        assertSame(values, cache.values());
        assertEquals(allocatedN, cache.allocatedN());
        assertEquals(computedN, cache.computedN());
        assertEquals(maxTime, cache.maxTime(), 0.0);
    }

    public void testInvalidationClearsComputedExtentButKeepsStorage() {
        SericolaCumulantMatrices cache = new SericolaCumulantMatrices(DIM);
        cache.ensureForTime(2.0, 5, TRANSITION_MATRIX, SORTED_ALPHA);

        double[] values = cache.values();
        int allocatedN = cache.allocatedN();

        cache.invalidateComputedExtent();

        assertSame(values, cache.values());
        assertEquals(allocatedN, cache.allocatedN());
        assertEquals(-1, cache.computedN());
        assertEquals(0.0, cache.maxTime(), 0.0);

        cache.ensureForTime(1.0, 4, TRANSITION_MATRIX, SORTED_ALPHA);
        assertSame(values, cache.values());
        assertEquals(4, cache.computedN());
    }

    public void testSecondDifferenceRowsMatchCumulantFormula() {
        SericolaCumulantMatrices cache = new SericolaCumulantMatrices(DIM);
        cache.ensureForTime(2.0, 5, TRANSITION_MATRIX, SORTED_ALPHA);
        cache.ensureSecondDifferenceCapacity();

        assertSecondDifferenceRow(cache, 1, 3);
        assertSecondDifferenceRow(cache, 2, 4);
    }

    public void testReverseUniformizationMatrixAdjointMatchesFiniteDifference() {
        final double lambda = 1.35;
        final double rewardProportion = 0.52;
        final double time = 0.75;
        final double step = 1.0e-6;

        double[] transitionAdjoint = reverseUniformizationAdjoint(
                TRANSITION_MATRIX,
                lambda,
                rewardProportion,
                time,
                UPSTREAM_DENSITY_ADJOINT);

        for (int uv = 0; uv < DIM2; uv++) {
            double[] plus = TRANSITION_MATRIX.clone();
            double[] minus = TRANSITION_MATRIX.clone();
            plus[uv] += step;
            minus[uv] -= step;

            double finiteDifference = (
                    rewardDensityObjective(plus, lambda, rewardProportion, time, UPSTREAM_DENSITY_ADJOINT) -
                            rewardDensityObjective(minus, lambda, rewardProportion, time, UPSTREAM_DENSITY_ADJOINT)) /
                    (2.0 * step);

            double tolerance = Math.max(5.0e-6, Math.abs(finiteDifference) * 5.0e-5);
            assertEquals("P entry " + uv, finiteDifference, transitionAdjoint[uv], tolerance);
        }
    }

    public void testReverseUniformizationRateAdjointMatchesFiniteDifference() {
        final double lambda = 1.35;
        final double rewardProportion = 0.52;
        final double time = 0.75;
        final double step = 1.0e-6;

        double lambdaAdjoint = reverseUniformizationRateAdjoint(
                TRANSITION_MATRIX,
                lambda,
                rewardProportion,
                time,
                UPSTREAM_DENSITY_ADJOINT);

        double finiteDifference = (
                rewardDensityObjective(TRANSITION_MATRIX, lambda + step, rewardProportion, time, UPSTREAM_DENSITY_ADJOINT) -
                        rewardDensityObjective(TRANSITION_MATRIX, lambda - step, rewardProportion, time, UPSTREAM_DENSITY_ADJOINT)) /
                (2.0 * step);

        double tolerance = Math.max(5.0e-6, Math.abs(finiteDifference) * 5.0e-5);
        assertEquals(finiteDifference, lambdaAdjoint, tolerance);
    }

    public void testReverseRewardRateAdjointMatchesFiniteDifference() {
        final double lambda = 1.35;
        final double rewardProportion = 0.52;
        final double time = 0.75;
        final double step = 1.0e-6;

        double[] rewardRateAdjoint = reverseRewardRateAdjoint(
                TRANSITION_MATRIX,
                lambda,
                rewardProportion,
                time,
                UPSTREAM_DENSITY_ADJOINT,
                SORTED_ALPHA);

        for (int i = 0; i < DIM; i++) {
            double[] plusAlpha = SORTED_ALPHA.clone();
            double[] minusAlpha = SORTED_ALPHA.clone();
            plusAlpha[i] += step;
            minusAlpha[i] -= step;

            double finiteDifference = (
                    rewardDensityObjective(TRANSITION_MATRIX, lambda, rewardProportion, time,
                            UPSTREAM_DENSITY_ADJOINT, plusAlpha) -
                            rewardDensityObjective(TRANSITION_MATRIX, lambda, rewardProportion, time,
                                    UPSTREAM_DENSITY_ADJOINT, minusAlpha)) /
                    (2.0 * step);

            double tolerance = Math.max(2.0e-5, Math.abs(finiteDifference) * 1.0e-4);
            assertEquals("reward rate " + i, finiteDifference, rewardRateAdjoint[i], tolerance);
        }
    }

    private static void assertSecondDifferenceRow(SericolaCumulantMatrices cache, int h, int n) {
        cache.prepareSecondDifferenceRow(h, n);

        double[] C = cache.values();
        double[] D2 = cache.secondDifferences();

        for (int k = 0; k <= n - 1; k++) {
            int d2Offset = cache.secondDifferenceOffset(h, n, k);
            int c0 = cache.offset(h, n + 1, k);
            int c1 = cache.offset(h, n + 1, k + 1);
            int c2 = cache.offset(h, n + 1, k + 2);

            for (int uv = 0; uv < DIM2; uv++) {
                double expected = C[c2 + uv] - 2.0 * C[c1 + uv] + C[c0 + uv];
                assertEquals("h=" + h + " n=" + n + " k=" + k + " uv=" + uv,
                        expected, D2[d2Offset + uv], TOL);
            }
        }
    }

    private static void assertSameCumulants(
            SericolaCumulantMatrices expected,
            SericolaCumulantMatrices observed,
            int maxN) {

        double[] expectedValues = snapshot(expected, maxN);
        assertSnapshotMatches(observed, maxN, expectedValues);
    }

    private static double[] snapshot(SericolaCumulantMatrices cache, int maxN) {
        double[] snapshot = new double[validBlockCount(maxN) * DIM2];
        double[] values = cache.values();
        int index = 0;

        for (int h = 1; h <= PHI; h++) {
            for (int n = 0; n <= maxN; n++) {
                for (int k = 0; k <= n; k++) {
                    int offset = cache.offset(h, n, k);
                    for (int uv = 0; uv < DIM2; uv++) {
                        snapshot[index++] = values[offset + uv];
                    }
                }
            }
        }

        return snapshot;
    }

    private static void assertSnapshotMatches(SericolaCumulantMatrices cache, int maxN, double[] expected) {
        assertEquals(validBlockCount(maxN) * DIM2, expected.length);

        double[] values = cache.values();
        int index = 0;

        for (int h = 1; h <= PHI; h++) {
            for (int n = 0; n <= maxN; n++) {
                for (int k = 0; k <= n; k++) {
                    int offset = cache.offset(h, n, k);
                    for (int uv = 0; uv < DIM2; uv++) {
                        assertEquals("h=" + h + " n=" + n + " k=" + k + " uv=" + uv,
                                expected[index++], values[offset + uv], TOL);
                    }
                }
            }
        }
    }

    private static int validBlockCount(int maxN) {
        return PHI * (maxN + 1) * (maxN + 2) / 2;
    }

    private static double[] reverseUniformizationAdjoint(
            double[] transitionMatrix,
            double lambda,
            double rewardProportion,
            double time,
            double[] densityAdjoint) {

        SericolaRewardDensityGradient gradient = computeReverseGradient(
                transitionMatrix,
                lambda,
                rewardProportion,
                time,
                densityAdjoint);
        return gradient.transitionMatrixAdjoint().clone();
    }

    private static double reverseUniformizationRateAdjoint(
            double[] transitionMatrix,
            double lambda,
            double rewardProportion,
            double time,
            double[] densityAdjoint) {

        SericolaCumulantMatrices cache = new SericolaCumulantMatrices(DIM);
        SericolaRewardDensityWorkspace workspace = new SericolaRewardDensityWorkspace(DIM, 1.0e-10);
        int h = workspace.prepareDerivative(rewardProportion, SORTED_ALPHA, INV_ALPHA_DIFF);
        int requiredN = workspace.determineNumberOfSteps(lambda, time) + 1;
        cache.ensureForTime(time, requiredN, transitionMatrix, SORTED_ALPHA);

        SericolaRewardDensityGradient gradient =
                new SericolaRewardDensityGradient(DIM, OUT_ROW_BASE_BY_SORTED, OUT_COL_BY_SORTED);
        return gradient.computeWrtUniformizationMatrixInto(
                densityAdjoint,
                h,
                cache.computedN() - 1,
                time,
                lambda,
                INV_ALPHA_DIFF[h],
                workspace.isZero(0),
                workspace.isOne(0),
                workspace.xh(0),
                false,
                transitionMatrix,
                SORTED_ALPHA,
                cache);
    }

    private static double[] reverseRewardRateAdjoint(
            double[] transitionMatrix,
            double lambda,
            double rewardProportion,
            double time,
            double[] densityAdjoint,
            double[] sortedAlpha) {

        double[] invAlphaDiff = inverseAlphaDiff(sortedAlpha);
        SericolaCumulantMatrices cache = new SericolaCumulantMatrices(DIM);
        SericolaRewardDensityWorkspace workspace = new SericolaRewardDensityWorkspace(DIM, 1.0e-10);
        int h = workspace.prepareDerivative(rewardProportion, sortedAlpha, invAlphaDiff);
        int requiredN = workspace.determineNumberOfSteps(lambda, time) + 1;
        cache.ensureForTime(time, requiredN, transitionMatrix, sortedAlpha);

        double[] rewardRateAdjoint = new double[DIM];
        SericolaRewardDensityGradient gradient =
                new SericolaRewardDensityGradient(DIM, OUT_ROW_BASE_BY_SORTED, OUT_COL_BY_SORTED);
        gradient.computeWrtRewardRatesInto(
                densityAdjoint,
                h,
                cache.computedN() - 1,
                time,
                lambda,
                invAlphaDiff[h],
                workspace.isZero(0),
                workspace.isOne(0),
                workspace.xh(0),
                false,
                transitionMatrix,
                sortedAlpha,
                cache,
                rewardRateAdjoint);
        return rewardRateAdjoint;
    }

    private static SericolaRewardDensityGradient computeReverseGradient(
            double[] transitionMatrix,
            double lambda,
            double rewardProportion,
            double time,
            double[] densityAdjoint) {

        SericolaCumulantMatrices cache = new SericolaCumulantMatrices(DIM);
        SericolaRewardDensityWorkspace workspace = new SericolaRewardDensityWorkspace(DIM, 1.0e-10);
        int h = workspace.prepareDerivative(rewardProportion, SORTED_ALPHA, INV_ALPHA_DIFF);
        int requiredN = workspace.determineNumberOfSteps(lambda, time) + 1;
        cache.ensureForTime(time, requiredN, transitionMatrix, SORTED_ALPHA);

        SericolaRewardDensityGradient gradient =
                new SericolaRewardDensityGradient(DIM, OUT_ROW_BASE_BY_SORTED, OUT_COL_BY_SORTED);
        gradient.computeWrtUniformizationMatrixInto(
                densityAdjoint,
                h,
                cache.computedN() - 1,
                time,
                lambda,
                INV_ALPHA_DIFF[h],
                workspace.isZero(0),
                workspace.isOne(0),
                workspace.xh(0),
                false,
                transitionMatrix,
                SORTED_ALPHA,
                cache);
        return gradient;
    }

    private static double rewardDensityObjective(
            double[] transitionMatrix,
            double lambda,
            double rewardProportion,
            double time,
            double[] densityAdjoint) {

        return rewardDensityObjective(transitionMatrix, lambda, rewardProportion, time,
                densityAdjoint, SORTED_ALPHA);
    }

    private static double rewardDensityObjective(
            double[] transitionMatrix,
            double lambda,
            double rewardProportion,
            double time,
            double[] densityAdjoint,
            double[] sortedAlpha) {

        final double[] invAlphaDiff = inverseAlphaDiff(sortedAlpha);

        SericolaCumulantMatrices cache = new SericolaCumulantMatrices(DIM);
        SericolaRewardDensityWorkspace workspace = new SericolaRewardDensityWorkspace(DIM, 1.0e-10);
        int requiredN = workspace.determineNumberOfSteps(lambda, time) + 1;
        cache.ensureForTime(time, requiredN, transitionMatrix, sortedAlpha);

        double[][] density = new double[][]{new double[DIM2]};
        workspace.preparePdf(
                new double[]{rewardProportion},
                new double[]{time},
                false,
                lambda,
                sortedAlpha,
                invAlphaDiff);
        new SericolaRewardDensityPdf(DIM, OUT_ROW_BASE_BY_SORTED, OUT_COL_BY_SORTED).accumulateInto(
                density,
                1,
                cache.computedN() - 1,
                lambda,
                invAlphaDiff,
                cache,
                workspace);

        return dot(densityAdjoint, density[0]);
    }

    private static double dot(double[] left, double[] right) {
        double sum = 0.0;
        for (int i = 0; i < left.length; i++) {
            sum += left[i] * right[i];
        }
        return sum;
    }

    private static double[] inverseAlphaDiff() {
        return inverseAlphaDiff(SORTED_ALPHA);
    }

    private static double[] inverseAlphaDiff(double[] sortedAlpha) {
        double[] inverse = new double[DIM];
        for (int h = 1; h < DIM; h++) {
            inverse[h] = 1.0 / (sortedAlpha[h] - sortedAlpha[h - 1]);
        }
        return inverse;
    }
}
