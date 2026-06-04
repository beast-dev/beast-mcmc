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
}
