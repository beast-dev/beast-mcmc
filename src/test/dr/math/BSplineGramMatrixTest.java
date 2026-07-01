/*
 * BSplineGramMatrixTest.java
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

package test.dr.math;

import dr.inference.model.Parameter;
import dr.math.BSplineGramMatrix;
import dr.math.SquaredCachedSplines;
import junit.framework.TestCase;

/**
 * Verifies BSplineGramMatrix against two independent references:
 *
 *   (1) getIntegral consistency: gramMatrix.quadraticForm(a,b,u) must equal
 *       splines.getIntegral(a,b) parametrized at the same u.
 *
 *   (2) Finite-difference gradient: gramMatrix.gradient(a,b,u) must match the
 *       numerical derivative of the quadratic form w.r.t. u.
 */
public class BSplineGramMatrixTest extends TestCase {

    private static final double[] INTERIOR_KNOTS = {0.25, 0.5, 0.75};
    private static final double   LOWER          = 0.0;
    private static final double   UPPER          = 1.0;
    private static final int      DEGREE         = 3;

    // dim = INTERIOR_KNOTS.length + DEGREE = 6 spline coefficients
    private static final double[] THETA     = {0.13, -0.27, 0.41, -0.11, 0.08, 0.29};
    private static final double   INTERCEPT = 0.21;

    // u = [intercept, theta_0, ..., theta_{dim-1}]
    private static final double[] U = {INTERCEPT, 0.13, -0.27, 0.41, -0.11, 0.08, 0.29};

    private static final double TOLERANCE = 5e-7;

    // -----------------------------------------------------------------------
    // Quadratic form consistency with getIntegral
    // -----------------------------------------------------------------------

    public void testQuadraticFormMatchesGetIntegralInsideBoundary() {
        checkQuadraticFormMatchesGetIntegral(0.1, 0.9);
    }

    public void testQuadraticFormMatchesGetIntegralCrossingBoundary() {
        checkQuadraticFormMatchesGetIntegral(0.35, 1.35);
    }

    public void testQuadraticFormMatchesGetIntegralFullyBeyondBoundary() {
        checkQuadraticFormMatchesGetIntegral(1.1, 1.8);
    }

    // -----------------------------------------------------------------------
    // Gradient correctness via finite differences
    // -----------------------------------------------------------------------

    public void testGradientMatchesFiniteDifferenceInsideBoundary() {
        checkGradientMatchesFiniteDifference(0.1, 0.9);
    }

    public void testGradientMatchesFiniteDifferenceCrossingBoundary() {
        checkGradientMatchesFiniteDifference(0.35, 1.35);
    }

    public void testGradientMatchesFiniteDifferenceFullyBeyondBoundary() {
        checkGradientMatchesFiniteDifference(1.1, 1.8);
    }

    // -----------------------------------------------------------------------
    // Cache behaviour
    // -----------------------------------------------------------------------

    public void testGradientIsZeroForEmptyInterval() {
        BSplineGramMatrix gram = buildGramMatrix();
        double[] grad = gram.gradient(0.5, 0.3, U);   // end < start
        for (double v : grad) {
            assertEquals(0.0, v, 0.0);
        }
    }

    public void testMatrixIsCachedAcrossCalls() {
        BSplineGramMatrix gram = buildGramMatrix();
        double[][] m1 = gram.getMatrix(0.1, 0.9);
        double[][] m2 = gram.getMatrix(0.1, 0.9);
        assertSame("Same interval should return cached object", m1, m2);
    }

    public void testClearCacheReturnsFreshMatrixWithSameValues() {
        BSplineGramMatrix gram = buildGramMatrix();
        double[][] m1 = gram.getMatrix(0.1, 0.9);
        gram.clearCache();
        double[][] m2 = gram.getMatrix(0.1, 0.9);
        assertNotSame("After clearCache a new object should be allocated", m1, m2);
        for (int i = 0; i < m1.length; i++) {
            for (int j = 0; j < m1[i].length; j++) {
                assertEquals(m1[i][j], m2[i][j], 0.0);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    private static SquaredCachedSplines buildSplines() {
        return new SquaredCachedSplines(
                new Parameter.Default(THETA),
                new Parameter.Default(INTERCEPT),
                INTERIOR_KNOTS, LOWER, UPPER, DEGREE);
    }

    private static BSplineGramMatrix buildGramMatrix() {
        return buildSplines().getBSplineGramMatrix();
    }

    private void checkQuadraticFormMatchesGetIntegral(double start, double end) {
        SquaredCachedSplines splines = buildSplines();
        BSplineGramMatrix gram = splines.getBSplineGramMatrix();

        double expected = splines.getIntegral(start, end);
        double actual   = gram.quadraticForm(start, end, U);

        assertEquals("quadraticForm vs getIntegral [" + start + "," + end + "]",
                expected, actual, TOLERANCE);
    }

    private void checkGradientMatchesFiniteDifference(double start, double end) {
        BSplineGramMatrix gram = buildGramMatrix();

        double[] u        = U.clone();
        double[] analytic = gram.gradient(start, end, u);
        double[] numeric  = finiteDifferenceGradient(gram, start, end, u);

        assertEquals("Gradient size", numeric.length, analytic.length);
        for (int i = 0; i < analytic.length; i++) {
            assertEquals("Gradient[" + i + "] [" + start + "," + end + "]",
                    numeric[i], analytic[i], TOLERANCE);
        }
    }

    private static double[] finiteDifferenceGradient(BSplineGramMatrix gram,
                                                      double start, double end,
                                                      double[] u) {
        final double h = 1e-6;
        double[] grad = new double[u.length];
        for (int i = 0; i < u.length; i++) {
            double orig = u[i];
            u[i] = orig + h;
            double plus = gram.quadraticForm(start, end, u);
            u[i] = orig - h;
            double minus = gram.quadraticForm(start, end, u);
            u[i] = orig;
            grad[i] = (plus - minus) / (2.0 * h);
        }
        return grad;
    }
}
