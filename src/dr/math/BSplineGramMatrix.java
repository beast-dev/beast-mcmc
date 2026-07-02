/*
 * BSplineGramMatrix.java
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

package dr.math;

import java.util.List;

/**
 * Computes the integrated basis-product Gram matrix M(a,b) for a B-spline basis.
 *
 * For augmented parameter vector u = [intercept, theta_0, ..., theta_p] and
 * augmented basis phi(t) = [1, b_0(t), ..., b_p(t)]:
 *
 *   integral_a^b (u' phi(t))^2 dt  =  u' M(a,b) u
 *
 * M(a,b) is theta-independent and computed exactly via polynomial integration.
 *
 * This class is a pure computation engine — it does not cache.  Callers that
 * need caching (e.g. the weighted Gram matrix G in the coalescent) are responsible
 * for their own storage strategy.
 *
 * @author Filippo Monti
 * @author Marc A. Suchard
 */
public class BSplineGramMatrix {

    private final List<BSpline.PPoly> basis;
    private final double[] expandedKnots;
    private final int degree;
    private final int dim;               // number of spline coefficients (excludes intercept)
    private final double upperBoundary;

    // phi(upperBoundary) = [1, b_0(ub), ..., b_{dim-1}(ub)] — for the tail correction.
    private final double[] boundaryBasis;

    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------

    public BSplineGramMatrix(List<BSpline.PPoly> basis,
                              double[] expandedKnots,
                              int degree,
                              int dim,
                              double upperBoundary) {
        this.basis         = basis;
        this.expandedKnots = expandedKnots;
        this.degree        = degree;
        this.dim           = dim;
        this.upperBoundary = upperBoundary;

        this.boundaryBasis = new double[dim + 1];
        this.boundaryBasis[0] = 1.0;
        for (int i = 0; i < dim; i++) {
            this.boundaryBasis[i + 1] = basis.get(i + 1).evaluate(upperBoundary);
        }
    }

    // -------------------------------------------------------------------------
    // Public API
    // -------------------------------------------------------------------------

    /**
     * Accumulates weight * M(a,b) into the existing matrix G in-place.
     *
     * G must be (dim+1) x (dim+1) and fully symmetric (both triangles filled).
     * On entry G is read and updated; the caller owns G and is responsible for
     * zeroing it before the first accumulation.
     */
    public void addScaledToMatrix(double a, double b, double weight, double[][] G) {
        if (!(b > a) || weight == 0.0) return;

        // Compute M(a,b) — allocates once per call.
        double[][] M = computeMatrix(a, b);
        int n = G.length;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                G[i][j] += weight * M[i][j];
    }

    /**
     * Returns u' M(a,b) u.  Allocates M(a,b) each call; intended for utility/test use.
     * For the hot path use addScaledToMatrix to build the aggregate G instead.
     */
    public double quadraticForm(double a, double b, double[] u) {
        if (!(b > a)) return 0.0;
        return symmetricQuadraticForm(computeMatrix(a, b), u);
    }

    // -------------------------------------------------------------------------
    // Private: matrix computation
    // -------------------------------------------------------------------------

    /** Computes and returns the full symmetric (dim+1)×(dim+1) Gram matrix for [a,b]. */
    public double[][] computeMatrix(double a, double b) {
        int n = dim + 1;
        double[][] M = new double[n][n];

        if (b <= upperBoundary) {
            addOnSplines(a, b, M);
        } else if (a < upperBoundary) {
            addOnSplines(a, upperBoundary, M);
            addTail(b - upperBoundary, M);
        } else {
            addTail(b - a, M);
        }

        // Mirror upper triangle to lower.
        for (int i = 1; i < n; i++)
            for (int j = 0; j < i; j++)
                M[i][j] = M[j][i];

        return M;
    }

    /**
     * Accumulates the span-by-span contributions to M over [a,b] ⊆ [lowerBoundary, upperBoundary].
     * Only the upper triangle (i <= j) is filled.
     */
    private void addOnSplines(double a, double b, double[][] M) {
        for (int k = 0; k < expandedKnots.length - 1; k++) {
            double left  = expandedKnots[k];
            double right = expandedKnots[k + 1];
            if (right <= a || left >= b) continue;
            double L = Math.max(a, left);
            double R = Math.min(b, right);
            if (L >= R) continue;

            // M[0][0]: integral of 1*1 = (R-L)
            M[0][0] += R - L;

            // M[0][j+1]: integral of 1 * b_j
            for (int j = 0; j < dim; j++)
                M[0][j + 1] += BSpline.polyIntegral(basis.get(j + 1).pieces[k], L, R);

            // M[i+1][j+1]: integral of b_i * b_j (banded, upper triangle only)
            for (int j = 0; j < dim; j++) {
                double[] pj = basis.get(j + 1).pieces[k];
                int iStart = Math.max(0, j - degree);
                for (int i = iStart; i <= j; i++) {
                    double[] pi = basis.get(i + 1).pieces[k];
                    M[i + 1][j + 1] += BSpline.polyIntegral(BSpline.polyMultiply(pi, pj), L, R);
                }
            }
        }
    }

    /** Rank-1 tail contribution for constant extension of length {@code length} beyond upperBoundary. */
    private void addTail(double length, double[][] M) {
        if (length <= 0.0) return;
        int n = dim + 1;
        for (int j = 0; j < n; j++) {
            double bj = boundaryBasis[j];
            for (int i = 0; i <= j; i++)
                M[i][j] += length * boundaryBasis[i] * bj;
        }
    }

    // -------------------------------------------------------------------------
    // Package-private linear algebra helpers (used by callers that have G)
    // -------------------------------------------------------------------------

    /** u' M u for a full symmetric M. */
    public static double symmetricQuadraticForm(double[][] M, double[] u) {
        int n = u.length;
        double total = 0.0;
        for (int i = 0; i < n; i++) {
            double ui = u[i];
            for (int j = 0; j < n; j++)
                total += ui * M[i][j] * u[j];
        }
        return total;
    }

    /**
     * Adds scale * M * u into out (M is full symmetric).
     * Reordering: M is in u-layout [intercept, theta_0, ..., theta_{dim-1}];
     * out is in rateParameter layout [theta_0, ..., theta_{dim-1}, intercept].
     */
    public static void addMatVecReordered(double[][] M, double[] u, double scale, double[] out) {
        int n = u.length;
        int coeffDim = n - 1;
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) sum += M[i][j] * u[j];
            double contrib = scale * sum;
            if (i == 0) out[coeffDim] += contrib;   // intercept → last position
            else        out[i - 1]    += contrib;    // theta_i → position i-1
        }
    }
}
