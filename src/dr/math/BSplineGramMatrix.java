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

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Precomputes and caches the integrated basis-product Gram matrix M(a,b) for a B-spline basis.
 *
 * For augmented parameter vector u = [intercept, theta_1, ..., theta_p] and
 * augmented basis phi(t) = [1, b_1(t), ..., b_p(t)], the squared-spline integral is
 *
 *   integral_a^b (u' phi(t))^2 dt = u' M(a,b) u
 *
 * and its gradient with respect to u is 2 M(a,b) u.
 *
 * M(a,b) depends only on the interval endpoints and the basis — not on u.  It is computed
 * exactly using polynomial products and power-rule antiderivatives, then cached keyed on (a,b)
 * so that multiple evaluations at different u (e.g. HMC leapfrog steps) pay no recomputation cost.
 *
 * Call clearCache() whenever the interval endpoints change (i.e. when tree node times change).
 *
 * @author Filippo Monti
 */
public class BSplineGramMatrix {

    // B-spline basis: basis.get(i+1) is the i-th coefficient's basis function (0-indexed), i=0..dim-1.
    // basis.get(0) is unused (convention from BSpline.bSplineBasis / IntegratedTransformedSplines).
    private final List<BSpline.PPoly> basis;
    private final double[] expandedKnots;
    private final int degree;
    private final int dim;                   // number of spline coefficients (excludes intercept)
    private final double upperBoundary;      // constant-extrapolation boundary

    // Augmented basis vector evaluated at the upper boundary: phi(upperBoundary) = [1, b_0(ub), ..., b_{p-1}(ub)].
    // Used for the rank-1 tail correction when an interval extends past upperBoundary.
    private final double[] boundaryBasis;

    // Interval-keyed hash cache: (a, b) -> (dim+1) x (dim+1) Gram matrix.
    // Used as a fallback when the indexed cache is not active (e.g. branch-rate models).
    private Map<IntervalKey, double[][]> cache;

    private static final double CACHE_PRECISION = 1e-5;

    // Index-based cache: indexedCache[globalIntervalIdx] -> Gram matrix.
    // Eliminates IntervalKey allocation on cache hits.  Activated by calling initIndexedCache().
    // Cleared to nulls (not reallocated) by clearCache() on tree change.
    private double[][][] indexedCache;

    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------

    /**
     * @param basis         result of BSpline.bSplineBasis(expandedKnots, degree); basis.get(i+1) is the i-th basis function
     * @param expandedKnots the full padded knot vector (length = dim + 2*(degree+1))
     * @param degree        B-spline degree
     * @param dim           number of spline coefficients (= expandedKnots.length - 1 - degree - 1)
     * @param upperBoundary right end of the spline domain; intervals beyond this use constant extrapolation
     */
    public BSplineGramMatrix(List<BSpline.PPoly> basis,
                              double[] expandedKnots,
                              int degree,
                              int dim,
                              double upperBoundary) {
        this.basis          = basis;
        this.expandedKnots  = expandedKnots;
        this.degree         = degree;
        this.dim            = dim;
        this.upperBoundary  = upperBoundary;
        this.cache          = new HashMap<>();

        // phi(upperBoundary): [1, b_0(ub), b_1(ub), ..., b_{dim-1}(ub)]
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
     * Returns the (dim+1) x (dim+1) augmented Gram matrix for interval [a, b], fetching from
     * cache or computing on demand.  The matrix is symmetric; only the upper triangle (i <= j)
     * is filled — use the provided quadraticForm and gradient helpers rather than reading it raw.
     */
    public double[][] getMatrix(double a, double b) {
        IntervalKey key = new IntervalKey(roundKey(a), roundKey(b));
        double[][] M = cache.get(key);
        if (M != null) return M;
        M = computeMatrix(a, b);
        cache.put(key, M);
        return M;
    }

    /**
     * Computes u' M(a,b) u — the squared-spline integral over [a, b].
     * u = [intercept, theta_0, ..., theta_{dim-1}], size dim+1.
     * Returns 0 if b <= a.
     */
    public double quadraticForm(double a, double b, double[] u) {
        if (!(b > a)) return 0.0;
        return symmetricQuadraticForm(getMatrix(a, b), u);
    }

    /**
     * Returns 2 M(a,b) u — the gradient of the squared-spline integral w.r.t. u.
     * Result order: [d/d_intercept, d/d_theta_0, ..., d/d_theta_{dim-1}], size dim+1.
     * Returns zero vector if b <= a.
     */
    public double[] gradient(double a, double b, double[] u) {
        double[] grad = new double[dim + 1];
        if (!(b > a)) return grad;
        addScaledMatVec(getMatrix(a, b), u, 2.0, grad);
        return grad;
    }

    /**
     * Adds 2 M(a,b) u into an existing gradient array (in-place, avoids allocation).
     * No-op if b <= a.
     */
    public void addGradient(double a, double b, double[] u, double[] gradOut) {
        if (!(b > a)) return;
        addScaledMatVec(getMatrix(a, b), u, 2.0, gradOut);
    }

    /**
     * Adds 2 * coef * M(a,b) * u into gradOut with intercept-first → intercept-last reordering.
     *
     * M*u layout (internal): [d/d_intercept, d/d_theta_0, ..., d/d_theta_{coeffDim-1}]
     * gradOut layout expected: [d/d_theta_0, ..., d/d_theta_{coeffDim-1}, d/d_intercept]
     *
     * This is the zero-allocation hot path for the coalescent gradient accumulation.
     * coef is typically negative (it equals -C(n,2) * ploidy, already applied by the caller).
     * No-op if b <= a.
     */
    public void addScaledGradientReordered(double a, double b, double[] u,
                                            double coef, double[] gradOut, int coeffDim) {
        if (!(b > a)) return;
        double[][] M = getMatrix(a, b);
        int n = u.length;  // dim+1
        double scale = 2.0 * coef;
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += M[i][j] * u[j];
            }
            double contrib = scale * sum;
            if (i == 0) {
                gradOut[coeffDim] += contrib;   // intercept → last position
            } else {
                gradOut[i - 1]    += contrib;   // theta_i → position i-1
            }
        }
    }

    /**
     * Activates the index-based cache for {@code totalIntervals} interval slots.
     * Called once from the coalescent likelihood constructor (interval count = 2n-2, fixed).
     * Subsequent clearCache() calls reset entries to null without reallocating the array.
     */
    public void initIndexedCache(int totalIntervals) {
        indexedCache = new double[totalIntervals][][];
    }

    /**
     * Clears the matrix cache.  Must be called whenever the interval endpoints may change,
     * i.e. on every tree topology or node-time update.  Need NOT be called when only u changes.
     */
    public void clearCache() {
        cache.clear();
        if (indexedCache != null) Arrays.fill(indexedCache, null);
    }

    // -------------------------------------------------------------------------
    // Index-based API (zero IntervalKey allocation; for coalescent use)
    // -------------------------------------------------------------------------

    /**
     * Returns M(a,b) via the index-based cache — no IntervalKey allocation.
     * Falls back to the hash cache if indexed cache is not active or idx is out of range.
     */
    private double[][] getMatrixByIndex(int idx, double a, double b) {
        if (indexedCache != null && idx >= 0 && idx < indexedCache.length) {
            double[][] M = indexedCache[idx];
            if (M == null) {
                M = computeMatrix(a, b);
                indexedCache[idx] = M;
            }
            return M;
        }
        return getMatrix(a, b);  // hash-cache fallback
    }

    /** Computes u' M(a,b) u via the indexed cache — no object allocation on a cache hit. */
    public double quadraticFormByIndex(int idx, double a, double b, double[] u) {
        if (!(b > a)) return 0.0;
        return symmetricQuadraticForm(getMatrixByIndex(idx, a, b), u);
    }

    /**
     * Adds 2 * coef * M(a,b) * u into gradOut (intercept-reordered) via the indexed cache.
     * Identical to addScaledGradientReordered but looks up M by index — zero object allocation.
     */
    public void addScaledGradientReorderedByIndex(int idx, double a, double b, double[] u,
                                                   double coef, double[] gradOut, int coeffDim) {
        if (!(b > a)) return;
        double[][] M = getMatrixByIndex(idx, a, b);
        int n = u.length;
        double scale = 2.0 * coef;
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += M[i][j] * u[j];
            }
            double contrib = scale * sum;
            if (i == 0) {
                gradOut[coeffDim] += contrib;
            } else {
                gradOut[i - 1]    += contrib;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Private: matrix computation
    // -------------------------------------------------------------------------

    private double[][] computeMatrix(double a, double b) {
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

        // Mirror upper triangle to lower triangle to produce a full symmetric matrix.
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                M[i][j] = M[j][i];
            }
        }

        return M;
    }

    /**
     * Accumulates M contributions from the interval [a, b] within the spline domain.
     * Only the upper triangle (i <= j) is written.
     *
     * Block structure of the (dim+1) x (dim+1) matrix:
     *   M[0][0]     = integral_a^b 1 * 1 dt  = b - a
     *   M[0][j+1]   = integral_a^b 1 * b_j dt   (constant x spline)
     *   M[i+1][j+1] = integral_a^b b_i * b_j dt  (spline x spline, banded: |i-j| <= degree)
     */
    private void addOnSplines(double a, double b, double[][] M) {
        for (int k = 0; k < expandedKnots.length - 1; k++) {
            double left  = expandedKnots[k];
            double right = expandedKnots[k + 1];
            if (right <= a || left >= b) continue;
            double L = Math.max(a, left);
            double R = Math.min(b, right);
            if (L >= R) continue;

            // M[0][0]: integral of constant 1 over [L, R]
            M[0][0] += R - L;

            // M[0][j+1]: integral of b_j over [L, R]
            for (int j = 0; j < dim; j++) {
                M[0][j + 1] += BSpline.polyIntegral(basis.get(j + 1).pieces[k], L, R);
            }

            // M[i+1][j+1]: integral of b_i * b_j over [L, R].
            // B-splines are banded: M[i+1][j+1] = 0 whenever |i - j| > degree,
            // so only compute the upper triangle within the band.
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

    /**
     * Accumulates the rank-1 tail contribution for a constant extension of length {@code length}
     * beyond upperBoundary:  length * phi(ub) phi(ub)'.
     * Only the upper triangle is written.
     */
    private void addTail(double length, double[][] M) {
        if (length <= 0.0) return;
        int n = dim + 1;
        for (int j = 0; j < n; j++) {
            double bj = boundaryBasis[j];
            for (int i = 0; i <= j; i++) {
                M[i][j] += length * boundaryBasis[i] * bj;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Private: linear algebra on the symmetric matrix
    // -------------------------------------------------------------------------

    /**
     * Computes u' M u using the full symmetric matrix M.
     */
    private static double symmetricQuadraticForm(double[][] M, double[] u) {
        int n = u.length;
        double total = 0.0;
        for (int i = 0; i < n; i++) {
            double ui = u[i];
            for (int j = 0; j < n; j++) {
                total += ui * M[i][j] * u[j];
            }
        }
        return total;
    }

    /**
     * Adds scale * M * u into out.  M is full symmetric.
     */
    private static void addScaledMatVec(double[][] M, double[] u, double scale, double[] out) {
        int n = u.length;
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += M[i][j] * u[j];
            }
            out[i] += scale * sum;
        }
    }

    // -------------------------------------------------------------------------
    // Cache key
    // -------------------------------------------------------------------------

    private static double roundKey(double x) {
        if (Double.isNaN(x) || Double.isInfinite(x)) return x;
        return Math.round(x / CACHE_PRECISION) * CACHE_PRECISION;
    }

    private static final class IntervalKey {
        private final double lower;
        private final double upper;

        IntervalKey(double lower, double upper) {
            this.lower = lower;
            this.upper = upper;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof IntervalKey)) return false;
            IntervalKey other = (IntervalKey) o;
            return Double.compare(lower, other.lower) == 0 &&
                    Double.compare(upper, other.upper) == 0;
        }

        @Override
        public int hashCode() {
            long b1 = Double.doubleToLongBits(lower);
            long b2 = Double.doubleToLongBits(upper);
            return (int)((b1 ^ (b1 >>> 32)) ^ ((b2 << 1) ^ (b2 >>> 31)));
        }
    }
}
