/*
 * SquaredCachedSplines.java
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

import dr.inference.model.Parameter;

import java.util.Arrays;
import java.util.List;

/**
 * Standalone squared B-spline model for the skyspline coalescent.
 *
 * Rate function: λ(t) = ε + f(t)²,   f(t) = intercept + θ'b(t)
 *
 * Integrated hazard (exact):
 *   ∫_a^b f(t)² dt  =  u' M(a,b) u
 * where u = [intercept, θ_0, ..., θ_{p-1}] and M(a,b) is the Gram matrix
 * cached by BSplineGramMatrix (theta-independent; invalidated only on tree changes).
 *
 * This class is intentionally independent of IntegratedTransformedSplines.
 * It owns its own B-spline setup and provides only the methods needed by
 * SquaredSplineCoalescentLikelihood:
 *   - evaluateSpline(t)               for log λ(t)
 *   - evaluateAugmentedBasisInPlace()  for ∇ log λ(t)
 *   - getIntegralForInterval()         for ∫ λ dt  (indexed cache, zero alloc)
 *   - addGradientWrtParametersInPlace() for ∇ ∫ λ dt (indexed cache, zero alloc)
 *   - clearGramCache() / initIndexedGramCache()  for tree-change invalidation
 *
 * All hot-path methods are zero-allocation on a cache hit.
 *
 * @author Filippo Monti
 * @author Marc A. Suchard
 */
public class SquaredCachedSplines {

    private final Parameter intercept;
    private final Parameter coefficient;
    private final double upperBoundary;
    private final double[] expandedKnots;
    private final int degree;
    private final List<BSpline.PPoly> basis;

    private final BSplineGramMatrix gramMatrix;

    // -------------------------------------------------------------------------
    // Hot-path buffers (allocated once at construction)
    // -------------------------------------------------------------------------

    /** [intercept, θ_0, ..., θ_{p-1}] — filled by getCombinedParamBuffer(). */
    private final double[] combinedParamBuffer;

    /** [1, b_0(t), ..., b_{p-1}(t)] — filled by evaluateAugmentedBasisInPlace(). */
    private final double[] augmentedBasisBuffer;

    /**
     * Pre-allocated piece arrays for the assembled f(t) = intercept + θ'b(t).
     * BSpline.PPoly does NOT clone pieces, so splineFunction.pieces == mutablePieces.
     * refreshSplineFunction() fills these in-place; evaluate() sees updates immediately.
     */
    private final double[][] mutablePieces;

    /** PPoly wrapping mutablePieces — pre-built once, never reassigned. */
    private final BSpline.PPoly splineFunction;

    /** Snapshot of [intercept, θ...] at the last cache miss. */
    private double[] cachedValues;

    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------

    /**
     * @param coefficient  spline coefficients θ (length = knots.length + degree)
     * @param intercept    scalar intercept parameter
     * @param knots        interior knot positions, sorted, in (lowerBoundary, upperBoundary)
     * @param lowerBoundary left end of the spline domain
     * @param upperBoundary right end; constant extension (f = f(upperBoundary)) beyond this
     * @param degree       B-spline polynomial degree (typically 3 = cubic)
     */
    public SquaredCachedSplines(Parameter coefficient,
                                 Parameter intercept,
                                 double[] knots,
                                 double lowerBoundary,
                                 double upperBoundary,
                                 int degree) {
        if (coefficient.getDimension() != knots.length + degree) {
            throw new IllegalArgumentException(
                    "Coefficient dimension must equal knots.length + degree");
        }

        this.intercept     = intercept;
        this.coefficient   = coefficient;
        this.upperBoundary = upperBoundary;
        this.degree        = degree;

        this.expandedKnots = buildExpandedKnots(knots, lowerBoundary, upperBoundary, degree);
        this.basis         = BSpline.bSplineBasis(expandedKnots, degree);
        this.gramMatrix    = new BSplineGramMatrix(
                basis, expandedKnots, degree, coefficient.getDimension(), upperBoundary);

        int dim    = coefficient.getDimension();
        int nSpans = expandedKnots.length - 1;
        this.combinedParamBuffer  = new double[1 + dim];
        this.augmentedBasisBuffer = new double[1 + dim];
        this.mutablePieces        = new double[nSpans][degree + 1];
        this.splineFunction       = new BSpline.PPoly(expandedKnots, mutablePieces);
        this.cachedValues         = null;
    }

    private static double[] buildExpandedKnots(double[] knots, double lower,
                                                double upper, int degree) {
        double[] ek = new double[knots.length + 2 * (degree + 1)];
        for (int i = 0; i <= degree; i++) {
            ek[i] = lower;
            ek[degree + knots.length + i + 1] = upper;
        }
        for (int i = 0; i < knots.length; i++) {
            ek[degree + i + 1] = knots[i];
        }
        return ek;
    }

    // -------------------------------------------------------------------------
    // Cache management
    // -------------------------------------------------------------------------

    /**
     * Activates the index-based Gram matrix cache for {@code totalIntervals} slots.
     * Call once after the tree intervals are known (interval count is fixed for a given dataset).
     * Subsequent clearGramCache() calls reset entries to null without reallocation.
     */
    public void initIndexedGramCache(int totalIntervals) {
        gramMatrix.initIndexedCache(totalIntervals);
    }

    /**
     * Invalidates all cached Gram matrices.  Call on every tree topology or node-time change.
     * No-op when only θ or ε change — M(a,b) is theta-independent.
     */
    public void clearGramCache() {
        gramMatrix.clearCache();
    }

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    public Parameter getCoefficients()          { return coefficient; }
    public Parameter getIntercept()             { return intercept; }
    public int getCoefficientDim()              { return coefficient.getDimension(); }
    public BSplineGramMatrix getBSplineGramMatrix() { return gramMatrix; }

    // -------------------------------------------------------------------------
    // General-purpose integral and gradient (hash-cache path)
    // -------------------------------------------------------------------------

    /** Returns u' M(start,end) u = ∫_start^end f(t)² dt via the hash-keyed gram matrix. */
    public double getIntegral(double start, double end) {
        return !(end > start) ? 0.0 :
                gramMatrix.quadraticForm(start, end, getCombinedParamBuffer());
    }

    /** Returns d/d_intercept ∫_start^end f² dt = 2 M(start,end) u at row 0. */
    public double getGradientWrtIntercept(double start, double end) {
        if (!(end > start)) return 0.0;
        return gramMatrix.gradient(start, end, getCombinedParamBuffer())[0];
    }

    /** Returns d/d_theta_i ∫_start^end f² dt = 2 M(start,end) u at row i+1. */
    public double getGradientWrtCoefficient(double start, double end, int i) {
        if (i < 0 || i >= coefficient.getDimension())
            throw new IllegalArgumentException("Coefficient index out of bounds: " + i);
        if (!(end > start)) return 0.0;
        return gramMatrix.gradient(start, end, getCombinedParamBuffer())[i + 1];
    }

    // -------------------------------------------------------------------------
    // Hot-path evaluation methods (called per coalescent event / interval)
    // -------------------------------------------------------------------------

    /**
     * Evaluates f(t) = intercept + θ'b(t).
     * Used by SquaredSplineCoalescentLikelihood.logHazard() to compute log(ε + f²).
     * Zero-allocation on a cache hit (parameters unchanged since last call).
     */
    public double evaluateSpline(double t) {
        refreshSplineFunction();
        return splineFunction.evaluate(t);
    }

    /**
     * Fills out[0..dim] with φ(t) = [1, b_0(t), ..., b_{p-1}(t)] in-place.
     * Used by addLogHazardGradient() to compute ∇ log(ε + f²) = (2f/(ε+f²)) · φ(t).
     * out must have length >= getCoefficientDim() + 1.  No allocation.
     */
    public void evaluateAugmentedBasisInPlace(double t, double[] out) {
        out[0] = 1.0;
        for (int i = 0; i < coefficient.getDimension(); i++) {
            out[i + 1] = basis.get(i + 1).evaluate(t);
        }
    }

    /**
     * Returns u' M(globalIdx, start, end) u — the squared-spline integral ∫_start^end f² dt.
     * Looks up M by global interval index: no IntervalKey allocation on a cache hit.
     * Requires initIndexedGramCache() to have been called.
     */
    public double getIntegralForInterval(int globalIdx, double start, double end) {
        return !(end > start) ? 0.0 :
                gramMatrix.quadraticFormByIndex(globalIdx, start, end, getCombinedParamBuffer());
    }

    /**
     * Adds coef * 2 M(globalIdx, start, end) u into grad in rateParameter layout.
     * grad layout: [d/d_θ_0, ..., d/d_θ_{p-1}, d/d_intercept].
     * Zero-allocation on a cache hit.  Requires initIndexedGramCache() to have been called.
     */
    public void addGradientWrtParametersInPlace(int globalIdx, double start, double end,
                                                 double coef, double[] grad) {
        if (!(end > start)) return;
        gramMatrix.addScaledGradientReorderedByIndex(
                globalIdx, start, end, getCombinedParamBuffer(), coef, grad,
                coefficient.getDimension());
    }

    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    /**
     * Fills combinedParamBuffer with [intercept, θ_0, ..., θ_{p-1}] and returns it.
     * Zero-allocation.  Contract: callers must not retain the reference across parameter changes.
     */
    private double[] getCombinedParamBuffer() {
        combinedParamBuffer[0] = intercept.getParameterValue(0);
        for (int i = 0; i < coefficient.getDimension(); i++) {
            combinedParamBuffer[i + 1] = coefficient.getParameterValue(i);
        }
        return combinedParamBuffer;
    }

    /**
     * Fills mutablePieces[k][j] in-place when parameters have changed.
     *
     * f_k(t) = intercept + Σ θ_i b_{i,k}(t)  on knot span k.
     *
     * splineFunction wraps mutablePieces (PPoly doesn't clone pieces), so evaluate()
     * sees the updated coefficients immediately — no PPoly rebuild, no allocation on hits.
     * One clone per miss to store a stable cachedValues snapshot.
     */
    private void refreshSplineFunction() {
        double[] current = getCombinedParamBuffer();
        if (Arrays.equals(current, cachedValues)) return;

        int nSpans = expandedKnots.length - 1;
        for (int k = 0; k < nSpans; k++) {
            double[] fk = mutablePieces[k];
            fk[0] = current[0];
            for (int j = 1; j <= degree; j++) fk[j] = 0.0;
            for (int i = 0; i < coefficient.getDimension(); i++) {
                double theta = current[i + 1];
                double[] bik = basis.get(i + 1).pieces[k];
                for (int j = 0; j < bik.length; j++) fk[j] += theta * bik[j];
            }
        }
        cachedValues = current.clone();
    }
}
