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
import dr.math.matrix.SymmetricMatrixAccumulator;

import java.util.Arrays;
import java.util.List;

/**
 * Standalone squared B-spline model for the skyspline coalescent.
 *
 * Rate function:  λ(t) = ε + f(t)²,    f(t) = intercept + θ'b(t)
 * Augmented parameter:  u = [intercept, θ_0, ..., θ_{p-1}]
 * Augmented basis:      φ(t) = [1, b_0(t), ..., b_{p-1}(t)]
 * Squared integral:     ∫_a^b f(t)² dt  =  u' M(a,b) u
 *
 * This class owns the B-spline setup and exposes the three primitives needed
 * by ExactTimeVaryingCoalescentLikelihood:
 *
 *   addScaledGramMatrix(start, end, weight, G)  — accumulate weight*M(a,b) into G
 *   evaluateSpline(t)                           — f(t) for log λ(t)
 *   evaluateAugmentedBasisInPlace(t, out)       — φ(t) for ∇ log λ(t)
 *
 * The aggregate weighted Gram matrix G = Σ_r A_r M_r is managed by the
 * coalescent likelihood, not here.
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
    // Buffers (allocated once)
    // -------------------------------------------------------------------------

    /** [intercept, θ_0, ..., θ_{p-1}] — filled by getCombinedParamBuffer(). */
    private final double[] combinedParamBuffer;

    /** [1, b_0(t), ..., b_{p-1}(t)] — filled by evaluateAugmentedBasisInPlace(). */
    private final double[] augmentedBasisBuffer;

    /**
     * Pre-allocated polynomial piece arrays for the assembled f(t).
     * BSpline.PPoly does NOT clone pieces, so splineFunction.pieces == mutablePieces.
     * refreshSplineFunction() fills these in-place; evaluate() sees the update immediately.
     */
    private final double[][] mutablePieces;

    /** PPoly wrapping mutablePieces — pre-built at construction, never reassigned. */
    private final BSpline.PPoly splineFunction;

    /** Snapshot of [intercept, θ...] at the last refresh miss. */
    private double[] cachedValues;

    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------

    public SquaredCachedSplines(Parameter coefficient,
                                 Parameter intercept,
                                 double[] knots,
                                 double lowerBoundary,
                                 double upperBoundary,
                                 int degree) {
        if (coefficient.getDimension() != knots.length + degree)
            throw new IllegalArgumentException(
                    "Coefficient dimension must equal knots.length + degree");

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
        for (int i = 0; i < knots.length; i++)
            ek[degree + i + 1] = knots[i];
        return ek;
    }

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    public Parameter getCoefficients()          { return coefficient; }
    public Parameter getIntercept()             { return intercept; }
    public int getCoefficientDim()              { return coefficient.getDimension(); }
    public BSplineGramMatrix getBSplineGramMatrix() { return gramMatrix; }

    // -------------------------------------------------------------------------
    // Core hot-path primitives (called by the coalescent likelihood)
    // -------------------------------------------------------------------------

    /**
     * Accumulates weight * M(start, end) into G in-place.
     * Called once per interval per tree change (when building G = Σ_r A_r M_r).
     * G must be (dim+1)×(dim+1) and zeroed before the first call.
     */
    public void addScaledGramMatrix(double start, double end, double weight, double[][] G) {
        gramMatrix.addScaledToMatrix(start, end, weight, G);
    }

    /**
     * Accumulates weight * M(start, end) into G (SymmetricMatrixAccumulator, augmented dim+1).
     * Intercept is at index 0; coefficient indices are 1..dim.
     */
    public void addScaledGramMatrix(double start, double end, double weight,
                                    SymmetricMatrixAccumulator G) {
        gramMatrix.addScaledToAccumulator(start, end, weight, G);
    }

    /**
     * Accumulates weight * M_base(start, end) into G (dim×dim, no intercept).
     */
    public void addScaledGramMatrixNoIntercept(double start, double end, double weight,
                                               SymmetricMatrixAccumulator G) {
        gramMatrix.addScaledToAccumulatorNoIntercept(start, end, weight, G);
    }

    /**
     * Fills out[0..dim-1] with [b_0(t), ..., b_{dim-1}(t)] (no intercept term).
     * out must have length >= getCoefficientDim().  No allocation.
     */
    public void evaluateBasisInPlace(double t, double[] out) {
        for (int i = 0; i < coefficient.getDimension(); i++) {
            out[i] = basis.get(i + 1).evaluate(t);
        }
    }

    /**
     * Fills out[0..dim] with [end-start, integral_start^end b_0 dt, ..., integral b_{dim-1} dt].
     * out must have length >= getCoefficientDim() + 1.
     */
    public void integrateAugmentedBasisInPlace(double start, double end, double[] out) {
        gramMatrix.integrateAugmentedBasisInPlace(start, end, out);
    }

    /**
     * Fills out[0..dim-1] with [integral_start^end b_0 dt, ..., integral b_{dim-1} dt].
     * out must have length >= getCoefficientDim().
     */
    public void integrateBasisInPlace(double start, double end, double[] out) {
        gramMatrix.integrateBasisInPlace(start, end, out);
    }

    /**
     * Returns f(t) = intercept + θ'b(t).
     * Used by the coalescent to evaluate log(ε + f²) at event times.
     * Zero-allocation on a cache hit (parameters unchanged).
     */
    public double evaluateSpline(double t) {
        refreshSplineFunction();
        return splineFunction.evaluate(t);
    }

    /**
     * Fills out[0..dim] with φ(t) = [1, b_0(t), ..., b_{dim-1}(t)] in-place.
     * Used by the coalescent gradient: ∇ log(ε + f²) = (2f/(ε+f²)) · φ(t).
     * out must have length >= getCoefficientDim() + 1.  No allocation.
     */
    public void evaluateAugmentedBasisInPlace(double t, double[] out) {
        out[0] = 1.0;
        for (int i = 0; i < coefficient.getDimension(); i++)
            out[i + 1] = basis.get(i + 1).evaluate(t);
    }

    /**
     * Returns the current u = [intercept, θ_0, ..., θ_{dim-1}] buffer.
     * Zero-allocation.  Contract: callers must not retain the reference across
     * any operation that could modify the spline parameters.
     */
    public double[] getCombinedParamBuffer() {
        combinedParamBuffer[0] = intercept.getParameterValue(0);
        for (int i = 0; i < coefficient.getDimension(); i++)
            combinedParamBuffer[i + 1] = coefficient.getParameterValue(i);
        return combinedParamBuffer;
    }

    // -------------------------------------------------------------------------
    // Utility methods (not in the HMC hot path; useful for debugging / tests)
    // -------------------------------------------------------------------------

    /** u' M(start,end) u — the squared-spline integral over [start,end]. */
    public double getIntegral(double start, double end) {
        return gramMatrix.quadraticForm(start, end, getCombinedParamBuffer());
    }

    /** d/d_intercept ∫ f² dt over [start,end]. */
    public double getGradientWrtIntercept(double start, double end) {
        if (!(end > start)) return 0.0;
        return gramMatrix.quadraticForm(start, end, getCombinedParamBuffer()) > 0 ?
                buildGradientFromMatrix(start, end)[0] : 0.0;
    }

    /** d/d_theta_i ∫ f² dt over [start,end]. */
    public double getGradientWrtCoefficient(double start, double end, int i) {
        if (i < 0 || i >= coefficient.getDimension())
            throw new IllegalArgumentException("Coefficient index out of bounds: " + i);
        if (!(end > start)) return 0.0;
        return buildGradientFromMatrix(start, end)[i + 1];
    }

    private double[] buildGradientFromMatrix(double start, double end) {
        double[][] M = gramMatrix.computeMatrix(start, end);
        double[]   u = getCombinedParamBuffer();
        double[]   g = new double[u.length];
        BSplineGramMatrix.addMatVecReordered(M, u, 2.0, g);
        // reorder back to u-layout for the caller
        double[] result = new double[u.length];
        result[0] = g[u.length - 1];           // intercept
        System.arraycopy(g, 0, result, 1, u.length - 1);
        return result;
    }

    // -------------------------------------------------------------------------
    // Private: zero-allocation spline function refresh
    // -------------------------------------------------------------------------

    /**
     * Fills mutablePieces[k][j] in-place when parameters have changed.
     * splineFunction wraps mutablePieces so evaluate() sees updates immediately.
     * One clone per miss; zero allocation on a cache hit.
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
