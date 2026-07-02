/*
 * SquaredSplineCoalescentLikelihood.java
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

package dr.evomodel.coalescent;

import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.SquaredCachedSplines;
import dr.util.Author;
import dr.util.Citation;

import java.util.Arrays;
import java.util.List;

/**
 * Coalescent likelihood for λ(t) = ε + f(t)²,  f(t) = intercept + θ'b(t).
 *
 * Implements the three primitives required by ExactTimeVaryingCoalescentLikelihood:
 *   accumulateToWeightedGramMatrix — delegates to SquaredCachedSplines.addScaledGramMatrix
 *   evaluateRateAt(t)              — ε + f(t)²
 *   addEventGradient(t, grad)      — (2f/λ) φ(t) accumulated into grad
 *
 * The aggregate weighted Gram matrix G = Σ_r A_r M_r is owned and cached by the
 * base class.  It is constant within an HMC trajectory and invalidated only when
 * the tree changes.
 *
 * @author Filippo Monti
 * @author Marc A. Suchard
 */
public class SquaredSplineCoalescentLikelihood extends ExactTimeVaryingCoalescentLikelihood {

    private final SquaredCachedSplines splines;
    private final Parameter epsilon;        // positivity floor ε
    private final Parameter rateParameter;  // compound [coefficients, intercept]

    // Pre-allocated for addEventGradient: φ(t) = [1, b_0(t), ..., b_{dim-1}(t)]
    private final double[] basisBuffer;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /**
     * @param rateParameter  The compound [coefficients, intercept] shared by the likelihood,
     *                       the HMC operator, and any prior gradient wrapper.
     *                       Must wrap exactly splines.getCoefficients() and splines.getIntercept().
     */
    public SquaredSplineCoalescentLikelihood(List<BigFastTreeIntervals> intervalsList,
                                              SquaredCachedSplines splines,
                                              Parameter epsilon,
                                              Parameter rateParameter,
                                              Parameter ploidyFactors) {
        super("SquaredSplineCoalescentLikelihood",
              intervalsList, ploidyFactors,
              splines.getCoefficientDim() + 1);  // augmentedDim = dim + 1

        if (epsilon.getDimension() != 1 || epsilon.getParameterValue(0) < 0.0)
            throw new IllegalArgumentException("epsilon must be a scalar ≥ 0");

        if (epsilon.getParameterValue(0) == 0.0)
            java.util.logging.Logger.getLogger(getClass().getName()).warning(
                    "epsilon=0: logHazard returns -Infinity when f(t)=0 at a coalescent event.");

        if (rateParameter.getDimension() != splines.getCoefficientDim() + 1)
            throw new IllegalArgumentException(
                    "rateParameter dimension must equal spline coefficients + 1 (intercept)");

        // Identity check: rateParameter must wrap the same objects as splines' coefficients/intercept.
        if (rateParameter instanceof CompoundParameter) {
            CompoundParameter cp = (CompoundParameter) rateParameter;
            boolean coefMatch = false, interceptMatch = false;
            for (int k = 0; k < cp.getParameterCount(); k++) {
                Parameter sub = cp.getParameter(k);
                if (sub == splines.getCoefficients()) coefMatch     = true;
                if (sub == splines.getIntercept())    interceptMatch = true;
            }
            if (!coefMatch || !interceptMatch)
                throw new IllegalArgumentException(
                        "rateParameter must wrap exactly splines.getCoefficients() " +
                        "and splines.getIntercept() (identity check failed).");
        }

        this.splines       = splines;
        this.epsilon       = epsilon;
        this.rateParameter = rateParameter;
        this.basisBuffer   = new double[rateParameter.getDimension()];  // dim+1

        addVariable(splines.getCoefficients());
        addVariable(splines.getIntercept());
        addVariable(epsilon);
    }

    // -----------------------------------------------------------------------
    // Template method implementations
    // -----------------------------------------------------------------------

    /**
     * Adds A_r * M(start, end) into G — one interval's contribution to the
     * weighted Gram matrix G = Σ_r A_r M_r.
     */
    @Override
    protected void accumulateToWeightedGramMatrix(double start, double end,
                                                   double weight, double[][] G) {
        splines.addScaledGramMatrix(start, end, weight, G);
    }

    /** λ(t) = ε + f(t)² */
    @Override
    protected double evaluateRateAt(double t) {
        double f = splines.evaluateSpline(t);
        return epsilon.getParameterValue(0) + f * f;
    }

    /**
     * Accumulates ∇_u log λ(t) = (2f/λ) φ(t) into grad (rateParameter layout).
     * grad layout: [d/dθ_0, ..., d/dθ_{dim-1}, d/d_intercept].
     */
    @Override
    protected void addEventGradient(double t, double[] grad) {
        double f     = splines.evaluateSpline(t);
        double rate  = epsilon.getParameterValue(0) + f * f;
        double scale = 2.0 * f / rate;
        splines.evaluateAugmentedBasisInPlace(t, basisBuffer);  // [1, b_0,...,b_{dim-1}]

        int dim = splines.getCoefficientDim();
        for (int i = 0; i < dim; i++)
            grad[i] += scale * basisBuffer[i + 1];   // d/dθ_i
        grad[dim] += scale * basisBuffer[0];           // d/d_intercept
    }

    /**
     * Accumulates ∇²_u log λ(t) = (2/λ − 4η²/λ²) φ φ' into H (u-layout).
     * H layout: [intercept, θ_0,...,θ_{dim-1}] × [intercept, θ...].
     */
    @Override
    protected void addEventHessian(double t, double[][] H) {
        double f     = splines.evaluateSpline(t);
        double rate  = epsilon.getParameterValue(0) + f * f;
        double eta   = f;
        double coeff = 2.0 / rate - 4.0 * eta * eta / (rate * rate);
        splines.evaluateAugmentedBasisInPlace(t, basisBuffer);  // reuse buffer

        int n = basisBuffer.length;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                H[i][j] += coeff * basisBuffer[i] * basisBuffer[j];
    }

    @Override
    protected double[] getRateParamBuffer() {
        return splines.getCombinedParamBuffer();
    }

    @Override
    protected double getEpsilonValue() {
        return epsilon.getParameterValue(0);
    }

    @Override
    public Parameter getRateParameter() { return rateParameter; }

    // -----------------------------------------------------------------------
    // Sign canonicalization on accept
    // -----------------------------------------------------------------------

    /**
     * After every accepted proposal, flip all of u if the largest-magnitude component
     * is negative.  L(u) = L(-u) exactly, so this is a zero-cost identifiability fix.
     * Uses setParameterValueQuietly to avoid re-dirtying the likelihood.
     */
    @Override
    protected void acceptState() {
        Parameter coeff     = splines.getCoefficients();
        Parameter intercept = splines.getIntercept();

        double anchor = intercept.getParameterValue(0);
        double maxAbs = Math.abs(anchor);
        for (int i = 0; i < coeff.getDimension(); i++) {
            double v = coeff.getParameterValue(i);
            double a = Math.abs(v);
            if (a > maxAbs) { maxAbs = a; anchor = v; }
        }

        if (anchor < 0.0) {
            intercept.setParameterValueQuietly(0, -intercept.getParameterValue(0));
            for (int i = 0; i < coeff.getDimension(); i++)
                coeff.setParameterValueQuietly(i, -coeff.getParameterValue(i));
        }
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    public SquaredCachedSplines getSplines() { return splines; }
    public Parameter getEpsilon()            { return epsilon; }

    // -----------------------------------------------------------------------
    // Citations
    // -----------------------------------------------------------------------

    @Override public String getDescription() { return "Squared-spline coalescent likelihood"; }

    @Override
    public List<Citation> getCitations() {
        return Arrays.asList(new Citation(
                new Author[]{ new Author("F", "Monti"), new Author("MA", "Suchard") },
                Citation.Status.IN_PREPARATION));
    }
}
