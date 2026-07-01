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
import dr.inference.model.Parameter;
import dr.math.SquaredCachedSplines;
import dr.util.Author;
import dr.util.Citation;

import java.util.Arrays;
import java.util.List;

/**
 * Coalescent likelihood for a squared-spline coalescent rate function:
 *
 *   λ(t) = ε + f(t)²,   f(t) = intercept + θ'b(t)
 *
 * The integrated hazard over [start, end] is the exact quadratic form
 *
 *   ∫λ dt = ε(end−start) + u' M(start,end) u
 *
 * where u = [intercept, θ_1,…,θ_p] and M(start,end) is the Gram matrix
 * cached in BSplineGramMatrix.  No numerical quadrature is needed.
 *
 * Cache invalidation:
 *   tree change  → clearGramCache() + likelihood dirty
 *   theta/ε change → likelihood dirty only (M(a,b) is theta-independent)
 *
 * @author Filippo Monti
 * @author Marc A. Suchard
 */
public class SquaredSplineCoalescentLikelihood extends ExactTimeVaryingCoalescentLikelihood {

    private final SquaredCachedSplines splines;
    private final Parameter epsilon;       // positivity floor ε ≥ 0, sampled by MH
    private final Parameter rateParameter; // compound [coefficients, intercept]

    // Pre-allocated for addLogHazardGradient: holds [1, b_0(t), ..., b_{dim-1}(t)].
    private final double[] basisBuffer;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /**
     * @param rateParameter  The compound parameter [coefficients, intercept] that HMC and
     *                       gradient providers operate on.  Must be supplied by the caller so
     *                       the same Java object is shared by the likelihood, the HMC operator,
     *                       and any prior gradient wrapper — preventing the anonymous-internal-
     *                       parameter identity mismatch that breaks JointGradient.
     */
    public SquaredSplineCoalescentLikelihood(List<BigFastTreeIntervals> intervalsList,
                                              SquaredCachedSplines splines,
                                              Parameter epsilon,
                                              Parameter rateParameter,
                                              Parameter ploidyFactors) {
        super("SquaredSplineCoalescentLikelihood", intervalsList, ploidyFactors);

        if (epsilon.getDimension() != 1 || epsilon.getParameterValue(0) < 0.0) {
            throw new IllegalArgumentException("epsilon must be a scalar ≥ 0");
        }
        if (epsilon.getParameterValue(0) == 0.0) {
            // epsilon=0 is mathematically allowed but dangerous: if f(t)=0 at a
            // coalescent event, logHazard returns log(0)=-Infinity, the likelihood
            // collapses to -Infinity, and HMC stalls.  Use a small positive floor
            // (e.g. 1e-10) unless the caller explicitly intends a degenerate model.
            java.util.logging.Logger.getLogger(getClass().getName()).warning(
                    "epsilon=0: logHazard returns -Infinity when f(t)=0 at a coalescent " +
                    "event. Consider epsilon > 0 to avoid likelihood singularities.");
        }
        if (rateParameter.getDimension() != splines.getCoefficientDim() + 1) {
            throw new IllegalArgumentException(
                    "rateParameter dimension must equal spline coefficients + 1 (intercept); got "
                    + rateParameter.getDimension() + " expected " + (splines.getCoefficientDim() + 1));
        }
        // Validate identity, not just dimension: rateParameter must wrap exactly the same
        // Parameter objects as splines.getCoefficients() and splines.getIntercept().
        // A CompoundParameter with different children of the same size would pass the
        // dimension check but silently decouple HMC updates from likelihood evaluation.
        if (rateParameter instanceof dr.inference.model.CompoundParameter) {
            dr.inference.model.CompoundParameter cp = (dr.inference.model.CompoundParameter) rateParameter;
            boolean coefMatch     = false;
            boolean interceptMatch = false;
            for (int k = 0; k < cp.getParameterCount(); k++) {
                dr.inference.model.Parameter sub = cp.getParameter(k);
                if (sub == splines.getCoefficients()) coefMatch = true;
                if (sub == splines.getIntercept())    interceptMatch = true;
            }
            if (!coefMatch || !interceptMatch) {
                throw new IllegalArgumentException(
                        "rateParameter must wrap exactly splines.getCoefficients() and " +
                        "splines.getIntercept() as sub-parameters (identity check failed).");
            }
        }

        this.splines       = splines;
        this.epsilon       = epsilon;
        this.rateParameter = rateParameter;
        this.basisBuffer   = new double[rateParameter.getDimension()]; // size = coeffDim+1

        // Activate index-based gram matrix cache.  The total number of intervals across all
        // trees is fixed (= 2*n_taxa - 2 per tree), so we allocate the array once here and
        // reset its entries to null on each tree change — no reallocation during MCMC.
        int totalIntervals = 0;
        for (BigFastTreeIntervals iv : intervalsList) {
            totalIntervals += iv.getIntervalCount();
        }
        splines.initIndexedGramCache(totalIntervals);

        addVariable(splines.getCoefficients());
        addVariable(splines.getIntercept());
        addVariable(epsilon);
    }

    // -----------------------------------------------------------------------
    // Template method implementations
    // -----------------------------------------------------------------------

    /**
     * ∫_start^end λ(t) dt = ε(end−start) + u'M(start,end)u
     * Delegates to the indexed gram matrix cache — no IntervalKey allocation on a cache hit.
     */
    @Override
    protected double integratedHazard(int globalIdx, double start, double end) {
        return epsilon.getParameterValue(0) * (end - start)
                + splines.getIntegralForInterval(globalIdx, start, end);
    }

    /**
     * log λ(t) = log(ε + f(t)²)
     */
    @Override
    protected double logHazard(double t) {
        double f    = splines.evaluateSpline(t);
        double rate = epsilon.getParameterValue(0) + f * f;
        return Math.log(rate);
    }

    /**
     * Adds coef * ∇_u ∫λ dt into grad, zero-allocation.
     * Uses the indexed gram matrix cache — no IntervalKey allocation on a cache hit.
     */
    @Override
    protected void addIntegratedHazardGradient(int globalIdx, double start, double end,
                                                double coef, double[] grad) {
        splines.addGradientWrtParametersInPlace(globalIdx, start, end, coef, grad);
    }

    /**
     * Adds ∇_u log(ε + f(t)²) into grad, zero-allocation.
     *
     * d/d_u log(ε + f²) = (2f / (ε + f²)) · phi(t)
     * phi(t) = [1, b_0(t), …, b_{p-1}(t)] filled into the pre-allocated basisBuffer.
     *
     * evaluateSpline() for the squared transform now bypasses the PPoly rebuild and
     * evaluates f(t) directly from the basis functions — no getCombinedParameter() alloc.
     */
    @Override
    protected void addLogHazardGradient(double t, double[] grad) {
        double f     = splines.evaluateSpline(t);
        double rate  = epsilon.getParameterValue(0) + f * f;
        double scale = 2.0 * f / rate;
        splines.evaluateAugmentedBasisInPlace(t, basisBuffer);  // fills basisBuffer in-place

        int dim = splines.getCoefficientDim();
        for (int i = 0; i < dim; i++) {
            grad[i] += scale * basisBuffer[i + 1];   // coefficient gradient
        }
        grad[dim] += scale * basisBuffer[0];          // intercept gradient (= scale)
    }

    /**
     * On tree change: clear gram matrix cache.
     * M(a,b) depends only on interval endpoints, not on theta — so this is only
     * needed when node times change, NOT when theta/ε change.
     */
    @Override
    protected void onTreeChanged() {
        splines.clearGramCache();
    }

    /**
     * After every accepted proposal, enforce the canonical sign convention:
     * find the component of u = [coefficients, intercept] with the largest absolute
     * value (the anchor); if it is negative, flip all signs.
     *
     * This is the correct place for sign canonicalization — not an MCMC operator.
     * An operator runs at an arbitrary frequency, leaving the chain in the non-canonical
     * form between invocations.  acceptState() is called after EVERY accepted step,
     * so the logged state is always canonical regardless of which operator proposed.
     *
     * The flip is a no-op on the likelihood (L(u) = L(-u) exactly) and on any
     * symmetric prior (p(u) = p(-u)), so it does not bias the posterior.
     *
     * Rule: argmax_i |u_i| is used as the anchor rather than u[0], matching the
     * Julia canonicalize_theta_sign! implementation: a near-zero entry's sign is
     * noise-dominated and a poor anchor; the largest entry's sign is well-determined.
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
            // Use setParameterValueQuietly to avoid firing handleVariableChangedEvent,
            // which would dirty likelihoodKnown immediately after acceptance and force
            // a redundant recomputation on the next getLogLikelihood() call.
            // The sign flip does not change the likelihood (L(u) = L(-u)) or the gram
            // matrix cache, so no invalidation is needed.
            intercept.setParameterValueQuietly(0, -intercept.getParameterValue(0));
            for (int i = 0; i < coeff.getDimension(); i++) {
                coeff.setParameterValueQuietly(i, -coeff.getParameterValue(i));
            }
        }
    }

    @Override
    public Parameter getRateParameter() {
        return rateParameter;
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    public SquaredCachedSplines getSplines() { return splines; }
    public Parameter getEpsilon()                    { return epsilon; }

    // -----------------------------------------------------------------------
    // Citations
    // -----------------------------------------------------------------------

    @Override
    public String getDescription() {
        return "Squared-spline continuous-time coalescent likelihood";
    }

    @Override
    public List<Citation> getCitations() {
        return Arrays.asList(
                new Citation(
                        new Author[]{
                                new Author("F", "Monti"),
                                new Author("MA", "Suchard")
                        },
                        Citation.Status.IN_PREPARATION
                )
        );
    }

}
