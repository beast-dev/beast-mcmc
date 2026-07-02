/*
 * ExactTimeVaryingCoalescentLikelihood.java
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

import dr.evolution.coalescent.IntervalType;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.*;
import dr.math.BSplineGramMatrix;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.Collections;
import java.util.List;

/**
 * Abstract base class for coalescent likelihoods driven by a squared-spline rate function.
 *
 * For the rate λ(t) = ε + f(t)² with f(t) = u'φ(t), the coalescent log-likelihood is
 *
 *   log L(u) = Σ_r d_r log λ(τ_r)  −  u'Gu  −  ε · c_ε
 *
 * where
 *   G       = Σ_r A_r M(a_r, b_r)      weighted Gram matrix (tree-dependent, θ-independent)
 *   c_ε     = Σ_r A_r (b_r − a_r)      scalar for the ε survival term
 *   A_r     = C(n_r,2) · ρ             lineage-pairs-at-risk × ploidy
 *   M(a,b)  = ∫_a^b φ(t)φ(t)' dt      per-interval Gram matrix
 *
 * G is constant within one HMC trajectory (it depends only on the tree, not on u or ε).
 * It is computed lazily on first use and invalidated whenever the tree changes.
 *
 * Gradient:
 *   ∇_u log L = Σ_r d_r (2η_r/λ_r) φ_r  −  2Gu
 *
 * Hessian:
 *   ∇²_u log L = Σ_r d_r (2/λ_r − 4η_r²/λ_r²) φ_r φ_r'  −  2G
 *
 * Subclasses implement three primitives:
 *   accumulateToWeightedGramMatrix  — adds A_r M(a,b) into G for one interval
 *   evaluateRateAt(t)               — λ(t) = ε + f(t)²
 *   addEventGradient(t, grad)       — ∇_u log λ(t) accumulated into grad
 *   addEventHessian(t, H)           — ∇²_u log λ(t) accumulated into H
 *   getRateParamBuffer()            — returns u = [intercept, θ...]
 *   getEpsilonValue()               — returns ε
 *
 * @author Filippo Monti
 * @author Marc A. Suchard
 */
public abstract class ExactTimeVaryingCoalescentLikelihood extends AbstractModelLikelihood
        implements HessianWrtParameterProvider, Citable, Reportable {

    private final List<BigFastTreeIntervals> intervalsList;
    private final Parameter ploidyFactors;
    private final int augmentedDim;  // = getCoefficientDim() + 1 = size of u

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown;

    // G = Σ_r A_r M_r — null when invalid (tree changed).
    // Constant within one HMC trajectory; invalidated only by tree changes, not θ/ε changes.
    private double[][] weightedGramMatrix;
    private double     weightedIntervalLength;  // c_ε = Σ_r A_r (b_r - a_r)

    // True when the tree changed since the last storeState().
    // G is invalidated in restoreState() only when this is true — a rejected tree
    // proposal leaves G stale for the proposed tree.  For θ/ε restores G is still valid.
    private boolean treeChangedSinceLastStore = false;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    public ExactTimeVaryingCoalescentLikelihood(String name,
                                                List<BigFastTreeIntervals> intervalsList,
                                                Parameter ploidyFactors,
                                                int augmentedDim) {
        super(name);
        this.intervalsList  = intervalsList;
        this.ploidyFactors  = ploidyFactors;
        this.augmentedDim   = augmentedDim;
        this.likelihoodKnown = false;

        if (ploidyFactors.getDimension() != intervalsList.size())
            throw new IllegalArgumentException(
                    "ploidyFactors dimension must equal number of trees");

        for (BigFastTreeIntervals intervals : intervalsList)
            addModel(intervals);
        addVariable(ploidyFactors);
    }

    // -----------------------------------------------------------------------
    // Abstract template methods — implement in each subclass
    // -----------------------------------------------------------------------

    /**
     * Accumulates A_r * M(start, end) into G for the interval [start, end]
     * with weight A_r = C(n,2) * ploidy.
     * Called once per interval per tree change (when building G).
     */
    protected abstract void accumulateToWeightedGramMatrix(
            double start, double end, double weight, double[][] G);

    /**
     * Returns λ(t) = ε + f(t)² — the coalescent rate at event time t.
     */
    protected abstract double evaluateRateAt(double t);

    /**
     * Accumulates ∇_u log λ(t) = (2f(t)/λ(t)) · φ(t) into grad.
     * grad is in rateParameter layout [θ_0,...,θ_{dim-1}, intercept].
     */
    protected abstract void addEventGradient(double t, double[] grad);

    /**
     * Accumulates ∇²_u log λ(t) = (2/λ − 4η²/λ²) φ φ' into H.
     * H is in u-layout [intercept, θ_0,...,θ_{dim-1}] × [intercept, θ...].
     */
    protected abstract void addEventHessian(double t, double[][] H);

    /**
     * Returns the current u = [intercept, θ_0, ..., θ_{dim-1}] buffer.
     * Zero-allocation; caller must not retain the reference across parameter changes.
     */
    protected abstract double[] getRateParamBuffer();

    /** Returns ε — the positivity floor. */
    protected abstract double getEpsilonValue();

    /** Returns the parameter HMC and gradient providers operate on. */
    public abstract Parameter getRateParameter();

    /**
     * Called when tree topology or node times change.
     * Default is a no-op; subclasses may override for additional invalidation.
     * G is already invalidated by the base class before this is called.
     */
    protected void onTreeChanged() { }

    // -----------------------------------------------------------------------
    // Weighted Gram matrix — lazy build, invalidated on tree change only
    // -----------------------------------------------------------------------

    private void ensureWeightedGramMatrix() {
        if (weightedGramMatrix != null) return;

        weightedGramMatrix     = new double[augmentedDim][augmentedDim];
        weightedIntervalLength = 0.0;

        for (int treeIdx = 0; treeIdx < intervalsList.size(); treeIdx++) {
            BigFastTreeIntervals intervals = intervalsList.get(treeIdx);
            double ploidy = ploidyFactors.getParameterValue(treeIdx);

            for (int i = 0; i < intervals.getIntervalCount(); i++) {
                int n = intervals.getLineageCount(i);
                if (n < 2) continue;
                double weight = 0.5 * n * (n - 1) * ploidy;
                double start  = intervals.getIntervalTime(i);
                double end    = start + intervals.getInterval(i);
                weightedIntervalLength += weight * (end - start);
                accumulateToWeightedGramMatrix(start, end, weight, weightedGramMatrix);
            }
        }
    }

    private void invalidateWeightedGramMatrix() {
        weightedGramMatrix = null;
    }

    // -----------------------------------------------------------------------
    // Likelihood
    // -----------------------------------------------------------------------

    @Override
    public final double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogCoalescentLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    protected double calculateLogCoalescentLikelihood() {
        ensureWeightedGramMatrix();
        double[] u   = getRateParamBuffer();
        double   eps = getEpsilonValue();

        // Survival — quadratic term (u'Gu) and epsilon term
        double logL = -BSplineGramMatrix.symmetricQuadraticForm(weightedGramMatrix, u)
                      - eps * weightedIntervalLength;

        // Event terms
        for (int treeIdx = 0; treeIdx < intervalsList.size(); treeIdx++) {
            BigFastTreeIntervals intervals = intervalsList.get(treeIdx);
            double ploidy = ploidyFactors.getParameterValue(treeIdx);

            for (int i = 0; i < intervals.getIntervalCount(); i++) {
                if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                    double t = intervals.getIntervalTime(i) + intervals.getInterval(i);
                    logL += Math.log(evaluateRateAt(t)) + Math.log(ploidy);
                }
            }
        }
        return logL;
    }

    // -----------------------------------------------------------------------
    // Gradient
    // -----------------------------------------------------------------------

    @Override
    public double[] getGradientLogDensity() {
        ensureWeightedGramMatrix();
        double[] u    = getRateParamBuffer();
        double[] grad = new double[getRateParameter().getDimension()];

        // Survival gradient: -2 G u  (reordered from u-layout to rateParameter layout)
        BSplineGramMatrix.addMatVecReordered(weightedGramMatrix, u, -2.0, grad);

        // Event gradient terms: Σ_r d_r * (2η_r/λ_r) φ_r
        for (int treeIdx = 0; treeIdx < intervalsList.size(); treeIdx++) {
            BigFastTreeIntervals intervals = intervalsList.get(treeIdx);
            for (int i = 0; i < intervals.getIntervalCount(); i++) {
                if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                    double t = intervals.getIntervalTime(i) + intervals.getInterval(i);
                    addEventGradient(t, grad);
                }
            }
        }
        return grad;
    }

    // -----------------------------------------------------------------------
    // HessianWrtParameterProvider — diagonal and full Hessian in rateParameter layout
    // -----------------------------------------------------------------------

    /**
     * Diagonal of ∇²_u log L in rateParameter layout [θ_0,...,θ_{dim-1}, intercept].
     * Used by DiagonalHessianPreconditioning to set the HMC mass matrix.
     */
    @Override
    public double[] getDiagonalHessianLogDensity() {
        double[][] H = getHessian();  // u-layout: [intercept, θ_0,...,θ_{dim-1}]
        int n = H.length;             // = augmentedDim
        int coeffDim = n - 1;
        double[] diag = new double[n];
        for (int i = 0; i < coeffDim; i++) diag[i] = H[i + 1][i + 1];  // θ_i → rateParam[i]
        diag[coeffDim] = H[0][0];                                          // intercept → rateParam[dim]
        return diag;
    }

    /**
     * Full ∇²_u log L in rateParameter layout [θ_0,...,θ_{dim-1}, intercept].
     * Used by FullHessianPreconditioning.
     */
    @Override
    public double[][] getHessianLogDensity() {
        double[][] H = getHessian();  // u-layout
        int n = H.length;
        int coeffDim = n - 1;
        double[][] result = new double[n][n];
        for (int i = 0; i < n; i++) {
            int ui = (i < coeffDim) ? i + 1 : 0;  // u-index for rateParam-index i
            for (int j = 0; j < n; j++) {
                int uj = (j < coeffDim) ? j + 1 : 0;
                result[i][j] = H[ui][uj];
            }
        }
        return result;
    }

    // -----------------------------------------------------------------------
    // Internal Hessian — in u-layout [intercept, θ...]; reordered by the two methods above
    // -----------------------------------------------------------------------

    /**
     * Returns the (augmented-dim × augmented-dim) Hessian in u-layout.
     *
     * ∇²_u log L = -2G + Σ_r d_r (2/λ_r - 4η_r²/λ_r²) φ_r φ_r'
     *
     * Use getDiagonalHessianLogDensity() / getHessianLogDensity() for the rateParameter layout.
     */
    public double[][] getHessian() {
        ensureWeightedGramMatrix();
        int n = augmentedDim;
        double[][] H = new double[n][n];

        // Survival: -2G
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                H[i][j] = -2.0 * weightedGramMatrix[i][j];

        // Event terms: Σ d_r (2/λ_r - 4η_r²/λ_r²) φ_r φ_r'
        for (int treeIdx = 0; treeIdx < intervalsList.size(); treeIdx++) {
            BigFastTreeIntervals intervals = intervalsList.get(treeIdx);
            for (int i = 0; i < intervals.getIntervalCount(); i++) {
                if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                    double t = intervals.getIntervalTime(i) + intervals.getInterval(i);
                    addEventHessian(t, H);
                }
            }
        }
        return H;
    }

    // -----------------------------------------------------------------------
    // GradientWrtParameterProvider
    // -----------------------------------------------------------------------

    @Override public final Parameter getParameter()  { return getRateParameter(); }
    @Override public final int getDimension()         { return getRateParameter().getDimension(); }
    @Override public final Likelihood getLikelihood() { return this; }

    // -----------------------------------------------------------------------
    // Model event handling
    // -----------------------------------------------------------------------

    @Override public final Model getModel()  { return this; }
    @Override public final void makeDirty()  { likelihoodKnown = false; }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model instanceof BigFastTreeIntervals) {
            invalidateWeightedGramMatrix();
            treeChangedSinceLastStore = true;
            onTreeChanged();
            likelihoodKnown = false;
        } else {
            throw new RuntimeException("Unknown model: " + model);
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // θ/ε changes dirty the likelihood but NOT G — G is tree-dependent only.
        likelihoodKnown = false;
    }

    // -----------------------------------------------------------------------
    // MCMC store / restore / accept
    // -----------------------------------------------------------------------

    @Override
    protected void storeState() {
        storedLogLikelihood       = logLikelihood;
        treeChangedSinceLastStore = false;  // reset: track changes from this point
    }

    @Override
    protected void restoreState() {
        logLikelihood   = storedLogLikelihood;
        likelihoodKnown = true;
        if (treeChangedSinceLastStore) {
            // A tree proposal was rejected: intervals restored to pre-proposal values
            // but G was built for the proposed tree.  Invalidate so any forced
            // recomputation (BEAST correctness check) uses the correct endpoints.
            invalidateWeightedGramMatrix();
            treeChangedSinceLastStore = false;
        }
        // If only θ/ε changed (e.g. HMC reject): G is still valid — tree unchanged.
    }

    @Override
    protected void acceptState() { }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    public List<BigFastTreeIntervals> getIntervalsList() { return intervalsList; }
    public Parameter getPloidyFactors()                  { return ploidyFactors; }

    // -----------------------------------------------------------------------
    // Citable / Reportable
    // -----------------------------------------------------------------------

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(
                this, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, null);
    }

    @Override public Citation.Category getCategory() { return Citation.Category.TREE_PRIORS; }

    @Override
    public String getDescription() { return "Exact time-varying coalescent likelihood"; }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(
                new Citation(
                        new Author[]{ new Author("F", "Monti"), new Author("MA", "Suchard") },
                        Citation.Status.IN_PREPARATION));
    }
}
