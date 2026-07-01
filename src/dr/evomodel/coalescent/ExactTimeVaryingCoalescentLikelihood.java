/*
 * ExactCoalescentLikelihood.java
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
import dr.inference.model.*;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.Collections;
import java.util.List;

/**
 * Abstract base class for coalescent likelihoods driven by a continuous-time rate function.
 *
 * The log-likelihood is
 *
 *   log L = sum_trees sum_intervals
 *              -C(n,2) * rho * integratedHazard(start, end)
 *            + [if COALESCENT event] logHazard(end) + log(rho)
 *
 * where C(n,2) = n(n-1)/2, rho is the ploidy factor, and the rate function 1/N(t) is
 * represented by the two abstract primitives integratedHazard and logHazard.
 *
 * The gradient w.r.t. the rate-function parameters is assembled from two matching abstract
 * primitives: addIntegratedHazardGradient and addLogHazardGradient.
 *
 * Subclasses need only implement these four methods plus getRateParameter() and onTreeChanged().
 * The interval loop, ploidy handling, multilocus support, gradient loop, model event wiring,
 * store/restore, and GradientWrtParameterProvider plumbing are all provided here.
 *
 * @author Filippo Monti
 */
public abstract class ExactTimeVaryingCoalescentLikelihood extends AbstractModelLikelihood
        implements GradientWrtParameterProvider, Citable, Reportable {

    private final List<BigFastTreeIntervals> intervalsList;
    private final Parameter ploidyFactors;

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    public ExactTimeVaryingCoalescentLikelihood(String name,
                                                List<BigFastTreeIntervals> intervalsList,
                                                Parameter ploidyFactors) {
        super(name);

        this.intervalsList  = intervalsList;
        this.ploidyFactors  = ploidyFactors;
        this.likelihoodKnown = false;

        if (ploidyFactors.getDimension() != intervalsList.size()) {
            throw new IllegalArgumentException(
                    "ploidyFactors dimension (" + ploidyFactors.getDimension() +
                    ") must equal number of trees (" + intervalsList.size() + ")");
        }

        for (BigFastTreeIntervals intervals : intervalsList) {
            addModel(intervals);
        }
        addVariable(ploidyFactors);
    }

    // -----------------------------------------------------------------------
    // Abstract template methods — implement in each subclass
    // -----------------------------------------------------------------------

    /**
     * Returns ∫_start^end λ(t) dt — the integrated hazard over an interval.
     * globalIdx is the sequential interval index across all trees in intervalsList
     * (tree 0: 0..n0-1, tree 1: n0..n0+n1-1, …).  Used by subclasses for index-based caching.
     */
    protected abstract double integratedHazard(int globalIdx, double start, double end);

    /**
     * Returns log λ(t) — the log hazard at a coalescent event time.
     * Called once per coalescent event per likelihood evaluation.
     */
    protected abstract double logHazard(double t);

    /**
     * Accumulates coef * ∇_θ ∫_start^end λ(t) dt into grad.
     * coef = -C(n,2) * ρ (already negative).
     * globalIdx: same sequential index as integratedHazard — for index-based gram cache.
     */
    protected abstract void addIntegratedHazardGradient(
            int globalIdx, double start, double end, double coef, double[] grad);

    /**
     * Accumulates ∇_θ log λ(t) into grad (one unit per coalescent event; no C(n,2) factor).
     */
    protected abstract void addLogHazardGradient(double t, double[] grad);

    /**
     * Returns the parameter that HMC / gradient operators act on.
     */
    public abstract Parameter getRateParameter();

    /**
     * Called when tree topology or node times change.
     * Subclasses should clear any interval-keyed caches here (e.g. gram matrix cache).
     */
    protected abstract void onTreeChanged();

    // -----------------------------------------------------------------------
    // Likelihood loop (concrete)
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
        double logL = 0.0;
        int globalIdx = 0;

        for (int treeIdx = 0; treeIdx < intervalsList.size(); treeIdx++) {
            BigFastTreeIntervals intervals = intervalsList.get(treeIdx);
            double ploidy = ploidyFactors.getParameterValue(treeIdx);

            for (int i = 0; i < intervals.getIntervalCount(); i++) {
                double start = intervals.getIntervalTime(i);
                double end   = start + intervals.getInterval(i);
                int    n     = intervals.getLineageCount(i);

                if (n > 1) {
                    double coef = -0.5 * n * (n - 1) * ploidy;
                    logL += coef * integratedHazard(globalIdx, start, end);
                }

                if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                    logL += logHazard(end) + Math.log(ploidy);
                }

                globalIdx++;
            }
        }

        return logL;
    }

    // -----------------------------------------------------------------------
    // Gradient loop (concrete)
    // -----------------------------------------------------------------------

    @Override
    public double[] getGradientLogDensity() {
        double[] grad = new double[getRateParameter().getDimension()];
        int globalIdx = 0;

        for (int treeIdx = 0; treeIdx < intervalsList.size(); treeIdx++) {
            BigFastTreeIntervals intervals = intervalsList.get(treeIdx);
            double ploidy = ploidyFactors.getParameterValue(treeIdx);

            for (int i = 0; i < intervals.getIntervalCount(); i++) {
                double start = intervals.getIntervalTime(i);
                double end   = start + intervals.getInterval(i);
                int    n     = intervals.getLineageCount(i);

                if (n > 1) {
                    double coef = -0.5 * n * (n - 1) * ploidy;
                    addIntegratedHazardGradient(globalIdx, start, end, coef, grad);
                }

                if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                    addLogHazardGradient(end, grad);
                }

                globalIdx++;
            }
        }

        return grad;
    }

    // -----------------------------------------------------------------------
    // GradientWrtParameterProvider (concrete)
    // -----------------------------------------------------------------------

    @Override
    public final Parameter getParameter() {
        return getRateParameter();
    }

    @Override
    public final int getDimension() {
        return getRateParameter().getDimension();
    }

    @Override
    public final Likelihood getLikelihood() {
        return this;
    }

    // -----------------------------------------------------------------------
    // Model event handling
    // -----------------------------------------------------------------------

    @Override
    public final Model getModel() { return this; }

    @Override
    public final void makeDirty() { likelihoodKnown = false; }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model instanceof BigFastTreeIntervals) {
            onTreeChanged();
            likelihoodKnown = false;
        } else {
            throw new RuntimeException("Unknown model: " + model);
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Any parameter change (ploidy, rate params, epsilon, …) dirties the likelihood.
        // Subclasses do NOT need to clear the gram cache here — M(a,b) is parameter-independent.
        likelihoodKnown = false;
    }

    // -----------------------------------------------------------------------
    // MCMC store / restore / accept
    // -----------------------------------------------------------------------

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        logLikelihood  = storedLogLikelihood;
        likelihoodKnown = true;
    }

    @Override
    protected void acceptState() { }

    // -----------------------------------------------------------------------
    // Accessors for subclasses
    // -----------------------------------------------------------------------

    public List<BigFastTreeIntervals> getIntervalsList() { return intervalsList; }
    public Parameter getPloidyFactors()                  { return ploidyFactors; }

    // -----------------------------------------------------------------------
    // Citable / Reportable
    // -----------------------------------------------------------------------

    @Override
    public String getReport() {
        // Produces "analytic: [...]\nnumeric : [...]" for use in TestXML assertEqual blocks.
        return GradientWrtParameterProvider.getReportAndCheckForError(
                this, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, null);
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Continuous-time coalescent likelihood";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(
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
