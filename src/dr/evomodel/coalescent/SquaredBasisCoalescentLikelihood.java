/*
 * SquaredBasisCoalescentLikelihood.java
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
import dr.evomodel.coalescent.basis.IntegratedBasisExpansion;
import dr.evomodel.coalescent.timeline.BasisContext;
import dr.evomodel.coalescent.timeline.CoalescentSegmentProvider;
import dr.evomodel.coalescent.timeline.SegmentConsumer;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.*;
import dr.math.matrix.DenseSymmetricMatrixAccumulator;
import dr.math.matrix.SymmetricMatrixAccumulator;
import dr.math.matrix.SymmetricMatrixAccumulatorFactory;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.Arrays;
import java.util.List;

/**
 * General squared-basis coalescent likelihood.
 *
 * Rate function: {@code lambda(t) = epsilon + eta(t)^2},
 * where {@code eta(t) = gamma' psi(t)} is a linear basis expansion.
 *
 * The likelihood is decomposed into:
 * <ul>
 *   <li>Survival: {@code - gamma' G gamma - epsilon * c_eps},
 *       where {@code G = sum_r A_r * integral_{a_r}^{b_r} psi psi' dt} (tree-dependent,
 *       gamma-independent) and {@code c_eps = sum_r A_r * (b_r - a_r)}.</li>
 *   <li>Event terms: {@code sum_r log lambda(tau_r) + log ploidy}.</li>
 * </ul>
 *
 * G is cached and rebuilt only when the tree, basis shape, covariate breakpoints,
 * or ploidy factors change (§13 of the refactor plan).
 *
 * Parameter layout: canonical {@code [intercept, base_0, ..., base_{p-1}, covariates...]}
 * matching {@code basis.getParameter()}.  No reordering anywhere.
 *
 * @author Filippo Monti
 */
public final class SquaredBasisCoalescentLikelihood
        extends AbstractModelLikelihood
        implements HessianWrtParameterProvider, Citable, Reportable {

    private final List<BigFastTreeIntervals> intervalsList;
    private final Parameter ploidyFactors;
    private final Parameter epsilon;

    private final IntegratedBasisExpansion basis;
    private final CoalescentSegmentProvider segmentProvider;
    private final SymmetricMatrixAccumulatorFactory gramFactory;

    private SymmetricMatrixAccumulator weightedGramMatrix; // null = invalid
    private double weightedIntervalLength;

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown;

    private boolean gramDirty = true;
    private boolean gramChangedSinceLastStore = false;

    private final double[] basisBuffer;   // workspace for evaluateBasis
    private final double[] gammaBuffer;   // owned copy of gamma
    private boolean gammaBufferKnown;

    private final GramConsumer gramConsumer;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    public SquaredBasisCoalescentLikelihood(List<BigFastTreeIntervals> intervalsList,
                                             IntegratedBasisExpansion basis,
                                             CoalescentSegmentProvider segmentProvider,
                                             Parameter epsilon,
                                             Parameter ploidyFactors,
                                             SymmetricMatrixAccumulatorFactory gramFactory) {
        super("SquaredBasisCoalescentLikelihood");

        if (epsilon.getDimension() != 1)
            throw new IllegalArgumentException("epsilon must be a scalar");
        if (ploidyFactors.getDimension() != intervalsList.size())
            throw new IllegalArgumentException(
                    "ploidyFactors dimension must equal number of trees");
        if (basis.getDimension() != basis.getParameter().getDimension())
            throw new IllegalArgumentException("basis dimension mismatch");

        this.intervalsList   = intervalsList;
        this.basis           = basis;
        this.segmentProvider = segmentProvider;
        this.epsilon         = epsilon;
        this.ploidyFactors   = ploidyFactors;
        this.gramFactory     = gramFactory;

        this.basisBuffer     = new double[basis.getDimension()];
        this.gammaBuffer     = new double[basis.getDimension()];
        this.gramConsumer    = new GramConsumer();

        // Listener wiring (§6.1): coefficient changes → handleVariableChangedEvent (no G invalidation)
        addVariable(basis.getParameter());
        addModel(basis);
        addModel(segmentProvider);
        addVariable(epsilon);
        addVariable(ploidyFactors);

        for (BigFastTreeIntervals intervals : intervalsList) {
            addModel(intervals);
        }
    }

    /** Convenience constructor using the default dense accumulator. */
    public SquaredBasisCoalescentLikelihood(List<BigFastTreeIntervals> intervalsList,
                                             IntegratedBasisExpansion basis,
                                             CoalescentSegmentProvider segmentProvider,
                                             Parameter epsilon,
                                             Parameter ploidyFactors) {
        this(intervalsList, basis, segmentProvider, epsilon, ploidyFactors,
             DenseSymmetricMatrixAccumulator::new);
    }

    // -----------------------------------------------------------------------
    // Inner GramConsumer — avoids per-interval lambda allocation
    // -----------------------------------------------------------------------

    private final class GramConsumer implements SegmentConsumer {
        double weight;

        @Override
        public void accept(double start, double end, BasisContext context) {
            weightedIntervalLength += weight * (end - start);
            basis.addWeightedGramMatrix(start, end, weight, context, weightedGramMatrix);
        }
    }

    // -----------------------------------------------------------------------
    // Gram matrix — lazy build, invalidated on tree/basis/covariate changes
    // -----------------------------------------------------------------------

    private void ensureWeightedGramMatrix() {
        if (!gramDirty && weightedGramMatrix != null) return;

        weightedGramMatrix     = gramFactory.create(basis.getDimension());
        weightedIntervalLength = 0.0;

        for (int treeIdx = 0; treeIdx < intervalsList.size(); treeIdx++) {
            BigFastTreeIntervals intervals = intervalsList.get(treeIdx);
            double ploidy = ploidyFactors.getParameterValue(treeIdx);

            for (int i = 0; i < intervals.getIntervalCount(); i++) {
                int n = intervals.getLineageCount(i);
                if (n < 2) continue;
                double start  = intervals.getIntervalTime(i);
                double end    = start + intervals.getInterval(i);
                gramConsumer.weight = 0.5 * n * (n - 1) * ploidy;
                segmentProvider.forEachSegment(treeIdx, start, end, gramConsumer);
            }
        }

        gramDirty = false;
    }

    private void invalidateGram() {
        gramDirty = true;
        gramChangedSinceLastStore = true;
    }

    // -----------------------------------------------------------------------
    // Gamma buffer — lazy, re-read after any parameter or restore event
    // -----------------------------------------------------------------------

    private double[] ensureGammaBuffer() {
        if (!gammaBufferKnown) {
            Parameter param = basis.getParameter();
            int dim = param.getDimension();
            for (int i = 0; i < dim; i++) {
                gammaBuffer[i] = param.getParameterValue(i);
            }
            gammaBufferKnown = true;
        }
        return gammaBuffer;
    }

    // -----------------------------------------------------------------------
    // Log-likelihood
    // -----------------------------------------------------------------------

    @Override
    public final double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood   = calculateLogCoalescentLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    private double calculateLogCoalescentLikelihood() {
        ensureWeightedGramMatrix();

        double[] gamma = ensureGammaBuffer();
        double eps     = epsilon.getParameterValue(0);

        double logL = -weightedGramMatrix.quadraticForm(gamma)
                      - eps * weightedIntervalLength;

        for (int treeIdx = 0; treeIdx < intervalsList.size(); treeIdx++) {
            BigFastTreeIntervals intervals = intervalsList.get(treeIdx);
            double ploidy = ploidyFactors.getParameterValue(treeIdx);

            for (int i = 0; i < intervals.getIntervalCount(); i++) {
                if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                    double t = intervals.getIntervalTime(i) + intervals.getInterval(i);
                    logL += Math.log(evaluateRateAt(treeIdx, t)) + Math.log(ploidy);
                }
            }
        }
        return logL;
    }

    private double evaluateRateAt(int treeIdx, double t) {
        BasisContext ctx = segmentProvider.getEventContext(treeIdx, t);
        double eta = basis.evaluateScalar(t, ctx);
        return epsilon.getParameterValue(0) + eta * eta;
    }

    // -----------------------------------------------------------------------
    // Gradient (canonical layout, no reordering)
    // -----------------------------------------------------------------------

    @Override
    public double[] getGradientLogDensity() {
        ensureWeightedGramMatrix();

        double[] gamma = ensureGammaBuffer();
        int dim        = basis.getDimension();
        double[] grad  = new double[dim];

        // Survival term: -2 G gamma
        weightedGramMatrix.addMatVec(gamma, -2.0, grad);

        // Event terms: sum_r (2 eta_r / lambda_r) psi_r
        for (int treeIdx = 0; treeIdx < intervalsList.size(); treeIdx++) {
            BigFastTreeIntervals intervals = intervalsList.get(treeIdx);
            for (int i = 0; i < intervals.getIntervalCount(); i++) {
                if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                    double t = intervals.getIntervalTime(i) + intervals.getInterval(i);
                    BasisContext ctx = segmentProvider.getEventContext(treeIdx, t);
                    basis.evaluateBasis(t, ctx, basisBuffer);

                    double eta    = dot(gamma, basisBuffer);
                    double lambda = epsilon.getParameterValue(0) + eta * eta;
                    double scale  = 2.0 * eta / lambda;

                    for (int k = 0; k < dim; k++) {
                        grad[k] += scale * basisBuffer[k];
                    }
                }
            }
        }
        return grad;
    }

    // -----------------------------------------------------------------------
    // Hessian
    // -----------------------------------------------------------------------

    @Override
    public double[][] getHessianLogDensity() {
        ensureWeightedGramMatrix();

        int dim        = basis.getDimension();
        double[] gamma = ensureGammaBuffer();
        double[][] H   = new double[dim][dim];

        // Survival term: -2 G
        weightedGramMatrix.addScaledToDense(-2.0, H);

        // Event terms: sum_r (2/lambda_r - 4*eta_r^2/lambda_r^2) psi_r psi_r'
        for (int treeIdx = 0; treeIdx < intervalsList.size(); treeIdx++) {
            BigFastTreeIntervals intervals = intervalsList.get(treeIdx);
            for (int i = 0; i < intervals.getIntervalCount(); i++) {
                if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                    double t = intervals.getIntervalTime(i) + intervals.getInterval(i);
                    BasisContext ctx = segmentProvider.getEventContext(treeIdx, t);
                    basis.evaluateBasis(t, ctx, basisBuffer);

                    double eta    = dot(gamma, basisBuffer);
                    double lambda = epsilon.getParameterValue(0) + eta * eta;
                    double coeff  = 2.0 / lambda - 4.0 * eta * eta / (lambda * lambda);

                    for (int r = 0; r < dim; r++) {
                        for (int c = 0; c <= r; c++) {
                            double v = coeff * basisBuffer[r] * basisBuffer[c];
                            H[r][c] += v;
                            if (r != c) H[c][r] += v;
                        }
                    }
                }
            }
        }
        return H;
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {
        double[][] H = getHessianLogDensity();
        double[] diag = new double[H.length];
        for (int i = 0; i < H.length; i++) diag[i] = H[i][i];
        return diag;
    }

    // -----------------------------------------------------------------------
    // HessianWrtParameterProvider / GradientWrtParameterProvider
    // -----------------------------------------------------------------------

    @Override public Parameter getParameter()  { return basis.getParameter(); }
    @Override public int getDimension()         { return basis.getParameter().getDimension(); }
    @Override public Likelihood getLikelihood() { return this; }

    // -----------------------------------------------------------------------
    // Model event handling (§13 of the refactor plan)
    // -----------------------------------------------------------------------

    @Override public Model getModel()  { return this; }
    @Override public void makeDirty()  { likelihoodKnown = false; }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index,
                                              Parameter.ChangeType type) {
        likelihoodKnown = false;

        if (variable == basis.getParameter()) {
            gammaBufferKnown = false;
            return; // coefficient change: do NOT invalidate G
        }
        if (variable == epsilon) {
            return; // epsilon does not affect G
        }
        if (variable == ploidyFactors) {
            invalidateGram(); // A_r includes ploidy; G = sum A_r M_r must be rebuilt
            return;
        }
        throw new RuntimeException(
                "Unexpected variable event: " + variable.getVariableName());
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;

        if (model instanceof BigFastTreeIntervals) {
            invalidateGram();
            return;
        }
        if (model == basis || model == segmentProvider) {
            // Basis-shape or covariate-grid change; G must be rebuilt.
            invalidateGram();
            return;
        }
        throw new RuntimeException(
                "Unexpected model event from: " + model.getModelName());
    }

    // -----------------------------------------------------------------------
    // MCMC store / restore / accept
    // -----------------------------------------------------------------------

    @Override
    protected void storeState() {
        storedLogLikelihood       = logLikelihood;
        storedLikelihoodKnown     = likelihoodKnown;
        gramChangedSinceLastStore = false;
    }

    @Override
    protected void restoreState() {
        logLikelihood    = storedLogLikelihood;
        likelihoodKnown  = storedLikelihoodKnown;
        gammaBufferKnown = false; // parameter values may have reverted; cheap re-read

        if (gramChangedSinceLastStore) {
            gramDirty = true;
            gramChangedSinceLastStore = false;
        }
    }

    @Override
    protected void acceptState() {
        canonicalizeSign();
    }

    /**
     * Sign canonicalization: (gamma' psi)^2 = ((-gamma)' psi)^2, so the posterior
     * has a sign symmetry.  Break it by flipping all of gamma when the
     * largest-magnitude component is negative.
     *
     * Must act on the FULL basis.getParameter() in canonical layout — including
     * covariate coefficients when CovariateAugmentedBasisExpansion is used.
     */
    private void canonicalizeSign() {
        Parameter param = basis.getParameter();
        int dim = param.getDimension();
        int anchor = 0;
        double maxAbs = 0.0;
        for (int i = 0; i < dim; i++) {
            double abs = Math.abs(param.getParameterValue(i));
            if (abs > maxAbs) { maxAbs = abs; anchor = i; }
        }
        if (param.getParameterValue(anchor) < 0.0) {
            for (int i = 0; i < dim; i++) {
                param.setParameterValueQuietly(i, -param.getParameterValue(i));
            }
            gammaBufferKnown = false;
            likelihoodKnown  = false;
            param.fireParameterChangedEvent();
        }
    }

    // -----------------------------------------------------------------------
    // Utility
    // -----------------------------------------------------------------------

    private static double dot(double[] a, double[] b) {
        double sum = 0.0;
        for (int i = 0; i < a.length; i++) sum += a[i] * b[i];
        return sum;
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    public List<BigFastTreeIntervals> getIntervalsList() { return intervalsList; }
    public Parameter getPloidyFactors()                  { return ploidyFactors; }
    public Parameter getEpsilon()                        { return epsilon; }
    public IntegratedBasisExpansion getBasis()           { return basis; }
    public CoalescentSegmentProvider getSegmentProvider() { return segmentProvider; }

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
    public String getDescription() { return "Squared-basis coalescent likelihood"; }

    @Override
    public List<Citation> getCitations() {
        return Arrays.asList(new Citation(
                new Author[]{ new Author("F", "Monti"), new Author("MA", "Suchard") },
                Citation.Status.IN_PREPARATION));
    }
}
