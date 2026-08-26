/*
 * MascotCore.java
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
 */

package dr.evomodel.coalescent.mascot;

import dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation;

import java.util.Arrays;
import java.util.List;

/**
 * Allocation-light engine for a MASCOT-style marginal structured coalescent
 * likelihood and the reverse-mode derivative of the same RK4 discretization.
 *
 * The flat parameter layout is per epoch:
 *
 * <pre>
 * m[0,1], m[0,2], ..., m[K-1,K-2], log N[0], ..., log N[K-1]
 * </pre>
 *
 * where rows are source states and columns are destination states.
 *
 * A {@code MascotCore} instance owns a reusable {@link Workspace} and a
 * per-epoch rate cache, both of which are updated in place across calls to
 * {@link #evaluate}. This makes a single instance non-thread-safe: callers
 * that need concurrent evaluation must use one instance per thread.
 */
public final class MascotCore {

    // Package-private (not private): also used by MascotForwardModeHelper's event-time
    // and epoch-boundary bookkeeping.
    static final double TIME_TOLERANCE = 1.0e-14;

    private final int stateCount;
    private final double[] boundaries;
    private final double maxStep;
    private final int migrationParametersPerEpoch;
    private final int parametersPerEpoch;
    private final int parameterCount;
    private final int epochCount;

    private final Workspace workspace = new Workspace();
    private final EpochRates[] epochRates;
    private final ActiveState activeState = new ActiveState();
    // Owns the forward RK4 integration and event application (including
    // recording the OperationTapeStore this instance's reverse pass replays);
    // see its own class doc for why this is a separate top-level class rather
    // than more private methods here.
    private final MascotForwardModeHelper forwardHelper;
    // Owns the backward replay over this instance's OperationTapeStore; see its
    // own class doc for why this is a separate top-level class rather than more
    // private methods here.
    private final MascotReverseModeHelper reverseHelper;

    public MascotCore(int stateCount, double[] boundaries, double maxStep) {
        if (stateCount < 2) {
            throw new IllegalArgumentException("stateCount must be at least 2");
        }
        if (boundaries == null || boundaries.length < 2) {
            throw new IllegalArgumentException("boundaries must include at least 0.0 and infinity");
        }
        if (Math.abs(boundaries[0]) > TIME_TOLERANCE) {
            throw new IllegalArgumentException("the first boundary must be 0.0");
        }
        for (int i = 1; i < boundaries.length; i++) {
            if (!(boundaries[i] > boundaries[i - 1])) {
                throw new IllegalArgumentException("boundaries must be strictly increasing");
            }
        }
        if (!Double.isInfinite(boundaries[boundaries.length - 1])) {
            throw new IllegalArgumentException("the final boundary must be positive infinity");
        }
        if (maxStep <= 0.0 || !Double.isFinite(maxStep)) {
            throw new IllegalArgumentException("maxStep must be positive and finite");
        }

        this.stateCount = stateCount;
        this.boundaries = boundaries.clone();
        this.maxStep = maxStep;
        this.migrationParametersPerEpoch = stateCount * (stateCount - 1);
        this.parametersPerEpoch = migrationParametersPerEpoch + stateCount;
        this.epochCount = this.boundaries.length - 1;
        this.parameterCount = epochCount * parametersPerEpoch;

        this.epochRates = new EpochRates[epochCount];
        for (int epoch = 0; epoch < epochCount; epoch++) {
            epochRates[epoch] = new EpochRates(stateCount);
        }

        this.forwardHelper = new MascotForwardModeHelper(stateCount, maxStep, this.boundaries, epochCount);
        this.reverseHelper = new MascotReverseModeHelper(stateCount, parametersPerEpoch,
                migrationParametersPerEpoch, parameterCount);
    }

    public int getStateCount() {
        return stateCount;
    }

    public int getEpochCount() {
        return epochCount;
    }

    public int getParameterCount() {
        return parameterCount;
    }

    public double logLikelihood(PreparedOperations prepared, double[] theta) {
        return evaluate(prepared, theta, null, false, false, false, false).logLikelihood;
    }

    public Result likelihoodAndGradient(PreparedOperations prepared, double[] theta) {
        return evaluate(prepared, theta, null, true, false, false, false);
    }

    /**
     * Likelihood plus adjoint node-state scores for every internal node (see
     * {@code Result#ancestralStateScores}), without a parameter gradient.
     */
    public Result likelihoodAndAncestralStates(PreparedOperations prepared, double[] theta, double[] branchRates,
                                               boolean checkProbabilities) {
        return evaluate(prepared, theta, branchRates, false, true, checkProbabilities, false);
    }

    /**
     * Sorts and validates an event array once so that repeated evaluations against
     * the same fixed tree (only parameters changing) do not repeat the sort or the
     * event-array clone. The returned object is immutable and safe to share across
     * many {@link #evaluate} (or {@link #evaluateInto}/{@link #evaluateLikelihood})
     * calls.
     */
    public static PreparedOperations prepareOperations(List<BranchIntervalOperation> operations,
                                                       List<Integer> intervalStarts,
                                                       double[][] tipPartials,
                                                       int nodeCount,
                                                       double initialTime) {
        if (operations == null || operations.isEmpty()) {
            throw new IllegalArgumentException("at least one branch interval operation is required");
        }
        if (intervalStarts == null || intervalStarts.size() < 2 ||
                intervalStarts.get(0) != 0 ||
                intervalStarts.get(intervalStarts.size() - 1) != operations.size()) {
            throw new IllegalArgumentException("invalid interval starts");
        }
        if (tipPartials == null || nodeCount <= 0 || tipPartials.length > nodeCount) {
            throw new IllegalArgumentException("invalid tip partials or node count");
        }
        if (initialTime < -TIME_TOLERANCE || !Double.isFinite(initialTime)) {
            throw new IllegalArgumentException("invalid initial time: " + initialTime);
        }
        BranchIntervalOperation[] operationArray =
                operations.toArray(new BranchIntervalOperation[operations.size()]);
        int[] starts = new int[intervalStarts.size()];
        for (int i = 0; i < starts.length; i++) {
            starts[i] = intervalStarts.get(i);
            if (i > 0 && starts[i] <= starts[i - 1]) {
                throw new IllegalArgumentException("empty or unordered operation interval at index " + (i - 1));
            }
        }
        double[][] partials = new double[tipPartials.length][];
        for (int i = 0; i < tipPartials.length; i++) {
            if (tipPartials[i] == null) {
                throw new IllegalArgumentException("missing tip partials for node " + i);
            }
            partials[i] = tipPartials[i].clone();
        }
        return new PreparedOperations(operationArray, starts, partials, nodeCount, initialTime);
    }

    /**
     * The one public, allocating, {@link Result}-returning evaluation entry
     * point. {@code branchRates}, when non-null, applies an optional
     * branch-specific clock-rate multiplier on the migration/transition
     * process only (never on the Ne-derived coalescent rate), matching BASTA's
     * {@code branchRateModel.getBranchRate(tree, child) * branchTime}
     * convention; it must be indexed by lineage id (the same ids used in
     * {@link Event#getLineage()}/child1/child2/parent) and sized larger than
     * the maximum lineage id among {@code prepared}'s events (the root's slot
     * is never read). {@code computeGradient} and {@code computeAncestralStates}
     * are independent: both, either, or neither may be requested in one call
     * (see {@code Result#ancestralStateScores} and {@code
     * MASCOT_ADJOINT_ANCESTRAL_RECONSTRUCTION.md}); requesting neither costs
     * nothing beyond the forward pass (no operation tape is built, no
     * ancestral array is allocated). {@code copyFinalState=false} skips
     * copying the final root probabilities and active-lineage ids into the
     * returned {@code Result}; callers that only need {@code
     * logLikelihood}/{@code gradient} should use {@code false} to avoid two
     * allocations and copies per call that nobody reads.
     * <p/>
     * Production BEAST code (see {@link
     * dr.evomodel.coalescent.mascot.MascotLikelihood}) does not use this --
     * it calls {@link #evaluateLikelihood}/{@link #evaluateInto} directly, to
     * avoid the {@link Result} allocation entirely. This method exists for
     * tests and other callers that want a self-contained {@code Result}
     * object rather than caller-owned output buffers.
     */
    public Result evaluate(PreparedOperations prepared, double[] theta, double[] branchRates, boolean computeGradient,
                           boolean computeAncestralStates, boolean checkProbabilities, boolean copyFinalState) {
        return evaluate(prepared, theta, branchRates, computeGradient, computeAncestralStates, checkProbabilities,
                copyFinalState, null);
    }

    /**
     * Test-only hook for finite-difference validation of the adjoint
     * node-state score (see {@code MASCOT_ADJOINT_ANCESTRAL_RECONSTRUCTION.md}
     * Section 13). Applies a hypothetical per-state log weight {@code
     * nodeLogWeights[nodeId * stateCount + state]} at each coalescent event's
     * parent node before normalizing that event's state distribution; {@code
     * nodeLogWeights == null} is exactly the production forward computation.
     * Package-private: not part of the public/XML-facing API.
     */
    Result evaluateWithNodeLogWeightsForTesting(PreparedOperations prepared, double[] theta, double[] branchRates,
                                                double[] nodeLogWeights, boolean checkProbabilities) {
        return evaluate(prepared, theta, branchRates, false, false, checkProbabilities, false, nodeLogWeights);
    }

    private Result evaluate(PreparedOperations prepared, double[] theta, double[] branchRates, boolean computeGradient,
                            boolean computeAncestralStates, boolean checkProbabilities, boolean copyFinalState,
                            double[] nodeLogWeights) {
        double[] gradientOut = computeGradient ? new double[parameterCount] : null;
        double[] clockGradientOut = computeGradient && branchRates != null
                ? new double[prepared.maxLineageId + 1] : null;
        double[] ancestralStateScoresOut = computeAncestralStates
                ? new double[(prepared.maxLineageId + 1) * stateCount] : null;

        double logLikelihood = evaluateCore(prepared, theta, branchRates, checkProbabilities, nodeLogWeights,
                gradientOut, clockGradientOut, ancestralStateScoresOut);

        ActiveState state = activeState;
        double[] rootProbabilities = copyFinalState ? state.copyProbabilities() : null;
        int[] activeLineages = copyFinalState ? state.copyActiveIds() : null;
        return new Result(logLikelihood, gradientOut, clockGradientOut, rootProbabilities, activeLineages,
                ancestralStateScoresOut);
    }

    /**
     * Caller-owned-output evaluation: writes into {@code gradientOut}/{@code
     * clockGradientOut}/{@code ancestralStateScoresOut} instead of allocating
     * them, and never constructs a {@link Result}. Each output array is
     * requested by passing a non-null, exactly-sized destination; passing
     * {@code null} skips exposing (not necessarily computing -- see below)
     * that output. This is the production path {@link
     * dr.evomodel.coalescent.mascot.MascotLikelihood} uses; {@link #evaluate}
     * and friends remain allocating compatibility wrappers over the same
     * {@link #evaluateCore} implementation, for tests and callers that want a
     * self-contained {@link Result}.
     * <p/>
     * The reverse traversal computes the parameter gradient and the
     * ancestral-state scores together (shared VJP), so requesting only one of
     * the two still runs the same reverse pass; the unrequested output is
     * written into a reused scratch buffer instead of the caller's array.
     */
    public double evaluateInto(PreparedOperations prepared, double[] theta, double[] branchRates,
                               double[] gradientOut, double[] clockGradientOut, double[] ancestralStateScoresOut,
                               boolean checkProbabilities) {
        return evaluateCore(prepared, theta, branchRates, checkProbabilities, null,
                gradientOut, clockGradientOut, ancestralStateScoresOut);
    }

    /** Scalar-only evaluation: never builds the reverse tape, never allocates a {@link Result}. */
    public double evaluateLikelihood(PreparedOperations prepared, double[] theta, double[] branchRates,
                                     boolean checkProbabilities) {
        return evaluateCore(prepared, theta, branchRates, checkProbabilities, null, null, null, null);
    }

    private double evaluateCore(PreparedOperations prepared, double[] theta, double[] branchRates,
                                boolean checkProbabilities, double[] nodeLogWeights,
                                double[] gradientOut, double[] clockGradientOut, double[] ancestralStateScoresOut) {
        if (prepared == null) {
            throw new IllegalArgumentException("prepared operations must not be null");
        }
        if (theta == null || theta.length != parameterCount) {
            throw new IllegalArgumentException("theta dimension " + (theta == null ? -1 : theta.length) +
                    " does not match expected dimension " + parameterCount);
        }
        if (branchRates != null && branchRates.length <= prepared.maxLineageId) {
            throw new IllegalArgumentException("branchRates dimension " + branchRates.length +
                    " does not cover the maximum lineage id " + prepared.maxLineageId);
        }
        if (nodeLogWeights != null) {
            int expectedLength = (prepared.maxLineageId + 1) * stateCount;
            if (nodeLogWeights.length != expectedLength) {
                throw new IllegalArgumentException("nodeLogWeights dimension " + nodeLogWeights.length +
                        " does not match expected dimension " + expectedLength);
            }
        }
        if (gradientOut != null && gradientOut.length != parameterCount) {
            throw new IllegalArgumentException("gradientOut dimension " + gradientOut.length +
                    " does not match expected dimension " + parameterCount);
        }
        int expectedLineageDimension = prepared.maxLineageId + 1;
        if (clockGradientOut != null) {
            if (branchRates == null) {
                throw new IllegalArgumentException("clockGradientOut was requested but branchRates is null");
            }
            if (clockGradientOut.length != expectedLineageDimension) {
                throw new IllegalArgumentException("clockGradientOut dimension " + clockGradientOut.length +
                        " does not match expected dimension " + expectedLineageDimension);
            }
        }
        if (ancestralStateScoresOut != null && ancestralStateScoresOut.length != expectedLineageDimension * stateCount) {
            throw new IllegalArgumentException("ancestralStateScoresOut dimension " + ancestralStateScoresOut.length +
                    " does not match expected dimension " + (expectedLineageDimension * stateCount));
        }

        for (int epoch = 0; epoch < epochCount; epoch++) {
            updateEpochRates(theta, epoch, epochRates[epoch]);
        }

        ActiveState state = activeState;
        boolean runReverse = gradientOut != null || ancestralStateScoresOut != null;
        OperationTapeStore operations = null;
        if (runReverse) {
            // One tape entry per event, plus at most one IntervalTape per epoch
            // transition between events. The store owns reusable tape objects and
            // stage arrays, so fixed-tree HMC overwrites the same storage instead
            // of rebuilding the reverse tape at every parameter evaluation.
            operations = workspace.operationTapes;
            operations.reset(prepared.operations.length + epochCount);
        }

        forwardHelper.forward(state, epochRates, workspace, prepared, branchRates, checkProbabilities,
                operations, nodeLogWeights);

        if (ancestralStateScoresOut != null) {
            Arrays.fill(ancestralStateScoresOut, 0, expectedLineageDimension * stateCount, Double.NaN);
        }

        if (runReverse) {
            double[] gradientBuffer = gradientOut;
            if (gradientBuffer == null) {
                workspace.gradientScratch = ensure(workspace.gradientScratch, parameterCount);
                gradientBuffer = workspace.gradientScratch;
            }
            double[] clockGradientBuffer = clockGradientOut;
            if (clockGradientBuffer == null && branchRates != null) {
                workspace.clockGradientScratch = ensure(workspace.clockGradientScratch, expectedLineageDimension);
                clockGradientBuffer = workspace.clockGradientScratch;
            }
            reverseHelper.reverse(workspace, epochRates, operations, state.activeCount, gradientBuffer,
                    clockGradientBuffer, ancestralStateScoresOut);
        }

        return state.logLikelihood;
    }

    // ------------------------------------------------------------------
    // Epoch rate cache
    // ------------------------------------------------------------------

    private void updateEpochRates(double[] theta, int epoch, EpochRates rates) {
        int thetaOffset = epoch * parametersPerEpoch;
        int index = 0;
        for (int source = 0; source < stateCount; source++) {
            double rowSum = 0.0;
            int row = source * stateCount;
            for (int sink = 0; sink < stateCount; sink++) {
                if (source == sink) {
                    continue;
                }
                double rate = theta[thetaOffset + index];
                // >= 0, not > 0: a rate of exactly zero (e.g. a BSSVS indicator
                // switched off) is mathematically fine here -- this value is only
                // ever written directly into migrationMatrix/summed for the
                // diagonal, never logged or divided into, so nothing downstream
                // requires strict positivity. Only negative/non-finite is a bug.
                if (!(rate >= 0.0) || !Double.isFinite(rate)) {
                    throw new NumericalException("invalid migration rate for epoch " + epoch +
                            ", source " + source + ", sink " + sink + ": " + rate);
                }
                rates.migrationRates[index] = rate;
                rates.migrationMatrix[row + sink] = rate;
                rowSum += rate;
                if (!Double.isFinite(rowSum)) {
                    throw new NumericalException("invalid migration row sum for epoch " + epoch +
                            ", source " + source + ": " + rowSum);
                }
                index++;
            }
            rates.migrationMatrix[row + source] = -rowSum;
        }

        int etaOffset = thetaOffset + migrationParametersPerEpoch;
        for (int state = 0; state < stateCount; state++) {
            double logPopulation = theta[etaOffset + state];
            double inversePopulation = Math.exp(-logPopulation);
            if (!Double.isFinite(logPopulation) || !Double.isFinite(inversePopulation)) {
                throw new NumericalException("invalid log population size for epoch " + epoch +
                        ", state " + state + ": " + logPopulation);
            }
            rates.inversePopulation[state] = inversePopulation;
        }
    }

    // ------------------------------------------------------------------
    // Small array helpers
    // ------------------------------------------------------------------

    // Package-private (not private): shared verbatim with MascotReverseModeHelper,
    // which has no MascotCore instance state of its own to call back through.
    static double[] ensure(double[] array, int size) {
        return (array == null || array.length < size) ? new double[size] : array;
    }

    private static int[] ensureInt(int[] array, int size) {
        return (array == null || array.length < size) ? new int[size] : array;
    }

    // Package-private (not private): also used by MascotForwardModeHelper's RK4 steps.
    static void addScaledInto(double[] x, double[] dx, double scale, double[] out, int n) {
        for (int i = 0; i < n; i++) {
            out[i] = x[i] + scale * dx[i];
        }
    }

    static void scaleInto(double[] x, double scale, double[] out, int n) {
        for (int i = 0; i < n; i++) {
            out[i] = scale * x[i];
        }
    }

    static void addInPlace(double[] destination, double[] source, int n) {
        for (int i = 0; i < n; i++) {
            destination[i] += source[i];
        }
    }

    static void addScaledInPlace(double[] destination, double[] source, double scale, int n) {
        for (int i = 0; i < n; i++) {
            destination[i] += scale * source[i];
        }
    }

    /**
     * An already-sorted, validated event sequence. Building one requires a full
     * clone and sort of the input array; evaluating against a {@code PreparedEvents}
     * does not. Callers that repeatedly evaluate the same fixed tree under changing
     * parameters (the common MCMC/HMC pattern) should build this once per tree
     * change and reuse it across evaluations.
     */
    public static final class PreparedOperations {
        // Package-private (not private): read directly by MascotForwardModeHelper.forward().
        final BranchIntervalOperation[] operations;
        final int[] intervalStarts;
        final double[][] tipPartials;
        final int nodeCount;
        final int maxLineageId;
        final double initialTime;

        private PreparedOperations(BranchIntervalOperation[] operations, int[] intervalStarts,
                                   double[][] tipPartials, int nodeCount, double initialTime) {
            this.operations = operations;
            this.intervalStarts = intervalStarts;
            this.tipPartials = tipPartials;
            this.nodeCount = nodeCount;
            this.maxLineageId = nodeCount - 1;
            this.initialTime = initialTime;
        }
    }

    public static final class Result {
        public final double logLikelihood;
        public final double[] gradient;
        /**
         * d(logLikelihood)/d(branchRate[lineageId]), indexed the same way as the
         * {@code branchRates} array passed into {@code evaluate(...)}. Null unless
         * both a gradient was requested and a non-null {@code branchRates} array
         * was supplied.
         */
        public final double[] clockGradient;
        public final double[] rootProbabilities;
        public final int[] activeLineages;
        /**
         * Flat, node-major adjoint node-state scores, indexed by {@code
         * nodeId * stateCount + state} (see
         * {@code MASCOT_ADJOINT_ANCESTRAL_RECONSTRUCTION.md}). {@code nodeId}
         * is the stable tree node id ({@code NodeRef.getNumber()}), not an
         * internal active-lineage array position. Tip rows and any unused
         * node-number slots are {@link Double#NaN}. Null unless ancestral
         * reconstruction was requested. These are exact sensitivities of the
         * discretized RK4 likelihood, not posterior probabilities: a
         * nonnegativity survey found a converged (step-size-independent)
         * negative score under highly unequal population sizes, so this
         * quantity is not always a valid probability distribution -- see
         * MASCOT_ADJOINT_ANCESTRAL_RECONSTRUCTION.md's "Probability
         * interpretation" note. Callers must not present these as
         * probabilities.
         */
        public final double[] ancestralStateScores;

        private Result(double logLikelihood, double[] gradient, double[] clockGradient,
                       double[] rootProbabilities, int[] activeLineages, double[] ancestralStateScores) {
            this.logLikelihood = logLikelihood;
            this.gradient = gradient;
            this.clockGradient = clockGradient;
            this.rootProbabilities = rootProbabilities;
            this.activeLineages = activeLineages;
            this.ancestralStateScores = ancestralStateScores;
        }
    }

    /**
     * Numerical failures caused by the current parameter proposal rather than by
     * malformed tree/event structure. BEAST likelihood wrappers may translate
     * this to {@code -Infinity}; gradient callers should fail loudly instead of
     * returning synthetic derivatives.
     */
    public static final class NumericalException extends RuntimeException {
        public NumericalException(String message) {
            super(message);
        }
    }

    /**
     * Reused across {@link #evaluate} calls (growth-only, like {@link Workspace})
     * rather than allocated fresh every call. Active/inactive status is tracked
     * with a generation stamp instead of a {@code -1}-filled map: {@link #reset}
     * bumps {@link #currentGeneration} instead of re-filling {@code
     * lineageGeneration}, so resetting costs O(1) instead of O(maxLineageId).
     * Class and members are package-private (not {@code private}) for the same
     * reason as {@link Workspace}'s: {@link MascotForwardModeHelper} advances
     * this state directly in its own hot loop.
     */
    static final class ActiveState {
        int[] activeIds;
        double[] probabilities;
        int activeCount;
        double logLikelihood;
        int stateCount;
        int[] lineageToActiveIndex;
        int[] lineageGeneration;
        int currentGeneration;

        private ActiveState() {
            this.activeIds = new int[1];
            this.probabilities = new double[1];
            this.lineageToActiveIndex = new int[1];
            this.lineageGeneration = new int[1];
            this.currentGeneration = 1;
        }

        void reset(int stateCount, int capacity, int maxLineageId) {
            this.stateCount = stateCount;
            if (activeIds.length < capacity) {
                activeIds = new int[capacity];
            }
            int probabilityCapacity = capacity * stateCount;
            if (probabilities.length < probabilityCapacity) {
                probabilities = new double[probabilityCapacity];
            }
            int lineageCapacity = maxLineageId + 1;
            if (lineageGeneration.length < lineageCapacity) {
                lineageToActiveIndex = new int[lineageCapacity];
                lineageGeneration = new int[lineageCapacity];
            }
            activeCount = 0;
            logLikelihood = 0.0;
            currentGeneration++;
            if (currentGeneration == Integer.MAX_VALUE) {
                Arrays.fill(lineageGeneration, 0);
                currentGeneration = 1;
            }
        }

        void ensureCapacity(int capacity) {
            if (activeIds.length < capacity) {
                activeIds = Arrays.copyOf(activeIds, capacity);
            }
            int probabilityCapacity = capacity * stateCount;
            if (probabilities.length < probabilityCapacity) {
                probabilities = Arrays.copyOf(probabilities, probabilityCapacity);
            }
        }

        boolean isActive(int lineageId) {
            return lineageId >= 0 && lineageId < lineageGeneration.length
                    && lineageGeneration[lineageId] == currentGeneration;
        }

        int activeIndexOf(int lineageId) {
            if (!isActive(lineageId)) {
                return -1;
            }
            return lineageToActiveIndex[lineageId];
        }

        void setActiveIndex(int lineageId, int index) {
            lineageGeneration[lineageId] = currentGeneration;
            lineageToActiveIndex[lineageId] = index;
        }

        void clearActiveIndex(int lineageId) {
            if (lineageId >= 0 && lineageId < lineageGeneration.length) {
                lineageGeneration[lineageId] = 0;
            }
        }

        double[] copyProbabilities() {
            return Arrays.copyOf(probabilities, activeCount * stateCount);
        }

        int[] copyActiveIds() {
            return Arrays.copyOf(activeIds, activeCount);
        }
    }

    /**
     * Per-instance reusable scratch memory. Growth-only: arrays are replaced by a
     * larger array when a call needs more capacity than is currently held, and are
     * never shrunk. Not thread-safe; see the class-level documentation.
     *
     * Fields are package-private (not {@code private}) so that MascotCore's own
     * hot-loop methods read/write them as direct field access rather than through
     * the synthetic accessor bridge methods javac must otherwise generate for
     * cross-nested-class private access; profiling showed those bridges taking
     * measurable self time (see MascotCoreProfileDriver). The class itself is
     * package-private for the same reason, one level up: {@link
     * MascotReverseModeHelper} is a separate top-level class (not nested in
     * MascotCore) that reads/writes these fields directly in its own hot loop.
     */
    static final class Workspace {
        double[] integrationState;
        double[] integrationOut;
        /** Per-active-lineage clock multiplier, rebuilt once per integrateSegment call. */
        double[] activeClockRates;

        double[] k1;
        double[] k2;
        double[] k3;
        double[] k4;
        double[] y2;
        double[] y3;
        double[] y4;

        double[] sums;
        double[] sumsSquares;
        double[] hValues;
        double[] rValues;
        double[] bValues;
        double[] cValues;
        double[] cSums;
        double[] gradQ;
        /** K*K accumulator; see the doc comment at its use site in rhsVjpInto. */
        double[] migrationGram;

        double[] adjointY0;
        double[] adjointK1;
        double[] adjointK2;
        double[] adjointK3;
        double[] adjointK4;
        double[] vjpY;

        double[] reverseCursorA;
        double[] reverseCursorB;

        /** Operation-level ping-pong buffers for reverse(); see the comment there. */
        double[] reverseOperationA;
        double[] reverseOperationB;

        double[] coalP1;
        double[] coalP2;
        double[] coalParent;
        boolean[] sampleInputs;
        double[] coalescentTimes;

        /**
         * Discardable reverse-pass output buffers, used only when a caller of
         * {@link #evaluateInto} wants ancestral-state scores but not the
         * parameter/clock gradient: the reverse traversal always computes both
         * together (shared VJP), so something must receive the unwanted one.
         * Reused/grown across calls like every other workspace array.
         */
        double[] gradientScratch;
        double[] clockGradientScratch;

        final OperationTapeStore operationTapes = new OperationTapeStore();
    }

    /**
     * Per-epoch cached migration/population transforms, refreshed once per
     * {@link #evaluate}. Fields are package-private for the same reason as
     * {@link Workspace}'s; the class itself is package-private so {@link
     * MascotReverseModeHelper} can read it directly too.
     */
    static final class EpochRates {
        final double[] migrationMatrix;
        final double[] migrationRates;
        final double[] inversePopulation;

        private EpochRates(int stateCount) {
            this.migrationMatrix = new double[stateCount * stateCount];
            this.migrationRates = new double[stateCount * (stateCount - 1)];
            this.inversePopulation = new double[stateCount];
        }
    }

    interface OperationTape {
    }

    /**
     * Reusable operation sequence for one {@link MascotCore} instance. The
     * operation order is stable during fixed-tree HMC, so each evaluation can
     * overwrite the previous tape objects and their grown backing arrays.
     */
    static final class OperationTapeStore {
        private OperationTape[] operations = new OperationTape[16];
        private int operationCount;

        private void reset(int expectedOperationCount) {
            ensureOperationCapacity(expectedOperationCount);
            operationCount = 0;
        }

        // Package-private (not private): read by MascotReverseModeHelper.reverse(),
        // which replays this store back-to-front from a separate top-level class.
        int size() {
            return operationCount;
        }

        OperationTape get(int index) {
            return operations[index];
        }

        // addInterval/addSample/addCoalescent are package-private (not private):
        // called by MascotForwardModeHelper.integrateSegment()/applySampleEvent()/
        // applyCoalescentEvent() from a separate top-level class.
        IntervalTape addInterval(int steps, int activeCount, int stateDimension, int epoch, double h,
                                 int[] activeIds, double[] activeClockRates) {
            int index = nextIndex();
            OperationTape operation = operations[index];
            IntervalTape tape;
            if (operation instanceof IntervalTape) {
                tape = (IntervalTape) operation;
            } else {
                tape = new IntervalTape();
                operations[index] = tape;
            }
            tape.reset(steps, activeCount, stateDimension, epoch, h, activeIds, activeClockRates);
            return tape;
        }

        void addSample(int sampleIndexAfter) {
            int index = nextIndex();
            OperationTape operation = operations[index];
            SampleTape tape;
            if (operation instanceof SampleTape) {
                tape = (SampleTape) operation;
            } else {
                tape = new SampleTape();
                operations[index] = tape;
            }
            tape.reset(sampleIndexAfter);
        }

        void addCoalescent(int epoch, int child1Index, int child2Index, int parentIndexAfter,
                          int movedFromIndexBefore, int movedToIndexAfter, int parentLineageId,
                          double[] p1, double[] p2,
                          double[] parentProbabilities, double lambda, int stateCount) {
            int index = nextIndex();
            OperationTape operation = operations[index];
            CoalescentTape tape;
            if (operation instanceof CoalescentTape) {
                tape = (CoalescentTape) operation;
            } else {
                tape = new CoalescentTape();
                operations[index] = tape;
            }
            tape.reset(epoch, child1Index, child2Index, parentIndexAfter,
                    movedFromIndexBefore, movedToIndexAfter, parentLineageId,
                    p1, p2, parentProbabilities, lambda, stateCount);
        }

        private int nextIndex() {
            ensureOperationCapacity(operationCount + 1);
            return operationCount++;
        }

        private void ensureOperationCapacity(int capacity) {
            if (operations.length < capacity) {
                int newCapacity = operations.length;
                while (newCapacity < capacity) {
                    newCapacity *= 2;
                }
                operations = Arrays.copyOf(operations, newCapacity);
            }
        }
    }

    /**
     * Flat, per-interval RK4 tape: one contiguous {@code steps * stateDimension}
     * array per RK4 stage rather than one small object per step. Fields are
     * package-private (not {@code private}) for the same reason as {@link
     * Workspace}'s: both MascotCore's forward taping methods and {@link
     * MascotReverseModeHelper}'s reverse replay read/write them directly.
     */
    static final class IntervalTape implements OperationTape {
        int steps;
        int activeCount;
        int stateDimension;
        int epoch;
        double h;
        double[] y0;
        double[] y2;
        double[] y3;
        double[] y4;
        /**
         * Per-active-index lineage id and clock multiplier, frozen at the moment
         * this interval was taped (the active set and branchRates are both fixed
         * for the segment's whole duration). Both null together when no clock was
         * supplied to evaluate(...). Snapshotted rather than read live from
         * Workspace/branchRates because later segments overwrite that live state
         * before reverse() replays this one.
         */
        int[] activeIds;
        double[] clockRates;

        private void reset(int steps, int activeCount, int stateDimension, int epoch, double h,
                           int[] sourceActiveIds, double[] sourceClockRates) {
            this.steps = steps;
            this.activeCount = activeCount;
            this.stateDimension = stateDimension;
            this.epoch = epoch;
            this.h = h;
            int storageSize = steps * stateDimension;
            y0 = ensure(y0, storageSize);
            y2 = ensure(y2, storageSize);
            y3 = ensure(y3, storageSize);
            y4 = ensure(y4, storageSize);
            if (sourceActiveIds != null) {
                activeIds = ensureInt(activeIds, activeCount);
                System.arraycopy(sourceActiveIds, 0, activeIds, 0, activeCount);
            } else {
                activeIds = null;
            }
            if (sourceClockRates != null) {
                clockRates = ensure(clockRates, activeCount);
                System.arraycopy(sourceClockRates, 0, clockRates, 0, activeCount);
            } else {
                clockRates = null;
            }
        }
    }

    /** Package-private fields/class for the same reason as {@link IntervalTape}'s. */
    static final class SampleTape implements OperationTape {
        int sampleIndexAfter;

        private void reset(int sampleIndexAfter) {
            this.sampleIndexAfter = sampleIndexAfter;
        }
    }

    /** Package-private fields/class for the same reason as {@link IntervalTape}'s. */
    static final class CoalescentTape implements OperationTape {
        int epoch;
        int child1Index;
        int child2Index;
        int parentIndexAfter;
        int movedFromIndexBefore;
        int movedToIndexAfter;
        /**
         * Stable tree node id (the coalescent parent's {@code event.parent}),
         * as opposed to {@link #parentIndexAfter}, which is a transient
         * active-lineage array position reused across events. Ancestral-state
         * reconstruction must index its output by this id, not by
         * {@code parentIndexAfter}.
         */
        int parentLineageId;
        double[] p1;
        double[] p2;
        double[] parentProbabilities;
        @SuppressWarnings("unused")
        private double lambda;

        private void reset(int epoch, int child1Index, int child2Index, int parentIndexAfter,
                           int movedFromIndexBefore, int movedToIndexAfter, int parentLineageId,
                           double[] p1, double[] p2,
                           double[] parentProbabilities, double lambda, int stateCount) {
            this.epoch = epoch;
            this.child1Index = child1Index;
            this.child2Index = child2Index;
            this.parentIndexAfter = parentIndexAfter;
            this.movedFromIndexBefore = movedFromIndexBefore;
            this.movedToIndexAfter = movedToIndexAfter;
            this.parentLineageId = parentLineageId;
            this.lambda = lambda;
            this.p1 = ensure(this.p1, stateCount);
            this.p2 = ensure(this.p2, stateCount);
            this.parentProbabilities = ensure(this.parentProbabilities, stateCount);
            System.arraycopy(p1, 0, this.p1, 0, stateCount);
            System.arraycopy(p2, 0, this.p2, 0, stateCount);
            System.arraycopy(parentProbabilities, 0, this.parentProbabilities, 0, stateCount);
        }
    }
}
