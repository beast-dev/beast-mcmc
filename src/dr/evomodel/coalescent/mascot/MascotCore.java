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

import java.util.Arrays;

/**
 * Allocation-light MASCOT core with reverse-mode derivatives for the same RK4
 * discretization. Per-epoch parameter order is row-major off-diagonal migration
 * rates followed by log population sizes. Instances are not thread-safe.
 */
public final class MascotCore implements MascotLikelihoodBackend {

    // Also used by MascotForwardModeHelper's event-time and epoch bookkeeping.
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
    private final MascotForwardModeHelper forwardHelper;
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

    public double logLikelihood(MascotPreparedInput prepared, double[] theta) {
        return evaluate(prepared, theta, null, false, false, false, false).logLikelihood;
    }

    public Result likelihoodAndGradient(MascotPreparedInput prepared, double[] theta) {
        return evaluate(prepared, theta, null, true, false, false, false);
    }

    public Result likelihoodAndAncestralStates(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                                               boolean checkProbabilities) {
        return evaluate(prepared, theta, branchRates, false, true, checkProbabilities, false);
    }

    /**
     * Allocating wrapper for tests and callers that want a {@link Result}.
     * {@code branchRates} index biological lineage ids and scale only migration.
     * {@code copyFinalState} controls root-probability and active-lineage copies.
     */
    public Result evaluate(MascotPreparedInput prepared, double[] theta, double[] branchRates, boolean computeGradient,
                           boolean computeAncestralStates, boolean checkProbabilities, boolean copyFinalState) {
        return evaluate(prepared, theta, branchRates, computeGradient, computeAncestralStates, checkProbabilities,
                copyFinalState, null);
    }

    /**
     * Test-only finite-difference hook: optional node/state log weights are applied
     * before coalescent-parent normalization. Package-private by design.
     */
    Result evaluateWithNodeLogWeightsForTesting(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                                                double[] nodeLogWeights, boolean checkProbabilities) {
        return evaluate(prepared, theta, branchRates, false, false, checkProbabilities, false, nodeLogWeights);
    }

    private Result evaluate(MascotPreparedInput prepared, double[] theta, double[] branchRates, boolean computeGradient,
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
     * Caller-owned-output production path. Null outputs are skipped except that
     * gradient and ancestral-state requests share one reverse pass.
     */
    @Override
    public double evaluateInto(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                               double[] gradientOut, double[] clockGradientOut, double[] ancestralStateScoresOut,
                               boolean checkProbabilities) {
        return evaluateCore(prepared, theta, branchRates, checkProbabilities, null,
                gradientOut, clockGradientOut, ancestralStateScoresOut);
    }

    /** Scalar-only evaluation: never builds the reverse tape, never allocates a {@link Result}. */
    @Override
    public double evaluateLikelihood(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                                     boolean checkProbabilities) {
        return evaluateCore(prepared, theta, branchRates, checkProbabilities, null, null, null, null);
    }

    private double evaluateCore(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                                boolean checkProbabilities, double[] nodeLogWeights,
                                double[] gradientOut, double[] clockGradientOut, double[] ancestralStateScoresOut) {
        if (prepared == null) {
            throw new IllegalArgumentException("prepared input must not be null");
        }
        if (prepared.tipData.stateCount != stateCount) {
            throw new IllegalArgumentException("tip data stateCount " + prepared.tipData.stateCount +
                    " does not match MascotCore stateCount " + stateCount);
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
            // Reuse tape objects across fixed-tree HMC evaluations.
            operations = workspace.operationTapes;
            operations.reset(prepared.schedule.getIntervalCount() * 2 + epochCount + 1);
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
                // Zero rates are valid (e.g. inactive BSSVS indicators);
                // only negative or non-finite rates break the equations.
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

    // Shared with MascotReverseModeHelper.
    static double[] ensure(double[] array, int size) {
        return (array == null || array.length < size) ? new double[size] : array;
    }

    private static int[] ensureInt(int[] array, int size) {
        return (array == null || array.length < size) ? new int[size] : array;
    }

    // Shared with MascotForwardModeHelper's RK4 steps.
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

    public static final class Result {
        public final double logLikelihood;
        public final double[] gradient;
        /**
         * d(logLikelihood)/d(branchRate[lineageId]); null unless a gradient and
         * branch rates were both requested.
         */
        public final double[] clockGradient;
        public final double[] rootProbabilities;
        public final int[] activeLineages;
        /**
         * Node-major adjoint scores by stable tree node id, not active-lineage slot.
         * Tip/unused rows are NaN. These are RK4 sensitivities, not posterior
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
     * Parameter-proposal numerical failure. Likelihood wrappers may map this to
     * {@code -Infinity}; gradient callers should fail instead.
     */
    public static final class NumericalException extends RuntimeException {
        public NumericalException(String message) {
            super(message);
        }
    }

    /**
     * Reused active-lineage state. Generation stamps make reset O(1);
     * package-private fields are shared with helper hot loops.
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
     * Per-instance growth-only scratch memory. Not thread-safe.
     * Package-private fields are shared with forward/reverse helper hot loops.
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
        /** K*K accumulator for reverse-mode migration terms. */
        double[] migrationGram;

        double[] adjointY0;
        double[] adjointK1;
        double[] adjointK2;
        double[] adjointK3;
        double[] adjointK4;
        double[] vjpY;

        double[] reverseCursorA;
        double[] reverseCursorB;

        /** Operation-level ping-pong buffers for reverse replay. */
        double[] reverseOperationA;
        double[] reverseOperationB;

        double[] coalP1;
        double[] coalP2;
        double[] coalParent;
        double[] coalescentTimes;

        /** Discardable outputs when callers request only one shared reverse product. */
        double[] gradientScratch;
        double[] clockGradientScratch;

        final OperationTapeStore operationTapes = new OperationTapeStore();
    }

    /**
     * Per-epoch migration/population transforms refreshed once per evaluation.
     * Package-private for direct helper access.
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
     * Reusable operation sequence for one core instance.
     * Evaluations overwrite existing tape objects and grown arrays.
     */
    static final class OperationTapeStore {
        private OperationTape[] operations = new OperationTape[16];
        private int operationCount;

        private void reset(int expectedOperationCount) {
            ensureOperationCapacity(expectedOperationCount);
            operationCount = 0;
        }

        int size() {
            return operationCount;
        }

        OperationTape get(int index) {
            return operations[index];
        }

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
     * Flat per-interval RK4 tape: one contiguous array per stage.
     * Package-private fields are shared with reverse replay.
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
        /** Frozen active lineage ids and clock multipliers for reverse replay. */
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

    static final class SampleTape implements OperationTape {
        int sampleIndexAfter;

        private void reset(int sampleIndexAfter) {
            this.sampleIndexAfter = sampleIndexAfter;
        }
    }

    static final class CoalescentTape implements OperationTape {
        int epoch;
        int child1Index;
        int child2Index;
        int parentIndexAfter;
        int movedFromIndexBefore;
        int movedToIndexAfter;
        /** Stable tree node id for ancestral-state output indexing. */
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
