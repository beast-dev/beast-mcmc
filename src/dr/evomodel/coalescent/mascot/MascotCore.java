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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Allocation-light engine for a MASCOT-style marginal structured coalescent
 * likelihood and the reverse-mode derivative of the same RK4 discretization.
 *
 * The flat parameter layout is per epoch:
 *
 * <pre>
 * log m[0,1], log m[0,2], ..., log m[K-1,K-2], log N[0], ..., log N[K-1]
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

    private static final double TIME_TOLERANCE = 1.0e-14;
    private static final double EPS = 1.0e-300;

    private final int stateCount;
    private final double[] boundaries;
    private final double maxStep;
    private final int migrationParametersPerEpoch;
    private final int parametersPerEpoch;
    private final int parameterCount;
    private final int epochCount;

    private final Workspace workspace = new Workspace();
    private final EpochRates[] epochRates;

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

    public double logLikelihood(Event[] events, double[] theta) {
        return evaluate(events, theta, false, false).logLikelihood;
    }

    public double logLikelihood(PreparedEvents events, double[] theta) {
        return evaluate(events, theta, false, false).logLikelihood;
    }

    public Result likelihoodAndGradient(Event[] events, double[] theta) {
        return evaluate(events, theta, true, false);
    }

    public Result likelihoodAndGradient(PreparedEvents events, double[] theta) {
        return evaluate(events, theta, true, false);
    }

    /**
     * Sorts and validates an event array once so that repeated evaluations against
     * the same fixed tree (only parameters changing) do not repeat the sort or the
     * event-array clone. The returned object is immutable and safe to share across
     * many {@link #evaluate(PreparedEvents, double[], boolean, boolean)} calls.
     */
    public static PreparedEvents prepareEvents(Event[] events) {
        if (events == null || events.length == 0) {
            throw new IllegalArgumentException("at least one tree event is required");
        }
        Event[] sorted = events.clone();
        Arrays.sort(sorted, EVENT_COMPARATOR);

        int maxLineageId = 0;
        for (Event event : events) {
            maxLineageId = Math.max(maxLineageId, event.lineage);
            maxLineageId = Math.max(maxLineageId, event.child1);
            maxLineageId = Math.max(maxLineageId, event.child2);
            maxLineageId = Math.max(maxLineageId, event.parent);
        }

        return new PreparedEvents(sorted, maxLineageId);
    }

    public Result evaluate(Event[] events, double[] theta, boolean computeGradient, boolean checkProbabilities) {
        return evaluate(prepareEvents(events), theta, computeGradient, checkProbabilities);
    }

    public Result evaluate(PreparedEvents prepared, double[] theta, boolean computeGradient, boolean checkProbabilities) {
        if (prepared == null) {
            throw new IllegalArgumentException("prepared events must not be null");
        }
        if (theta == null || theta.length != parameterCount) {
            throw new IllegalArgumentException("theta dimension " + (theta == null ? -1 : theta.length) +
                    " does not match expected dimension " + parameterCount);
        }

        for (int epoch = 0; epoch < epochCount; epoch++) {
            updateEpochRates(theta, epoch, epochRates[epoch]);
        }

        Event[] sortedEvents = prepared.sortedEvents;
        ActiveState state = new ActiveState(sortedEvents.length, stateCount, prepared.maxLineageId);
        double currentTime = 0.0;
        List<OperationTape> operations = computeGradient ? new ArrayList<OperationTape>() : null;

        for (Event event : sortedEvents) {
            if (event.time < currentTime - TIME_TOLERANCE) {
                throw new IllegalArgumentException("events must be sorted by nondecreasing time");
            }

            while (event.time > currentTime + TIME_TOLERANCE) {
                if (state.activeCount == 0) {
                    currentTime = event.time;
                    break;
                }
                double segmentEnd = nextBoundaryAfter(currentTime, event.time);
                int epoch = epochAt(currentTime + TIME_TOLERANCE);
                integrateSegment(state, epochRates[epoch], currentTime, segmentEnd, epoch, operations);
                currentTime = segmentEnd;
                if (checkProbabilities) {
                    checkProbabilities(state);
                }
            }

            currentTime = event.time;
            if (event.type == EventType.SAMPLE) {
                applySampleEvent(state, event, operations);
            } else {
                int epoch = epochAt(event.time);
                applyCoalescentEvent(state, event, epoch, epochRates[epoch], operations);
            }

            if (checkProbabilities) {
                checkProbabilities(state);
            }
        }

        double[] gradient = null;
        if (computeGradient) {
            gradient = reverse(operations, state.activeCount);
        }

        return new Result(state.logLikelihood, gradient, state.copyProbabilities(), state.copyActiveIds());
    }

    // ------------------------------------------------------------------
    // Forward integration
    // ------------------------------------------------------------------

    private void integrateSegment(ActiveState state, EpochRates rates, double start, double end, int epoch,
                                  List<OperationTape> operations) {
        if (!(end > start)) {
            throw new IllegalArgumentException("empty integration segment");
        }

        int steps = Math.max(1, (int) Math.ceil((end - start) / maxStep));
        double h = (end - start) / steps;
        int activeCount = state.activeCount;
        int dim = activeCount * stateCount + 1;

        workspace.integrationState = ensure(workspace.integrationState, dim);
        double[] y = workspace.integrationState;
        packStateInto(state, y);

        workspace.integrationOut = ensure(workspace.integrationOut, dim);
        double[] yOut = workspace.integrationOut;

        if (operations == null) {
            for (int i = 0; i < steps; i++) {
                rk4StepInto(y, activeCount, rates, h, yOut, workspace);
                System.arraycopy(yOut, 0, y, 0, dim);
            }
        } else {
            IntervalTape tape = new IntervalTape(steps, activeCount, dim, epoch, h);
            for (int i = 0; i < steps; i++) {
                rk4StepWithTapeInto(y, activeCount, rates, h, yOut, workspace, tape, i);
                System.arraycopy(yOut, 0, y, 0, dim);
            }
            operations.add(tape);
        }

        unpackStateFrom(y, state);
    }

    private void rk4StepInto(double[] y, int activeCount, EpochRates rates, double h, double[] yOut, Workspace w) {
        int dim = activeCount * stateCount + 1;

        w.k1 = ensure(w.k1, dim);
        rhsInto(y, activeCount, rates, w.k1, w);
        w.y2 = ensure(w.y2, dim);
        addScaledInto(y, w.k1, 0.5 * h, w.y2, dim);
        w.k2 = ensure(w.k2, dim);
        rhsInto(w.y2, activeCount, rates, w.k2, w);
        w.y3 = ensure(w.y3, dim);
        addScaledInto(y, w.k2, 0.5 * h, w.y3, dim);
        w.k3 = ensure(w.k3, dim);
        rhsInto(w.y3, activeCount, rates, w.k3, w);
        w.y4 = ensure(w.y4, dim);
        addScaledInto(y, w.k3, h, w.y4, dim);
        w.k4 = ensure(w.k4, dim);
        rhsInto(w.y4, activeCount, rates, w.k4, w);

        for (int i = 0; i < dim; i++) {
            yOut[i] = y[i] + (h / 6.0) * (w.k1[i] + 2.0 * w.k2[i] + 2.0 * w.k3[i] + w.k4[i]);
        }
    }

    private void rk4StepWithTapeInto(double[] y, int activeCount, EpochRates rates, double h, double[] yOut,
                                     Workspace w, IntervalTape tape, int stepIndex) {
        int dim = tape.stateDimension;
        int offset = stepIndex * dim;
        System.arraycopy(y, 0, tape.y0, offset, dim);

        w.k1 = ensure(w.k1, dim);
        rhsInto(y, activeCount, rates, w.k1, w);
        w.y2 = ensure(w.y2, dim);
        addScaledInto(y, w.k1, 0.5 * h, w.y2, dim);
        System.arraycopy(w.y2, 0, tape.y2, offset, dim);

        w.k2 = ensure(w.k2, dim);
        rhsInto(w.y2, activeCount, rates, w.k2, w);
        w.y3 = ensure(w.y3, dim);
        addScaledInto(y, w.k2, 0.5 * h, w.y3, dim);
        System.arraycopy(w.y3, 0, tape.y3, offset, dim);

        w.k3 = ensure(w.k3, dim);
        rhsInto(w.y3, activeCount, rates, w.k3, w);
        w.y4 = ensure(w.y4, dim);
        addScaledInto(y, w.k3, h, w.y4, dim);
        System.arraycopy(w.y4, 0, tape.y4, offset, dim);

        w.k4 = ensure(w.k4, dim);
        rhsInto(w.y4, activeCount, rates, w.k4, w);

        for (int i = 0; i < dim; i++) {
            yOut[i] = y[i] + (h / 6.0) * (w.k1[i] + 2.0 * w.k2[i] + 2.0 * w.k3[i] + w.k4[i]);
        }
    }

    private void rhsInto(double[] y, int activeCount, EpochRates rates, double[] out, Workspace w) {
        int K = stateCount;
        int stateSize = activeCount * K;

        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            double p0 = y[offset];
            for (int sink = 0; sink < K; sink++) {
                out[offset + sink] = p0 * rates.migrationMatrix[sink];
            }
            for (int source = 1; source < K; source++) {
                double p = y[offset + source];
                int row = source * K;
                for (int sink = 0; sink < K; sink++) {
                    out[offset + sink] += p * rates.migrationMatrix[row + sink];
                }
            }
        }

        w.sums = ensure(w.sums, K);
        w.sumsSquares = ensure(w.sumsSquares, K);
        Arrays.fill(w.sums, 0, K, 0.0);
        Arrays.fill(w.sumsSquares, 0, K, 0.0);
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            for (int state = 0; state < K; state++) {
                double p = y[offset + state];
                w.sums[state] += p;
                w.sumsSquares[state] += p * p;
            }
        }

        double hazard = 0.0;
        for (int state = 0; state < K; state++) {
            hazard += 0.5 * rates.inversePopulation[state] * (w.sums[state] * w.sums[state] - w.sumsSquares[state]);
        }

        w.hValues = ensure(w.hValues, stateSize);
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            double r = 0.0;
            for (int state = 0; state < K; state++) {
                double h = (w.sums[state] - y[offset + state]) * rates.inversePopulation[state];
                w.hValues[offset + state] = h;
                r += y[offset + state] * h;
            }
            for (int state = 0; state < K; state++) {
                out[offset + state] += -y[offset + state] * (w.hValues[offset + state] - r);
            }
        }

        out[stateSize] = -hazard;
    }

    // ------------------------------------------------------------------
    // Reverse pass
    // ------------------------------------------------------------------

    private double[] reverse(List<OperationTape> operations, int finalActiveCount) {
        double[] gradient = new double[parameterCount];
        int dim = finalActiveCount * stateCount + 1;
        double[] adjointAfter = new double[dim];
        adjointAfter[dim - 1] = 1.0;

        for (int opIndex = operations.size() - 1; opIndex >= 0; opIndex--) {
            OperationTape operation = operations.get(opIndex);
            if (operation instanceof IntervalTape) {
                adjointAfter = reverseInterval((IntervalTape) operation, adjointAfter, gradient);
            } else if (operation instanceof CoalescentTape) {
                adjointAfter = reverseCoalescent((CoalescentTape) operation, adjointAfter, gradient);
            } else if (operation instanceof SampleTape) {
                adjointAfter = reverseSample((SampleTape) operation, adjointAfter);
            } else {
                throw new IllegalArgumentException("unknown tape operation: " + operation.getClass());
            }
        }

        return gradient;
    }

    private double[] reverseInterval(IntervalTape tape, double[] adjointAfter, double[] gradient) {
        int dim = tape.stateDimension;
        EpochRates rates = epochRates[tape.epoch];

        workspace.reverseCursorA = ensure(workspace.reverseCursorA, dim);
        workspace.reverseCursorB = ensure(workspace.reverseCursorB, dim);
        System.arraycopy(adjointAfter, 0, workspace.reverseCursorA, 0, dim);
        double[] cursor = workspace.reverseCursorA;
        double[] next = workspace.reverseCursorB;

        for (int step = tape.steps - 1; step >= 0; step--) {
            int offset = step * dim;
            reverseStepInto(tape, offset, rates, cursor, next, gradient, workspace);
            double[] swap = cursor;
            cursor = next;
            next = swap;
        }

        double[] result = new double[dim];
        System.arraycopy(cursor, 0, result, 0, dim);
        return result;
    }

    private void reverseStepInto(IntervalTape tape, int offset, EpochRates rates, double[] adjointAfter,
                                 double[] adjointBeforeOut, double[] gradient, Workspace w) {
        int dim = tape.stateDimension;
        int activeCount = tape.activeCount;
        int epoch = tape.epoch;
        double h = tape.h;

        w.adjointY0 = ensure(w.adjointY0, dim);
        System.arraycopy(adjointAfter, 0, w.adjointY0, 0, dim);

        w.adjointK1 = ensure(w.adjointK1, dim);
        scaleInto(adjointAfter, h / 6.0, w.adjointK1, dim);
        w.adjointK2 = ensure(w.adjointK2, dim);
        scaleInto(adjointAfter, h / 3.0, w.adjointK2, dim);
        w.adjointK3 = ensure(w.adjointK3, dim);
        scaleInto(adjointAfter, h / 3.0, w.adjointK3, dim);
        w.adjointK4 = ensure(w.adjointK4, dim);
        scaleInto(adjointAfter, h / 6.0, w.adjointK4, dim);

        w.tapeSlice = ensure(w.tapeSlice, dim);
        w.vjpY = ensure(w.vjpY, dim);

        System.arraycopy(tape.y4, offset, w.tapeSlice, 0, dim);
        rhsVjpInto(w.tapeSlice, activeCount, rates, epoch, w.adjointK4, w.vjpY, gradient, w);
        addInPlace(w.adjointY0, w.vjpY, dim);
        addScaledInPlace(w.adjointK3, w.vjpY, h, dim);

        System.arraycopy(tape.y3, offset, w.tapeSlice, 0, dim);
        rhsVjpInto(w.tapeSlice, activeCount, rates, epoch, w.adjointK3, w.vjpY, gradient, w);
        addInPlace(w.adjointY0, w.vjpY, dim);
        addScaledInPlace(w.adjointK2, w.vjpY, 0.5 * h, dim);

        System.arraycopy(tape.y2, offset, w.tapeSlice, 0, dim);
        rhsVjpInto(w.tapeSlice, activeCount, rates, epoch, w.adjointK2, w.vjpY, gradient, w);
        addInPlace(w.adjointY0, w.vjpY, dim);
        addScaledInPlace(w.adjointK1, w.vjpY, 0.5 * h, dim);

        System.arraycopy(tape.y0, offset, w.tapeSlice, 0, dim);
        rhsVjpInto(w.tapeSlice, activeCount, rates, epoch, w.adjointK1, w.vjpY, gradient, w);
        addInPlace(w.adjointY0, w.vjpY, dim);

        System.arraycopy(w.adjointY0, 0, adjointBeforeOut, 0, dim);
    }

    private void rhsVjpInto(double[] y, int activeCount, EpochRates rates, int epoch, double[] adjointRhs,
                            double[] adjointYOut, double[] gradient, Workspace w) {
        int K = stateCount;
        int stateSize = activeCount * K;
        // adjointYOut[stateSize] (the VJP component for the input's log-likelihood
        // slot) is always zero: the forward RHS's out[stateSize] = -hazard never
        // reads y[stateSize], so nothing below ever writes this index again.
        adjointYOut[stateSize] = 0.0;

        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            for (int source = 0; source < K; source++) {
                double v = 0.0;
                int row = source * K;
                for (int sink = 0; sink < K; sink++) {
                    v += adjointRhs[offset + sink] * rates.migrationMatrix[row + sink];
                }
                adjointYOut[offset + source] = v;
            }
        }

        int thetaOffset = epoch * parametersPerEpoch;
        int rateIndex = 0;
        for (int source = 0; source < K; source++) {
            for (int sink = 0; sink < K; sink++) {
                if (source == sink) {
                    continue;
                }
                double rate = rates.migrationRates[rateIndex];
                double contribution = 0.0;
                for (int lineage = 0; lineage < activeCount; lineage++) {
                    int offset = lineage * K;
                    contribution += y[offset + source] * (adjointRhs[offset + sink] - adjointRhs[offset + source]);
                }
                gradient[thetaOffset + rateIndex] += rate * contribution;
                rateIndex++;
            }
        }

        w.sums = ensure(w.sums, K);
        w.sumsSquares = ensure(w.sumsSquares, K);
        Arrays.fill(w.sums, 0, K, 0.0);
        Arrays.fill(w.sumsSquares, 0, K, 0.0);
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            for (int state = 0; state < K; state++) {
                double p = y[offset + state];
                w.sums[state] += p;
                w.sumsSquares[state] += p * p;
            }
        }

        w.hValues = ensure(w.hValues, stateSize);
        w.rValues = ensure(w.rValues, activeCount);
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            double r = 0.0;
            for (int state = 0; state < K; state++) {
                double h = (w.sums[state] - y[offset + state]) * rates.inversePopulation[state];
                w.hValues[offset + state] = h;
                r += y[offset + state] * h;
            }
            w.rValues[lineage] = r;
        }

        w.bValues = ensure(w.bValues, activeCount);
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            double b = 0.0;
            for (int state = 0; state < K; state++) {
                b += adjointRhs[offset + state] * y[offset + state];
            }
            w.bValues[lineage] = b;
        }

        w.cSums = ensure(w.cSums, K);
        w.cValues = ensure(w.cValues, stateSize);
        Arrays.fill(w.cSums, 0, K, 0.0);
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            double b = w.bValues[lineage];
            for (int state = 0; state < K; state++) {
                double upstream = adjointRhs[offset + state];
                adjointYOut[offset + state] += upstream * (w.rValues[lineage] - w.hValues[offset + state]) +
                        b * w.hValues[offset + state];
                double c = y[offset + state] * (b - upstream);
                w.cValues[offset + state] = c;
                w.cSums[state] += c;
            }
        }

        w.gradQ = ensure(w.gradQ, K);
        Arrays.fill(w.gradQ, 0, K, 0.0);
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            for (int state = 0; state < K; state++) {
                double c = w.cValues[offset + state];
                adjointYOut[offset + state] += rates.inversePopulation[state] * (w.cSums[state] - c);
                w.gradQ[state] += c * (w.sums[state] - y[offset + state]);
            }
        }

        double ellAdjoint = adjointRhs[stateSize];
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            for (int state = 0; state < K; state++) {
                adjointYOut[offset + state] += -ellAdjoint * w.hValues[offset + state];
            }
        }
        for (int state = 0; state < K; state++) {
            double pairSum = 0.5 * (w.sums[state] * w.sums[state] - w.sumsSquares[state]);
            w.gradQ[state] += -ellAdjoint * pairSum;
        }

        int etaOffset = thetaOffset + migrationParametersPerEpoch;
        for (int state = 0; state < K; state++) {
            gradient[etaOffset + state] += -rates.inversePopulation[state] * w.gradQ[state];
        }
    }

    // ------------------------------------------------------------------
    // Events
    // ------------------------------------------------------------------

    private void applySampleEvent(ActiveState state, Event event, List<OperationTape> operations) {
        if (event.lineage < 0) {
            throw new IllegalArgumentException("sample event has no lineage id");
        }
        if (state.isActive(event.lineage)) {
            throw new IllegalArgumentException("lineage " + event.lineage + " is already active");
        }

        double[] probabilities = normalizedSampleProbabilities(event);
        int index = state.activeCount;
        state.ensureCapacity(index + 1);
        state.activeIds[index] = event.lineage;
        System.arraycopy(probabilities, 0, state.probabilities, index * stateCount, stateCount);
        state.setActiveIndex(event.lineage, index);
        state.activeCount++;

        if (operations != null) {
            operations.add(new SampleTape(index));
        }
    }

    private void applyCoalescentEvent(ActiveState state, Event event, int epoch, EpochRates rates,
                                      List<OperationTape> operations) {
        int first = state.activeIndexOf(event.child1);
        int second = state.activeIndexOf(event.child2);
        if (first < 0 || second < 0) {
            throw new IllegalArgumentException("coalescent children are not both active: " +
                    event.child1 + ", " + event.child2);
        }
        if (first == second) {
            throw new IllegalArgumentException("coalescent children must be distinct");
        }
        if (state.isActive(event.parent)) {
            throw new IllegalArgumentException("coalescent parent is already active: " + event.parent);
        }

        boolean recordTape = operations != null;
        double[] q = rates.inversePopulation;

        workspace.coalP1 = ensure(workspace.coalP1, stateCount);
        workspace.coalP2 = ensure(workspace.coalP2, stateCount);
        double[] p1 = workspace.coalP1;
        double[] p2 = workspace.coalP2;
        System.arraycopy(state.probabilities, first * stateCount, p1, 0, stateCount);
        System.arraycopy(state.probabilities, second * stateCount, p2, 0, stateCount);

        double[] parentProbabilities;
        if (recordTape) {
            parentProbabilities = new double[stateCount];
        } else {
            workspace.coalParent = ensure(workspace.coalParent, stateCount);
            parentProbabilities = workspace.coalParent;
        }

        double lambda = 0.0;
        for (int s = 0; s < stateCount; s++) {
            parentProbabilities[s] = p1[s] * p2[s] * q[s];
            lambda += parentProbabilities[s];
        }
        if (!(lambda > 0.0) || !Double.isFinite(lambda)) {
            throw new IllegalArgumentException("invalid coalescent rate: " + lambda);
        }
        for (int s = 0; s < stateCount; s++) {
            parentProbabilities[s] /= lambda;
        }

        int beforeCount = state.activeCount;
        int afterCount = beforeCount - 1;
        int[] keepBeforeIndices = recordTape ? new int[afterCount - 1] : null;

        int out = 0;
        for (int i = 0; i < beforeCount; i++) {
            if (i == first || i == second) {
                continue;
            }
            if (out != i) {
                state.activeIds[out] = state.activeIds[i];
                System.arraycopy(state.probabilities, i * stateCount, state.probabilities, out * stateCount, stateCount);
            }
            state.setActiveIndex(state.activeIds[out], out);
            if (recordTape) {
                keepBeforeIndices[out] = i;
            }
            out++;
        }

        int parentIndexAfter = out;
        state.activeIds[parentIndexAfter] = event.parent;
        System.arraycopy(parentProbabilities, 0, state.probabilities, parentIndexAfter * stateCount, stateCount);
        state.setActiveIndex(event.parent, parentIndexAfter);
        state.clearActiveIndex(event.child1);
        state.clearActiveIndex(event.child2);

        state.activeCount = afterCount;
        state.logLikelihood += Math.log(lambda);

        if (operations != null) {
            operations.add(new CoalescentTape(epoch, first, second, parentIndexAfter, keepBeforeIndices,
                    p1.clone(), p2.clone(), parentProbabilities, lambda));
        }
    }

    private double[] normalizedSampleProbabilities(Event event) {
        double[] probabilities = new double[stateCount];
        if (event.stateProbabilities == null) {
            if (event.state < 0 || event.state >= stateCount) {
                throw new IllegalArgumentException("sample state out of range: " + event.state);
            }
            probabilities[event.state] = 1.0;
            return probabilities;
        }
        if (event.stateProbabilities.length != stateCount) {
            throw new IllegalArgumentException("sample probability dimension mismatch");
        }
        double sum = 0.0;
        for (int s = 0; s < stateCount; s++) {
            if (event.stateProbabilities[s] < 0.0) {
                throw new IllegalArgumentException("sample probabilities must be nonnegative");
            }
            sum += event.stateProbabilities[s];
        }
        if (!(sum > 0.0)) {
            throw new IllegalArgumentException("sample probabilities must have positive sum");
        }
        for (int s = 0; s < stateCount; s++) {
            probabilities[s] = event.stateProbabilities[s] / sum;
        }
        return probabilities;
    }

    private double[] reverseSample(SampleTape tape, double[] adjointAfter) {
        int afterCount = (adjointAfter.length - 1) / stateCount;
        int beforeCount = afterCount - 1;
        double[] adjointBefore = new double[beforeCount * stateCount + 1];

        int out = 0;
        for (int i = 0; i < afterCount; i++) {
            if (i == tape.sampleIndexAfter) {
                continue;
            }
            System.arraycopy(adjointAfter, i * stateCount, adjointBefore, out * stateCount, stateCount);
            out++;
        }
        adjointBefore[adjointBefore.length - 1] = adjointAfter[adjointAfter.length - 1];
        return adjointBefore;
    }

    private double[] reverseCoalescent(CoalescentTape tape, double[] adjointAfter, double[] gradient) {
        int afterCount = (adjointAfter.length - 1) / stateCount;
        int beforeCount = afterCount + 1;
        double[] adjointBefore = new double[beforeCount * stateCount + 1];

        for (int afterIndex = 0; afterIndex < tape.keepBeforeIndices.length; afterIndex++) {
            int beforeIndex = tape.keepBeforeIndices[afterIndex];
            System.arraycopy(adjointAfter, afterIndex * stateCount,
                    adjointBefore, beforeIndex * stateCount, stateCount);
        }

        int parentOffset = tape.parentIndexAfter * stateCount;
        double dot = 0.0;
        for (int s = 0; s < stateCount; s++) {
            dot += adjointAfter[parentOffset + s] * tape.parentProbabilities[s];
        }

        int thetaOffset = tape.epoch * parametersPerEpoch + migrationParametersPerEpoch;
        double ellAdjoint = adjointAfter[adjointAfter.length - 1];
        for (int s = 0; s < stateCount; s++) {
            double p1 = Math.max(tape.p1[s], EPS);
            double p2 = Math.max(tape.p2[s], EPS);
            double centered = adjointAfter[parentOffset + s] - dot;
            double bar = tape.parentProbabilities[s] * centered;

            adjointBefore[tape.child1Index * stateCount + s] += bar / p1;
            adjointBefore[tape.child2Index * stateCount + s] += bar / p2;
            gradient[thetaOffset + s] += -bar;

            adjointBefore[tape.child1Index * stateCount + s] += ellAdjoint * tape.parentProbabilities[s] / p1;
            adjointBefore[tape.child2Index * stateCount + s] += ellAdjoint * tape.parentProbabilities[s] / p2;
            gradient[thetaOffset + s] += -ellAdjoint * tape.parentProbabilities[s];
        }

        adjointBefore[adjointBefore.length - 1] = ellAdjoint;
        return adjointBefore;
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
                double rate = Math.exp(theta[thetaOffset + index]);
                rates.migrationRates[index] = rate;
                rates.migrationMatrix[row + sink] = rate;
                rowSum += rate;
                index++;
            }
            rates.migrationMatrix[row + source] = -rowSum;
        }

        int etaOffset = thetaOffset + migrationParametersPerEpoch;
        for (int state = 0; state < stateCount; state++) {
            rates.inversePopulation[state] = Math.exp(-theta[etaOffset + state]);
        }
    }

    // ------------------------------------------------------------------
    // Epoch/time bookkeeping
    // ------------------------------------------------------------------

    private int epochAt(double t) {
        if (t < -TIME_TOLERANCE) {
            throw new IllegalArgumentException("time is before zero: " + t);
        }
        int epoch = 0;
        while (epoch + 1 < boundaries.length && t >= boundaries[epoch + 1] - TIME_TOLERANCE) {
            epoch++;
        }
        if (epoch >= epochCount) {
            return epochCount - 1;
        }
        return epoch;
    }

    private double nextBoundaryAfter(double start, double stop) {
        for (int i = 1; i < boundaries.length; i++) {
            double boundary = boundaries[i];
            if (boundary > start + TIME_TOLERANCE) {
                if (boundary <= stop + TIME_TOLERANCE) {
                    return Math.min(boundary, stop);
                }
                return stop;
            }
        }
        return stop;
    }

    // ------------------------------------------------------------------
    // Pack/unpack and small array helpers
    // ------------------------------------------------------------------

    private void packStateInto(ActiveState state, double[] y) {
        int stateSize = state.activeCount * stateCount;
        System.arraycopy(state.probabilities, 0, y, 0, stateSize);
        y[stateSize] = state.logLikelihood;
    }

    private void unpackStateFrom(double[] y, ActiveState state) {
        int stateSize = state.activeCount * stateCount;
        state.ensureCapacity(state.activeCount);
        System.arraycopy(y, 0, state.probabilities, 0, stateSize);
        state.logLikelihood = y[stateSize];
    }

    private static double[] ensure(double[] array, int size) {
        return (array == null || array.length < size) ? new double[size] : array;
    }

    private static void addScaledInto(double[] x, double[] dx, double scale, double[] out, int n) {
        for (int i = 0; i < n; i++) {
            out[i] = x[i] + scale * dx[i];
        }
    }

    private static void scaleInto(double[] x, double scale, double[] out, int n) {
        for (int i = 0; i < n; i++) {
            out[i] = scale * x[i];
        }
    }

    private static void addInPlace(double[] destination, double[] source, int n) {
        for (int i = 0; i < n; i++) {
            destination[i] += source[i];
        }
    }

    private static void addScaledInPlace(double[] destination, double[] source, double scale, int n) {
        for (int i = 0; i < n; i++) {
            destination[i] += scale * source[i];
        }
    }

    private void checkProbabilities(ActiveState state) {
        for (int lineage = 0; lineage < state.activeCount; lineage++) {
            double sum = 0.0;
            int offset = lineage * stateCount;
            for (int s = 0; s < stateCount; s++) {
                double p = state.probabilities[offset + s];
                if (p < -1.0e-8 || !Double.isFinite(p)) {
                    throw new IllegalStateException("invalid probability for lineage " +
                            state.activeIds[lineage] + ": " + p);
                }
                sum += p;
            }
            if (Math.abs(sum - 1.0) > 1.0e-6) {
                throw new IllegalStateException("lineage probabilities do not sum to one for lineage " +
                        state.activeIds[lineage] + ": " + sum);
            }
        }
    }

    private static final Comparator<Event> EVENT_COMPARATOR = new Comparator<Event>() {
        @Override
        public int compare(Event first, Event second) {
            int time = Double.compare(first.time, second.time);
            if (time != 0) {
                return time;
            }
            int type = first.type.compareRank - second.type.compareRank;
            if (type != 0) {
                return type;
            }
            return first.sortId() - second.sortId();
        }
    };

    public enum EventType {
        SAMPLE(0),
        COALESCENCE(1);

        private final int compareRank;

        EventType(int compareRank) {
            this.compareRank = compareRank;
        }
    }

    public static final class Event {
        private final double time;
        private final EventType type;
        private final int lineage;
        private final int state;
        private final double[] stateProbabilities;
        private final int child1;
        private final int child2;
        private final int parent;

        private Event(double time, EventType type, int lineage, int state, double[] stateProbabilities,
                      int child1, int child2, int parent) {
            if (time < -TIME_TOLERANCE || !Double.isFinite(time)) {
                throw new IllegalArgumentException("invalid event time: " + time);
            }
            this.time = time;
            this.type = type;
            this.lineage = lineage;
            this.state = state;
            this.stateProbabilities = stateProbabilities == null ? null : stateProbabilities.clone();
            this.child1 = child1;
            this.child2 = child2;
            this.parent = parent;
        }

        public static Event sample(double time, int lineage, int state) {
            return new Event(time, EventType.SAMPLE, lineage, state, null, -1, -1, -1);
        }

        public static Event sample(double time, int lineage, double[] stateProbabilities) {
            return new Event(time, EventType.SAMPLE, lineage, -1, stateProbabilities, -1, -1, -1);
        }

        public static Event coalescence(double time, int child1, int child2, int parent) {
            return new Event(time, EventType.COALESCENCE, -1, -1, null, child1, child2, parent);
        }

        public double getTime() {
            return time;
        }

        public EventType getType() {
            return type;
        }

        public int getLineage() {
            return lineage;
        }

        public int getState() {
            return state;
        }

        public int getChild1() {
            return child1;
        }

        public int getChild2() {
            return child2;
        }

        public int getParent() {
            return parent;
        }

        private int sortId() {
            if (type == EventType.SAMPLE) {
                return lineage;
            }
            return parent;
        }
    }

    /**
     * An already-sorted, validated event sequence. Building one requires a full
     * clone and sort of the input array; evaluating against a {@code PreparedEvents}
     * does not. Callers that repeatedly evaluate the same fixed tree under changing
     * parameters (the common MCMC/HMC pattern) should build this once per tree
     * change and reuse it across evaluations.
     */
    public static final class PreparedEvents {
        private final Event[] sortedEvents;
        private final int maxLineageId;

        private PreparedEvents(Event[] sortedEvents, int maxLineageId) {
            this.sortedEvents = sortedEvents;
            this.maxLineageId = maxLineageId;
        }
    }

    public static final class Result {
        public final double logLikelihood;
        public final double[] gradient;
        public final double[] rootProbabilities;
        public final int[] activeLineages;

        private Result(double logLikelihood, double[] gradient, double[] rootProbabilities, int[] activeLineages) {
            this.logLikelihood = logLikelihood;
            this.gradient = gradient;
            this.rootProbabilities = rootProbabilities;
            this.activeLineages = activeLineages;
        }
    }

    private static final class ActiveState {
        private int[] activeIds;
        private double[] probabilities;
        private int activeCount;
        private double logLikelihood;
        private final int stateCount;
        private int[] lineageToActiveIndex;

        private ActiveState(int capacity, int stateCount, int maxLineageId) {
            this.stateCount = stateCount;
            this.activeIds = new int[Math.max(1, capacity)];
            this.probabilities = new double[Math.max(1, capacity * stateCount)];
            this.activeCount = 0;
            this.logLikelihood = 0.0;
            this.lineageToActiveIndex = new int[Math.max(1, maxLineageId + 1)];
            Arrays.fill(this.lineageToActiveIndex, -1);
        }

        private void ensureCapacity(int capacity) {
            if (activeIds.length < capacity) {
                activeIds = Arrays.copyOf(activeIds, capacity);
            }
            int probabilityCapacity = capacity * stateCount;
            if (probabilities.length < probabilityCapacity) {
                probabilities = Arrays.copyOf(probabilities, probabilityCapacity);
            }
        }

        private boolean isActive(int lineageId) {
            return lineageId >= 0 && lineageId < lineageToActiveIndex.length && lineageToActiveIndex[lineageId] >= 0;
        }

        private int activeIndexOf(int lineageId) {
            if (lineageId < 0 || lineageId >= lineageToActiveIndex.length) {
                return -1;
            }
            return lineageToActiveIndex[lineageId];
        }

        private void setActiveIndex(int lineageId, int index) {
            lineageToActiveIndex[lineageId] = index;
        }

        private void clearActiveIndex(int lineageId) {
            lineageToActiveIndex[lineageId] = -1;
        }

        private double[] copyProbabilities() {
            return Arrays.copyOf(probabilities, activeCount * stateCount);
        }

        private int[] copyActiveIds() {
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
     * measurable self time (see MascotCoreProfileDriver).
     */
    private static final class Workspace {
        double[] integrationState;
        double[] integrationOut;

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

        double[] adjointY0;
        double[] adjointK1;
        double[] adjointK2;
        double[] adjointK3;
        double[] adjointK4;
        double[] tapeSlice;
        double[] vjpY;

        double[] reverseCursorA;
        double[] reverseCursorB;

        double[] coalP1;
        double[] coalP2;
        double[] coalParent;
    }

    /**
     * Per-epoch cached migration/population transforms, refreshed once per
     * {@link #evaluate}. Fields are package-private for the same reason as
     * {@link Workspace}'s.
     */
    private static final class EpochRates {
        final double[] migrationMatrix;
        final double[] migrationRates;
        final double[] inversePopulation;

        private EpochRates(int stateCount) {
            this.migrationMatrix = new double[stateCount * stateCount];
            this.migrationRates = new double[stateCount * (stateCount - 1)];
            this.inversePopulation = new double[stateCount];
        }
    }

    private interface OperationTape {
    }

    /**
     * Flat, per-interval RK4 tape: one contiguous {@code steps * stateDimension}
     * array per RK4 stage rather than one small object per step.
     */
    private static final class IntervalTape implements OperationTape {
        private final int steps;
        private final int activeCount;
        private final int stateDimension;
        private final int epoch;
        private final double h;
        private final double[] y0;
        private final double[] y2;
        private final double[] y3;
        private final double[] y4;

        private IntervalTape(int steps, int activeCount, int stateDimension, int epoch, double h) {
            this.steps = steps;
            this.activeCount = activeCount;
            this.stateDimension = stateDimension;
            this.epoch = epoch;
            this.h = h;
            this.y0 = new double[steps * stateDimension];
            this.y2 = new double[steps * stateDimension];
            this.y3 = new double[steps * stateDimension];
            this.y4 = new double[steps * stateDimension];
        }
    }

    private static final class SampleTape implements OperationTape {
        private final int sampleIndexAfter;

        private SampleTape(int sampleIndexAfter) {
            this.sampleIndexAfter = sampleIndexAfter;
        }
    }

    private static final class CoalescentTape implements OperationTape {
        private final int epoch;
        private final int child1Index;
        private final int child2Index;
        private final int parentIndexAfter;
        private final int[] keepBeforeIndices;
        private final double[] p1;
        private final double[] p2;
        private final double[] parentProbabilities;
        @SuppressWarnings("unused")
        private final double lambda;

        private CoalescentTape(int epoch, int child1Index, int child2Index, int parentIndexAfter,
                               int[] keepBeforeIndices, double[] p1, double[] p2,
                               double[] parentProbabilities, double lambda) {
            this.epoch = epoch;
            this.child1Index = child1Index;
            this.child2Index = child2Index;
            this.parentIndexAfter = parentIndexAfter;
            this.keepBeforeIndices = keepBeforeIndices;
            this.p1 = p1;
            this.p2 = p2;
            this.parentProbabilities = parentProbabilities;
            this.lambda = lambda;
        }
    }
}
