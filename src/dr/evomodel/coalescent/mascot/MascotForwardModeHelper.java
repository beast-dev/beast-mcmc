/*
 * MascotForwardModeHelper.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation;

import java.util.Arrays;

/**
 * The forward RK4 integration and event application for one {@link
 * MascotCore} instance: walks a {@link MascotCore.PreparedEvents} sequence
 * forward in time, advancing a {@link MascotCore.ActiveState} through
 * integration segments, samples, and coalescences. When a non-null {@link
 * MascotCore.OperationTapeStore} is supplied, this class also records the
 * tape {@link MascotReverseModeHelper} later replays backward -- taping is
 * necessarily a forward-pass concern, since the intermediate values it
 * captures cannot be reconstructed after the fact, so this class (not {@link
 * MascotReverseModeHelper}) owns both the integration and the recording.
 * {@link MascotCore} itself owns only {@code evaluate(...)}/{@code
 * evaluateInto(...)} orchestration: per-call validation, the shared {@link
 * MascotCore.EpochRates} setup ({@code updateEpochRates}, read identically by
 * both this class and {@link MascotReverseModeHelper}), and dispatching to
 * the reverse pass afterward.
 * <p/>
 * Extracted purely for readability, matching {@link MascotReverseModeHelper}:
 * every per-evaluation input this class needs -- the {@link
 * MascotCore.ActiveState}, the {@link MascotCore.EpochRates} array, the
 * {@link MascotCore.Workspace}, and the {@link MascotCore.OperationTapeStore}
 * -- is passed explicitly into {@link #forward}, rather than this class
 * reaching back into a shared {@code MascotCore} instance for them.
 * <p/>
 * Constructed once per {@code MascotCore} instance with that instance's fixed
 * dimensional constants and epoch boundaries (mirroring {@code MascotCore}'s
 * own fields), including the {@code epochCursor} scan position: that state is
 * exclusively a forward-pass concern (see {@code epochAt}/{@code
 * nextBoundaryAfter}'s doc), so it lives here now rather than on {@code
 * MascotCore}. Not thread-safe to call concurrently against the same {@code
 * ActiveState}/{@code Workspace}, for the same reason {@link MascotCore}
 * itself isn't (see its class doc).
 */
final class MascotForwardModeHelper {

    private final int stateCount;
    private final double maxStep;
    private final double[] boundaries;
    private final int epochCount;
    private int epochCursor;

    MascotForwardModeHelper(int stateCount, double maxStep, double[] boundaries, int epochCount) {
        this.stateCount = stateCount;
        this.maxStep = maxStep;
        this.boundaries = boundaries;
        this.epochCount = epochCount;
    }

    void forward(MascotCore.ActiveState state, MascotCore.EpochRates[] epochRates, MascotCore.Workspace workspace,
                MascotCore.PreparedOperations prepared, double[] branchRates, boolean checkProbabilities,
                MascotCore.OperationTapeStore operations, double[] nodeLogWeights) {
        BranchIntervalOperation[] branchOperations = prepared.operations;
        state.reset(stateCount, prepared.nodeCount, prepared.maxLineageId);
        if (workspace.coalescentTimes == null || workspace.coalescentTimes.length < prepared.nodeCount) {
            workspace.coalescentTimes = new double[prepared.nodeCount];
        }
        Arrays.fill(workspace.coalescentTimes, 0, prepared.nodeCount, Double.NaN);
        epochCursor = 0;
        double currentTime = prepared.initialTime;

        for (int interval = 0; interval < prepared.intervalStarts.length - 1; interval++) {
            int start = prepared.intervalStarts[interval];
            int end = prepared.intervalStarts[interval + 1];
            double intervalLength = branchOperations[start].intervalLength;
            if (intervalLength < 0.0 || !Double.isFinite(intervalLength)) {
                throw new IllegalArgumentException("invalid interval length: " + intervalLength);
            }
            for (int i = start + 1; i < end; i++) {
                if (Math.abs(branchOperations[i].intervalLength - intervalLength) > MascotCore.TIME_TOLERANCE) {
                    throw new IllegalArgumentException("operations in one interval have different lengths");
                }
            }
            for (int i = start; i < end; i++) {
                BranchIntervalOperation operation = branchOperations[i];
                if (operation.inputBuffer2 < 0 &&
                        nodeFromBuffer(operation.outputBuffer, prepared.nodeCount) !=
                                nodeFromBuffer(operation.inputBuffer1, prepared.nodeCount)) {
                    throw new IllegalArgumentException("propagation operation changed lineage identity");
                }
            }

            addNewSamples(state, prepared, branchOperations, start, end, operations, workspace);

            double intervalEnd = currentTime + intervalLength;
            while (intervalEnd > currentTime + MascotCore.TIME_TOLERANCE) {
                double segmentEnd = nextBoundaryAfter(currentTime, intervalEnd);
                int epoch = epochAt(currentTime + MascotCore.TIME_TOLERANCE);
                integrateSegment(state, epochRates[epoch], currentTime, segmentEnd, epoch, branchRates, operations,
                        workspace);
                currentTime = segmentEnd;
                if (checkProbabilities) {
                    checkProbabilities(state);
                }
            }
            currentTime = intervalEnd;

            BranchIntervalOperation coalescent = null;
            for (int i = start; i < end; i++) {
                BranchIntervalOperation operation = branchOperations[i];
                if (operation.inputBuffer2 >= 0) {
                    if (coalescent != null) {
                        throw new IllegalArgumentException("multiple coalescent operations in one interval");
                    }
                    coalescent = operation;
                }
            }
            if (coalescent != null) {
                int epoch = epochAt(currentTime);
                int child1 = nodeFromBuffer(coalescent.inputBuffer1, prepared.nodeCount);
                int child2 = nodeFromBuffer(coalescent.inputBuffer2, prepared.nodeCount);
                int parent = nodeFromBuffer(coalescent.outputBuffer, prepared.nodeCount);
                if (sameTime(workspace.coalescentTimes[child1], currentTime) ||
                        sameTime(workspace.coalescentTimes[child2], currentTime)) {
                    throw new IllegalArgumentException("dependent coalescent events at the same time are not " +
                            "currently supported");
                }
                applyCoalescentEvent(state, child1, child2, parent,
                        epoch, epochRates[epoch], operations, nodeLogWeights, workspace);
                workspace.coalescentTimes[parent] = currentTime;
            }

            if (checkProbabilities) {
                checkProbabilities(state);
            }
        }
        if (state.activeCount != 1) {
            throw new IllegalArgumentException("operation traversal ended with " + state.activeCount +
                    " active lineages");
        }
    }

    // ------------------------------------------------------------------
    // Forward integration
    // ------------------------------------------------------------------

    private void integrateSegment(MascotCore.ActiveState state, MascotCore.EpochRates rates, double start,
                                  double end, int epoch, double[] branchRates,
                                  MascotCore.OperationTapeStore operations, MascotCore.Workspace workspace) {
        if (!(end > start)) {
            throw new IllegalArgumentException("empty integration segment");
        }

        int steps = Math.max(1, (int) Math.ceil((end - start) / maxStep));
        double h = (end - start) / steps;
        int activeCount = state.activeCount;
        int dim = activeCount * stateCount + 1;

        workspace.integrationState = MascotCore.ensure(workspace.integrationState, dim);
        double[] y = workspace.integrationState;
        packStateInto(state, y);

        workspace.integrationOut = MascotCore.ensure(workspace.integrationOut, dim);
        double[] yOut = workspace.integrationOut;

        // The active lineage set (hence which clock rate applies to each ODE
        // slice) is fixed for the whole segment, so this snapshot is built once
        // here rather than once per RK4 stage.
        double[] activeClockRates = null;
        if (branchRates != null) {
            workspace.activeClockRates = MascotCore.ensure(workspace.activeClockRates, activeCount);
            activeClockRates = workspace.activeClockRates;
            for (int i = 0; i < activeCount; i++) {
                activeClockRates[i] = branchRates[state.activeIds[i]];
            }
        }

        if (operations == null) {
            for (int i = 0; i < steps; i++) {
                rk4StepInto(y, activeCount, rates, activeClockRates, h, yOut, workspace);
                double[] swap = y;
                y = yOut;
                yOut = swap;
            }
        } else {
            MascotCore.IntervalTape tape = operations.addInterval(steps, activeCount, dim, epoch, h,
                    state.activeIds, activeClockRates);
            for (int i = 0; i < steps; i++) {
                rk4StepWithTapeInto(y, activeCount, rates, activeClockRates, h, yOut, workspace, tape, i);
                double[] swap = y;
                y = yOut;
                yOut = swap;
            }
        }

        unpackStateFrom(y, state);
    }

    private void rk4StepInto(double[] y, int activeCount, MascotCore.EpochRates rates, double[] activeClockRates,
                             double h, double[] yOut, MascotCore.Workspace w) {
        int dim = activeCount * stateCount + 1;

        w.k1 = MascotCore.ensure(w.k1, dim);
        rhsInto(y, activeCount, rates, activeClockRates, w.k1, w);
        w.y2 = MascotCore.ensure(w.y2, dim);
        MascotCore.addScaledInto(y, w.k1, 0.5 * h, w.y2, dim);
        w.k2 = MascotCore.ensure(w.k2, dim);
        rhsInto(w.y2, activeCount, rates, activeClockRates, w.k2, w);
        w.y3 = MascotCore.ensure(w.y3, dim);
        MascotCore.addScaledInto(y, w.k2, 0.5 * h, w.y3, dim);
        w.k3 = MascotCore.ensure(w.k3, dim);
        rhsInto(w.y3, activeCount, rates, activeClockRates, w.k3, w);
        w.y4 = MascotCore.ensure(w.y4, dim);
        MascotCore.addScaledInto(y, w.k3, h, w.y4, dim);
        w.k4 = MascotCore.ensure(w.k4, dim);
        rhsInto(w.y4, activeCount, rates, activeClockRates, w.k4, w);

        for (int i = 0; i < dim; i++) {
            yOut[i] = y[i] + (h / 6.0) * (w.k1[i] + 2.0 * w.k2[i] + 2.0 * w.k3[i] + w.k4[i]);
        }
    }

    private void rk4StepWithTapeInto(double[] y, int activeCount, MascotCore.EpochRates rates,
                                     double[] activeClockRates, double h, double[] yOut, MascotCore.Workspace w,
                                     MascotCore.IntervalTape tape, int stepIndex) {
        int dim = tape.stateDimension;
        int offset = stepIndex * dim;
        System.arraycopy(y, 0, tape.y0, offset, dim);

        w.k1 = MascotCore.ensure(w.k1, dim);
        rhsInto(y, activeCount, rates, activeClockRates, w.k1, w);
        w.y2 = MascotCore.ensure(w.y2, dim);
        MascotCore.addScaledInto(y, w.k1, 0.5 * h, w.y2, dim);
        System.arraycopy(w.y2, 0, tape.y2, offset, dim);

        w.k2 = MascotCore.ensure(w.k2, dim);
        rhsInto(w.y2, activeCount, rates, activeClockRates, w.k2, w);
        w.y3 = MascotCore.ensure(w.y3, dim);
        MascotCore.addScaledInto(y, w.k2, 0.5 * h, w.y3, dim);
        System.arraycopy(w.y3, 0, tape.y3, offset, dim);

        w.k3 = MascotCore.ensure(w.k3, dim);
        rhsInto(w.y3, activeCount, rates, activeClockRates, w.k3, w);
        w.y4 = MascotCore.ensure(w.y4, dim);
        MascotCore.addScaledInto(y, w.k3, h, w.y4, dim);
        System.arraycopy(w.y4, 0, tape.y4, offset, dim);

        w.k4 = MascotCore.ensure(w.k4, dim);
        rhsInto(w.y4, activeCount, rates, activeClockRates, w.k4, w);

        for (int i = 0; i < dim; i++) {
            yOut[i] = y[i] + (h / 6.0) * (w.k1[i] + 2.0 * w.k2[i] + 2.0 * w.k3[i] + w.k4[i]);
        }
    }

    private void rhsInto(double[] y, int activeCount, MascotCore.EpochRates rates, double[] activeClockRates,
                         double[] out, MascotCore.Workspace w) {
        int K = stateCount;
        int stateSize = activeCount * K;

        w.sums = MascotCore.ensure(w.sums, K);
        w.sumsSquares = MascotCore.ensure(w.sumsSquares, K);

        // sums/sumsSquares accumulation is fused into the migration loop below
        // (same lineage-outer traversal, reusing the already-loaded p0/p values)
        // rather than a separate pass over activeCount * K. The migration loop's
        // own zero-fill-avoidance ("=" for source 0, "+=" for source 1..K-1 into
        // out) is preserved unchanged. sums/sumsSquares themselves are never
        // clock-scaled: they feed the coalescent-hazard term below, which is a
        // function of Ne only, per BASTA's convention that a branch clock scales
        // the migration/transition process but never the Ne-derived rate.
        double c0 = activeClockRates == null ? 1.0 : activeClockRates[0];
        double p0 = y[0];
        w.sums[0] = p0;
        w.sumsSquares[0] = p0 * p0;
        for (int sink = 0; sink < K; sink++) {
            out[sink] = c0 * p0 * rates.migrationMatrix[sink];
        }
        for (int source = 1; source < K; source++) {
            double p = y[source];
            w.sums[source] = p;
            w.sumsSquares[source] = p * p;
            int row = source * K;
            for (int sink = 0; sink < K; sink++) {
                out[sink] += c0 * p * rates.migrationMatrix[row + sink];
            }
        }

        for (int lineage = 1; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            double c = activeClockRates == null ? 1.0 : activeClockRates[lineage];
            p0 = y[offset];
            w.sums[0] += p0;
            w.sumsSquares[0] += p0 * p0;
            for (int sink = 0; sink < K; sink++) {
                out[offset + sink] = c * p0 * rates.migrationMatrix[sink];
            }
            for (int source = 1; source < K; source++) {
                double p = y[offset + source];
                w.sums[source] += p;
                w.sumsSquares[source] += p * p;
                int row = source * K;
                for (int sink = 0; sink < K; sink++) {
                    out[offset + sink] += c * p * rates.migrationMatrix[row + sink];
                }
            }
        }

        double hazard = 0.0;
        for (int state = 0; state < K; state++) {
            hazard += 0.5 * rates.inversePopulation[state] * (w.sums[state] * w.sums[state] - w.sumsSquares[state]);
        }

        w.hValues = MascotCore.ensure(w.hValues, stateSize);
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
    // Events
    // ------------------------------------------------------------------

    private void addNewSamples(MascotCore.ActiveState state, MascotCore.PreparedOperations prepared,
                               BranchIntervalOperation[] branchOperations, int start, int end,
                               MascotCore.OperationTapeStore operations, MascotCore.Workspace workspace) {
        if (workspace.sampleInputs == null || workspace.sampleInputs.length < prepared.nodeCount) {
            workspace.sampleInputs = new boolean[prepared.nodeCount];
        }
        Arrays.fill(workspace.sampleInputs, 0, prepared.nodeCount, false);
        for (int i = start; i < end; i++) {
            BranchIntervalOperation operation = branchOperations[i];
            workspace.sampleInputs[nodeFromBuffer(operation.inputBuffer1, prepared.nodeCount)] = true;
            if (operation.inputBuffer2 >= 0) {
                workspace.sampleInputs[nodeFromBuffer(operation.inputBuffer2, prepared.nodeCount)] = true;
            }
        }
        for (int node = 0; node < prepared.nodeCount; node++) {
            if (!workspace.sampleInputs[node] || state.isActive(node)) {
                continue;
            }
            if (node >= prepared.tipPartials.length) {
                throw new IllegalArgumentException("inactive internal lineage appears as operation input: " + node);
            }
            applySampleEvent(state, node, prepared.tipPartials[node], operations);
        }
    }

    private static int nodeFromBuffer(int buffer, int nodeCount) {
        if (buffer < 0) {
            throw new IllegalArgumentException("negative operation buffer: " + buffer);
        }
        return buffer % nodeCount;
    }

    private static boolean sameTime(double first, double second) {
        return !Double.isNaN(first) && Math.abs(first - second) <= MascotCore.TIME_TOLERANCE;
    }

    private void applySampleEvent(MascotCore.ActiveState state, int lineage, double[] stateProbabilities,
                                  MascotCore.OperationTapeStore operations) {
        int index = state.activeCount;
        state.ensureCapacity(index + 1);
        writeNormalizedSampleProbabilities(stateProbabilities, state.probabilities, index * stateCount);
        state.activeIds[index] = lineage;
        state.setActiveIndex(lineage, index);
        state.activeCount++;

        if (operations != null) {
            operations.addSample(index);
        }
    }

    private void applyCoalescentEvent(MascotCore.ActiveState state, int child1, int child2, int parent, int epoch,
                                      MascotCore.EpochRates rates, MascotCore.OperationTapeStore operations,
                                      double[] nodeLogWeights, MascotCore.Workspace workspace) {
        int first = state.activeIndexOf(child1);
        int second = state.activeIndexOf(child2);
        if (first < 0 || second < 0) {
            throw new IllegalArgumentException("coalescent children are not both active: " +
                    child1 + ", " + child2);
        }
        if (first == second) {
            throw new IllegalArgumentException("coalescent children must be distinct");
        }
        if (state.isActive(parent)) {
            throw new IllegalArgumentException("coalescent parent is already active: " + parent);
        }

        double[] q = rates.inversePopulation;

        workspace.coalP1 = MascotCore.ensure(workspace.coalP1, stateCount);
        workspace.coalP2 = MascotCore.ensure(workspace.coalP2, stateCount);
        double[] p1 = workspace.coalP1;
        double[] p2 = workspace.coalP2;
        System.arraycopy(state.probabilities, first * stateCount, p1, 0, stateCount);
        System.arraycopy(state.probabilities, second * stateCount, p2, 0, stateCount);

        workspace.coalParent = MascotCore.ensure(workspace.coalParent, stateCount);
        double[] parentProbabilities = workspace.coalParent;

        double lambda = 0.0;
        if (nodeLogWeights == null) {
            for (int s = 0; s < stateCount; s++) {
                parentProbabilities[s] = p1[s] * p2[s] * q[s];
                lambda += parentProbabilities[s];
            }
        } else {
            // Test-only hook (see evaluateWithNodeLogWeightsForTesting): applies a
            // hypothetical per-state log weight at this event's parent node, for
            // finite-difference validation of the adjoint node-state score. Never
            // exercised by the production (XML-facing) evaluate(...) overloads,
            // which always pass nodeLogWeights == null and take the branch above
            // without paying for exp(0).
            int weightOffset = parent * stateCount;
            for (int s = 0; s < stateCount; s++) {
                double weighted = p1[s] * p2[s] * q[s] * Math.exp(nodeLogWeights[weightOffset + s]);
                parentProbabilities[s] = weighted;
                lambda += weighted;
            }
        }
        if (!(lambda > 0.0) || !Double.isFinite(lambda)) {
            throw new MascotCore.NumericalException("invalid coalescent rate: " + lambda);
        }
        for (int s = 0; s < stateCount; s++) {
            parentProbabilities[s] /= lambda;
        }

        int beforeCount = state.activeCount;
        int afterCount = beforeCount - 1;
        int lastBefore = beforeCount - 1;
        int parentIndexAfter = first == lastBefore ? second : first;
        int removedIndex = parentIndexAfter == first ? second : first;
        int movedFromIndexBefore = -1;
        int movedToIndexAfter = -1;

        // Active-lineage order carries no probability meaning. Keep the parent in
        // one removed child slot and fill only the other hole with the old last
        // lineage when needed.
        if (removedIndex != lastBefore) {
            movedFromIndexBefore = lastBefore;
            movedToIndexAfter = removedIndex;
            int movedLineage = state.activeIds[lastBefore];
            state.activeIds[movedToIndexAfter] = movedLineage;
            System.arraycopy(state.probabilities, lastBefore * stateCount,
                    state.probabilities, movedToIndexAfter * stateCount, stateCount);
            state.setActiveIndex(movedLineage, movedToIndexAfter);
        }

        state.activeIds[parentIndexAfter] = parent;
        System.arraycopy(parentProbabilities, 0, state.probabilities, parentIndexAfter * stateCount, stateCount);
        state.setActiveIndex(parent, parentIndexAfter);
        state.clearActiveIndex(child1);
        state.clearActiveIndex(child2);

        state.activeCount = afterCount;
        state.logLikelihood += Math.log(lambda);

        if (operations != null) {
            operations.addCoalescent(epoch, first, second, parentIndexAfter,
                    movedFromIndexBefore, movedToIndexAfter, parent,
                    p1, p2, parentProbabilities, lambda, stateCount);
        }
    }

    /**
     * Writes the normalized sample-state probability vector directly into {@code
     * out} at {@code offset}, avoiding the per-sample-event {@code double[]}
     * allocation an intermediate array would cost.
     */
    private void writeNormalizedSampleProbabilities(double[] stateProbabilities, double[] out, int offset) {
        if (stateProbabilities.length != stateCount) {
            throw new IllegalArgumentException("sample probability dimension mismatch");
        }
        double sum = 0.0;
        for (int s = 0; s < stateCount; s++) {
            if (stateProbabilities[s] < 0.0) {
                throw new IllegalArgumentException("sample probabilities must be nonnegative");
            }
            sum += stateProbabilities[s];
        }
        if (!(sum > 0.0)) {
            throw new IllegalArgumentException("sample probabilities must have positive sum");
        }
        for (int s = 0; s < stateCount; s++) {
            out[offset + s] = stateProbabilities[s] / sum;
        }
    }

    // ------------------------------------------------------------------
    // Epoch/time bookkeeping
    // ------------------------------------------------------------------

    /**
     * Both {@link #epochAt} and {@link #nextBoundaryAfter} scan {@code
     * boundaries} starting from {@link #epochCursor} instead of from the
     * beginning: {@link #forward} only ever queries these with a
     * nondecreasing sequence of times (both {@code currentTime} and {@code
     * event.time} only increase over one {@code forward} call), so {@code
     * boundaries[epochCursor] <= any time queried so far} is an invariant, and
     * every boundary at or before {@code epochCursor} is therefore already known
     * to fail both methods' "is this boundary after the query time" tests. This
     * changes only the scan's starting point, not either method's comparison
     * expressions (still using the exact same TIME_TOLERANCE arithmetic), so it
     * cannot change which epoch/boundary is returned for a given time -- only how
     * many array entries are checked to find it.
     */
    private int epochAt(double t) {
        if (t < -MascotCore.TIME_TOLERANCE) {
            throw new IllegalArgumentException("time is before zero: " + t);
        }
        while (epochCursor + 1 < boundaries.length && t >= boundaries[epochCursor + 1] - MascotCore.TIME_TOLERANCE) {
            epochCursor++;
        }
        if (epochCursor >= epochCount) {
            return epochCount - 1;
        }
        return epochCursor;
    }

    private double nextBoundaryAfter(double start, double stop) {
        for (int i = epochCursor + 1; i < boundaries.length; i++) {
            double boundary = boundaries[i];
            if (boundary > start + MascotCore.TIME_TOLERANCE) {
                if (boundary <= stop + MascotCore.TIME_TOLERANCE) {
                    return Math.min(boundary, stop);
                }
                return stop;
            }
        }
        return stop;
    }

    // ------------------------------------------------------------------
    // Pack/unpack and validation
    // ------------------------------------------------------------------

    private void packStateInto(MascotCore.ActiveState state, double[] y) {
        int stateSize = state.activeCount * stateCount;
        System.arraycopy(state.probabilities, 0, y, 0, stateSize);
        y[stateSize] = state.logLikelihood;
    }

    private void unpackStateFrom(double[] y, MascotCore.ActiveState state) {
        int stateSize = state.activeCount * stateCount;
        state.ensureCapacity(state.activeCount);
        System.arraycopy(y, 0, state.probabilities, 0, stateSize);
        state.logLikelihood = y[stateSize];
    }

    private void checkProbabilities(MascotCore.ActiveState state) {
        for (int lineage = 0; lineage < state.activeCount; lineage++) {
            double sum = 0.0;
            int offset = lineage * stateCount;
            for (int s = 0; s < stateCount; s++) {
                double p = state.probabilities[offset + s];
                if (p < -1.0e-8 || !Double.isFinite(p)) {
                    throw new MascotCore.NumericalException("invalid probability for lineage " +
                            state.activeIds[lineage] + ": " + p);
                }
                sum += p;
            }
            if (Math.abs(sum - 1.0) > 1.0e-6) {
                throw new MascotCore.NumericalException("lineage probabilities do not sum to one for lineage " +
                        state.activeIds[lineage] + ": " + sum);
            }
        }
    }
}
