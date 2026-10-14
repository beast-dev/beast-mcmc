/*
 * MascotStochasticMapper.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.inference.markovjumps.StateChange;
import dr.inference.markovjumps.StateHistory;
import dr.math.MathUtils;

import java.util.Arrays;

/**
 * Shared machinery of the MASCOT stochastic mappers, which draw every node's location and every branch's complete
 * history jointly from MASCOT's frozen-environment surrogate: lineage i is a killed continuous-time Markov chain in
 * height with generator c_i M_e - diag(lambda_i(t)), lambda_i(s, t) = invN_e(s) sum_{j != i} p_j(s, t) from MASCOT's
 * own forward solution; a coalescence in state s carries weight invN_e(s); tips carry their tip partials; there is
 * no root prior. Branch histories are conditionally independent given the node states.
 * <p>
 * {@link #prepare} replays a taped forward pass (one cell per lineage per RK4 step) and copies everything a draw
 * needs, so later likelihood or gradient evaluations that reuse the tape do not affect {@link #draw}. A draw picks
 * the root state, then for each branch from the top down calls the sampler's {@link #drawBranch}, and assembles one
 * {@link StateHistory} per branch in increasing height.
 */
abstract class MascotStochasticMapper {

    protected final int stateCount;
    protected final MascotBranchCells cells;
    protected final KilledUniformizedBridge bridge;
    protected final double[] generator;

    protected double[] migration = new double[0];
    protected double[] inversePopulation = new double[0];
    protected int nodeCount;
    protected int rootNode;
    protected double[] nodeHeights = new double[0];
    protected double maxRateTimesLength;

    private int[] parentOf = new int[0];
    private int coalescentCount;
    private int[] coalescentParent = new int[0];
    private int[] coalescentChild1 = new int[0];
    private int[] coalescentChild2 = new int[0];
    private double logLikelihood;
    private boolean prepared;

    private final KilledUniformizedBridge.Path branchChanges = new KilledUniformizedBridge.Path();
    private double[] sortedTimes = new double[16];
    private int[] sortedStates = new int[16];

    MascotStochasticMapper(int stateCount, int payloadWidth) {
        this.stateCount = stateCount;
        this.cells = new MascotBranchCells(payloadWidth);
        this.bridge = new KilledUniformizedBridge(stateCount);
        this.generator = new double[stateCount * stateCount];
    }

    final boolean isPrepared() {
        return prepared;
    }

    /** Marks any earlier preparation as unusable, e.g. before a forward pass that may fail. */
    final void invalidate() {
        prepared = false;
    }

    /**
     * Replays a forward pass recorded on {@code tape}.
     *
     * @param nodeHeights   tree height of every node, by node number
     * @param logLikelihood the forward pass's log likelihood, reported with each draw
     */
    final void prepare(MascotAdjointTape.Store tape, MascotRuntime.EpochRates[] epochRates, int epochCount,
                       MascotPreparedInput input, double[] nodeHeights, double logLikelihood) {
        prepared = false;
        int K = stateCount;
        if (nodeHeights.length != input.nodeCount) {
            throw new IllegalArgumentException("nodeHeights length " + nodeHeights.length +
                    " does not match node count " + input.nodeCount);
        }

        migration = ensure(migration, epochCount * K * K);
        inversePopulation = ensure(inversePopulation, epochCount * K);
        for (int e = 0; e < epochCount; e++) {
            System.arraycopy(epochRates[e].migrationMatrix, 0, migration, e * K * K, K * K);
            System.arraycopy(epochRates[e].inversePopulation, 0, inversePopulation, e * K, K);
        }

        nodeCount = input.nodeCount;
        this.nodeHeights = ensure(this.nodeHeights, nodeCount);
        System.arraycopy(nodeHeights, 0, this.nodeHeights, 0, nodeCount);
        parentOf = ensureInt(parentOf, nodeCount);
        Arrays.fill(parentOf, 0, nodeCount, -1);
        coalescentParent = ensureInt(coalescentParent, nodeCount);
        coalescentChild1 = ensureInt(coalescentChild1, nodeCount);
        coalescentChild2 = ensureInt(coalescentChild2, nodeCount);
        coalescentCount = 0;
        maxRateTimesLength = 0.0;
        this.logLikelihood = logLikelihood;
        cells.reset();
        beginReplay();

        int[] slotLineage = new int[nodeCount];
        int activeCount = 0;
        for (int index = 0; index < tape.size(); index++) {
            MascotAdjointTape.Operation operation = tape.get(index);
            if (operation instanceof MascotAdjointTape.Sample) {
                MascotAdjointTape.Sample sample = (MascotAdjointTape.Sample) operation;
                if (sample.sampleIndexAfter != activeCount) {
                    throw new IllegalStateException("sample was not appended at the final active slot");
                }
                slotLineage[activeCount++] = sample.lineageId;
                onSample(sample.lineageId, input);
            } else if (operation instanceof MascotAdjointTape.Interval) {
                MascotAdjointTape.Interval interval = (MascotAdjointTape.Interval) operation;
                if (interval.activeCount != activeCount) {
                    throw new IllegalStateException("tape interval has " + interval.activeCount +
                            " active lineages; replay has " + activeCount);
                }
                if (interval.activeIds != null) {
                    for (int i = 0; i < activeCount; i++) {
                        if (interval.activeIds[i] != slotLineage[i]) {
                            throw new IllegalStateException("tape and replay disagree on active slot " + i);
                        }
                    }
                }
                onInterval(interval, slotLineage);
            } else if (operation instanceof MascotAdjointTape.Coalescent) {
                MascotAdjointTape.Coalescent coalescent = (MascotAdjointTape.Coalescent) operation;
                int child1 = slotLineage[coalescent.child1Index];
                int child2 = slotLineage[coalescent.child2Index];
                int parent = coalescent.parentLineageId;
                parentOf[child1] = parent;
                parentOf[child2] = parent;
                coalescentParent[coalescentCount] = parent;
                coalescentChild1[coalescentCount] = child1;
                coalescentChild2[coalescentCount] = child2;
                coalescentCount++;
                onCoalescent(parent, child1, child2, coalescent);
                if (coalescent.movedFromIndexBefore >= 0) {
                    slotLineage[coalescent.movedToIndexAfter] = slotLineage[coalescent.movedFromIndexBefore];
                }
                slotLineage[coalescent.parentIndexAfter] = parent;
                activeCount--;
            } else {
                throw new IllegalStateException("unknown tape operation: " + operation.getClass());
            }
        }
        if (activeCount != 1 || coalescentCount == 0) {
            throw new IllegalStateException("tape replay ended with " + activeCount + " active lineages");
        }
        rootNode = coalescentParent[coalescentCount - 1];
        cells.group(nodeCount, this.nodeHeights, parentOf, rootNode);
        prepared = true;
    }

    /** One joint draw of all node states and branch histories. */
    final MascotStochasticMap draw() {
        if (!prepared) {
            throw new IllegalStateException("prepare must be called before draw");
        }
        int[] states = new int[nodeCount];
        Arrays.fill(states, -1);
        StateHistory[] histories = new StateHistory[nodeCount];
        states[rootNode] = drawRootState();

        int totalJumps = 0;
        for (int i = coalescentCount - 1; i >= 0; i--) {
            for (int child : new int[]{coalescentChild1[i], coalescentChild2[i]}) {
                histories[child] = drawBranchHistory(child, states);
                totalJumps += histories[child].getNumberOfJumps();
            }
        }
        return new MascotStochasticMap(stateCount, states, histories, logLikelihood, maxRateTimesLength,
                cells.count(), totalJumps);
    }

    /** Called once before the tape is replayed. */
    protected void beginReplay() {
    }

    protected abstract void onSample(int lineage, MascotPreparedInput input);

    protected abstract void onInterval(MascotAdjointTape.Interval interval, int[] slotLineage);

    protected abstract void onCoalescent(int parent, int child1, int child2, MascotAdjointTape.Coalescent coalescent);

    protected abstract int drawRootState();

    /**
     * Draws the path down the branch above {@code node}, whose top is in {@code topState}, and returns the state at
     * the node. Every jump is added to {@code changes} as (height, state just above that height), in any order.
     * Called only for branches with at least one cell.
     */
    protected abstract int drawBranch(int node, int topState, KilledUniformizedBridge.Path changes);

    private StateHistory drawBranchHistory(int node, int[] states) {
        int parent = parentOf[node];
        branchChanges.clear();
        int bottom = cells.countOf(node) == 0 ? states[parent] : drawBranch(node, states[parent], branchChanges);
        states[node] = bottom;

        int m = branchChanges.size();
        if (sortedTimes.length < m) {
            sortedTimes = new double[2 * m];
            sortedStates = new int[2 * m];
        }
        for (int j = 0; j < m; j++) {
            sortedTimes[j] = branchChanges.getTime(j);
        }
        Arrays.sort(sortedTimes, 0, m);
        Arrays.fill(sortedStates, 0, m, -1);
        for (int j = 0; j < m; j++) {
            int position = Arrays.binarySearch(sortedTimes, 0, m, branchChanges.getTime(j));
            if (sortedStates[position] != -1) {
                throw new IllegalStateException("two jumps at the same height on the branch above node " + node);
            }
            sortedStates[position] = branchChanges.getState(j);
        }

        StateHistory history = new StateHistory(nodeHeights[node], bottom, stateCount);
        int current = bottom;
        for (int j = 0; j < m; j++) {
            if (sortedStates[j] == current) {
                throw new IllegalStateException("self-jump on the branch above node " + node);
            }
            current = sortedStates[j];
            history.addChange(new StateChange(sortedTimes[j], current));
        }
        if (current != states[parent]) {
            throw new IllegalStateException("branch above node " + node + " ends in state " + current +
                    " but its parent was drawn in state " + states[parent]);
        }
        history.addEndingState(new StateChange(nodeHeights[parent], states[parent]));
        return history;
    }

    /** generator = clockRate M_epoch - diag(killing at {@code killingOffset}). */
    protected final void buildGenerator(int epoch, double clockRate, double[] killing, int killingOffset) {
        int K = stateCount;
        int mOffset = epoch * K * K;
        for (int s = 0; s < K; s++) {
            for (int t = 0; t < K; t++) {
                generator[s * K + t] = clockRate * migration[mOffset + s * K + t];
            }
            generator[s * K + s] -= killing[killingOffset + s];
        }
    }

    /** Largest total exit rate (migration plus killing) of {@code clockRate M_epoch - diag(killing)}. */
    protected final double maxExitRate(int epoch, double clockRate, double[] killing, int killingOffset) {
        int K = stateCount;
        int mOffset = epoch * K * K;
        double mu = 0.0;
        for (int s = 0; s < K; s++) {
            mu = Math.max(mu, -clockRate * migration[mOffset + s * K + s] + killing[killingOffset + s]);
        }
        return mu;
    }

    /** Fills {@code killing[offset + s]} with the RK4-weighted environment of slot {@code slot} over one step. */
    protected final void stepKilling(MascotAdjointTape.Interval interval, int slot, double[] sums,
                                     double[] averages, double[] killing, int offset) {
        int K = stateCount;
        int qOffset = interval.epoch * K;
        for (int s = 0; s < K; s++) {
            killing[offset + s] = Math.max(0.0, (sums[s] - averages[slot * K + s]) * inversePopulation[qOffset + s]);
        }
    }

    /** RK4-weighted stage averages (y0 + 2 y2 + 2 y3 + y4) / 6 of every slot over one step, and their sums. */
    protected final void stageAverages(MascotAdjointTape.Interval interval, int step, double[] averages,
                                       double[] sums) {
        int K = stateCount;
        int offset = step * interval.stateDimension;
        Arrays.fill(sums, 0, K, 0.0);
        for (int i = 0; i < interval.activeCount * K; i++) {
            double average = (interval.y0[offset + i] + 2.0 * (interval.y2[offset + i] + interval.y3[offset + i]) +
                    interval.y4[offset + i]) / 6.0;
            averages[i] = average;
            sums[i % K] += average;
        }
    }

    protected final int drawCategorical(double[] weights, int offset) {
        double total = 0.0;
        for (int s = 0; s < stateCount; s++) {
            total += weights[offset + s];
        }
        if (!(total > 0.0) || !Double.isFinite(total)) {
            throw new IllegalStateException("cannot draw from weights with total " + total);
        }
        double u = MathUtils.nextDouble() * total;
        double cumulative = 0.0;
        int chosen = -1;
        for (int s = 0; s < stateCount; s++) {
            if (weights[offset + s] > 0.0) {
                chosen = s;
                cumulative += weights[offset + s];
                if (u < cumulative) {
                    break;
                }
            }
        }
        return chosen;
    }

    protected final void normalize(double[] values, int offset, String what) {
        double sum = 0.0;
        for (int s = 0; s < stateCount; s++) {
            sum += values[offset + s];
        }
        if (!(sum > 0.0) || !Double.isFinite(sum)) {
            throw new IllegalStateException("filter has total mass " + sum + " at the " + what);
        }
        for (int s = 0; s < stateCount; s++) {
            values[offset + s] /= sum;
        }
    }

    protected static double[] ensure(double[] array, int size) {
        return array.length >= size ? array : new double[size];
    }

    protected static int[] ensureInt(int[] array, int size) {
        return array.length >= size ? array : new int[size];
    }
}
