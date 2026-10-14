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
 * Draws node locations and complete branch histories jointly from Pi_h: MASCOT's frozen-environment surrogate,
 * discretized on MASCOT's own RK4 grid.
 * <p>
 * Lineage i is a killed continuous-time Markov chain in height whose generator on the RK4 step (cell) [t, t + h] is
 * c_i M_e - diag(lambda_i), where lambda_i(s) = invN_e(s) sum_{j != i} ybar_j(s) and ybar_j is lineage j's
 * RK4-weighted stage average (y0 + 2 y2 + 2 y3 + y4) / 6 over that step. A coalescence in state s carries weight
 * invN_e(s); tips carry their tip partials; there is no root prior. The normalized post-order filter of this model is
 * MASCOT's lineage probabilities up to O(h^2), and the model converges to the continuous frozen-environment
 * surrogate at O(h^2) as maxStep shrinks.
 * <p>
 * Within a cell the generator is constant, so a draw is exact for Pi_h: a forward filter over cell boundaries
 * (post-order, tips to root), the root drawn from its filter, then for each branch from the top down a joint draw
 * of the state at each cell's lower boundary and the path through the cell ({@link KilledUniformizedBridge}).
 * <p>
 * {@link #prepare} replays a taped forward pass and copies everything the draw needs, so later likelihood or
 * gradient evaluations that reuse the tape do not affect {@link #draw}.
 */
final class MascotStochasticMapper {

    private static final double SNAP_TOLERANCE = 1.0e-9;

    private final int stateCount;
    private final KilledUniformizedBridge bridge;
    private final KilledUniformizedBridge.Path cellPath = new KilledUniformizedBridge.Path();
    private final double[] generator;

    private double[] migration = new double[0];
    private double[] inversePopulation = new double[0];

    private int nodeCount;
    private int rootNode;
    private double[] nodeHeights = new double[0];
    private int[] parentOf = new int[0];
    private double[] nodeFilter = new double[0];
    private int coalescentCount;
    private int[] coalescentParent = new int[0];
    private int[] coalescentChild1 = new int[0];
    private int[] coalescentChild2 = new int[0];

    private int cellCount;
    private int[] branchCellStart = new int[0];
    private int[] branchCellCount = new int[0];
    private int[] cellEpoch = new int[0];
    private double[] cellStartTime = new double[0];
    private double[] cellEndTime = new double[0];
    private double[] cellLength = new double[0];
    private double[] cellClockRate = new double[0];
    private double[] cellKilling = new double[0];
    private double[] cellFilterStart = new double[0];

    // Cells in tape order, before grouping by branch.
    private int rawCount;
    private int[] rawBranch = new int[0];
    private int[] rawEpoch = new int[0];
    private double[] rawStartTime = new double[0];
    private double[] rawLength = new double[0];
    private double[] rawClockRate = new double[0];
    private double[] rawKilling = new double[0];
    private double[] rawFilterStart = new double[0];

    private double logLikelihood;
    private double maxRateTimesLength;
    private boolean prepared;

    private double[] stageWork = new double[0];
    private final double[] sums;
    private final double[] killing;
    private final KilledUniformizedBridge.Path branchPath = new KilledUniformizedBridge.Path();
    private int[] chunkStarts = new int[16];

    MascotStochasticMapper(int stateCount) {
        this.stateCount = stateCount;
        this.bridge = new KilledUniformizedBridge(stateCount);
        this.generator = new double[stateCount * stateCount];
        this.sums = new double[stateCount];
        this.killing = new double[stateCount];
    }

    boolean isPrepared() {
        return prepared;
    }

    /** Marks any earlier preparation as unusable, e.g. before a forward pass that may fail. */
    void invalidate() {
        prepared = false;
    }

    /**
     * Replays a forward pass recorded on {@code tape} into per-branch cells and runs the Pi_h post-order filter.
     *
     * @param nodeHeights   tree height of every node, by node number
     * @param logLikelihood the forward pass's log likelihood, reported with each draw
     */
    void prepare(MascotAdjointTape.Store tape, MascotRuntime.EpochRates[] epochRates, int epochCount,
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
        nodeFilter = ensure(nodeFilter, nodeCount * K);
        coalescentParent = ensureInt(coalescentParent, nodeCount);
        coalescentChild1 = ensureInt(coalescentChild1, nodeCount);
        coalescentChild2 = ensureInt(coalescentChild2, nodeCount);
        coalescentCount = 0;
        rawCount = 0;
        maxRateTimesLength = 0.0;
        this.logLikelihood = logLikelihood;

        // filter[b K + s]: current Pi_h filter of lineage b, normalized; slotLineage mirrors MASCOT's active slots.
        double[] filter = new double[nodeCount * K];
        int[] slotLineage = new int[nodeCount];
        int activeCount = 0;

        for (int index = 0; index < tape.size(); index++) {
            MascotAdjointTape.Operation operation = tape.get(index);
            if (operation instanceof MascotAdjointTape.Sample) {
                MascotAdjointTape.Sample sample = (MascotAdjointTape.Sample) operation;
                if (sample.sampleIndexAfter != activeCount) {
                    throw new IllegalStateException("sample was not appended at the final active slot");
                }
                int lineage = sample.lineageId;
                slotLineage[activeCount++] = lineage;
                input.tipData.writeTipPartials(lineage, filter, lineage * K);
                normalize(filter, lineage * K, "tip partial of node " + lineage);
                System.arraycopy(filter, lineage * K, nodeFilter, lineage * K, K);
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
                replayInterval(interval, slotLineage, filter);
            } else if (operation instanceof MascotAdjointTape.Coalescent) {
                MascotAdjointTape.Coalescent coalescent = (MascotAdjointTape.Coalescent) operation;
                int child1 = slotLineage[coalescent.child1Index];
                int child2 = slotLineage[coalescent.child2Index];
                int parent = coalescent.parentLineageId;
                int qOffset = coalescent.epoch * K;
                for (int s = 0; s < K; s++) {
                    filter[parent * K + s] = filter[child1 * K + s] * filter[child2 * K + s] *
                            inversePopulation[qOffset + s];
                }
                normalize(filter, parent * K, "coalescence at node " + parent);
                System.arraycopy(filter, parent * K, nodeFilter, parent * K, K);
                parentOf[child1] = parent;
                parentOf[child2] = parent;
                coalescentParent[coalescentCount] = parent;
                coalescentChild1[coalescentCount] = child1;
                coalescentChild2[coalescentCount] = child2;
                coalescentCount++;
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

        groupCellsByBranch();
        prepared = true;
    }

    /** One joint draw of all node states and branch histories. */
    MascotStochasticMap draw() {
        if (!prepared) {
            throw new IllegalStateException("prepare must be called before draw");
        }
        int[] states = new int[nodeCount];
        Arrays.fill(states, -1);
        StateHistory[] histories = new StateHistory[nodeCount];
        states[rootNode] = drawCategorical(nodeFilter, rootNode * stateCount);

        int totalJumps = 0;
        for (int i = coalescentCount - 1; i >= 0; i--) {
            histories[coalescentChild1[i]] = drawBranch(coalescentChild1[i], states);
            histories[coalescentChild2[i]] = drawBranch(coalescentChild2[i], states);
            totalJumps += histories[coalescentChild1[i]].getNumberOfJumps() +
                    histories[coalescentChild2[i]].getNumberOfJumps();
        }
        return new MascotStochasticMap(stateCount, states, histories, logLikelihood, maxRateTimesLength,
                cellCount, totalJumps);
    }

    private void replayInterval(MascotAdjointTape.Interval interval, int[] slotLineage, double[] filter) {
        int K = stateCount;
        int activeCount = interval.activeCount;
        int dim = interval.stateDimension;
        int epoch = interval.epoch;
        int mOffset = epoch * K * K;
        int qOffset = epoch * K;
        double h = interval.h;
        stageWork = ensure(stageWork, activeCount * K);
        double[] averages = stageWork;

        for (int step = 0; step < interval.steps; step++) {
            int offset = step * dim;
            Arrays.fill(sums, 0.0);
            for (int i = 0; i < activeCount * K; i++) {
                double average = (interval.y0[offset + i] + 2.0 * (interval.y2[offset + i] + interval.y3[offset + i]) +
                        interval.y4[offset + i]) / 6.0;
                averages[i] = average;
                sums[i % K] += average;
            }
            double stepStart = interval.startTime + step * h;

            for (int slot = 0; slot < activeCount; slot++) {
                int lineage = slotLineage[slot];
                double clockRate = interval.clockRates == null ? 1.0 : interval.clockRates[slot];
                double mu = 0.0;
                for (int s = 0; s < K; s++) {
                    killing[s] = Math.max(0.0, (sums[s] - averages[slot * K + s]) * inversePopulation[qOffset + s]);
                    mu = Math.max(mu, -clockRate * migration[mOffset + s * K + s] + killing[s]);
                }
                buildGenerator(mOffset, clockRate, killing, 0);
                double muh = mu * h;
                maxRateTimesLength = Math.max(maxRateTimesLength, muh);
                int pieces = muh > KilledUniformizedBridge.MAX_RECOMMENDED_RATE_TIMES_LENGTH ?
                        (int) Math.ceil(muh / KilledUniformizedBridge.MAX_RECOMMENDED_RATE_TIMES_LENGTH) : 1;
                double pieceLength = h / pieces;
                bridge.setGenerator(generator, pieceLength);
                for (int piece = 0; piece < pieces; piece++) {
                    addRawCell(lineage, epoch, stepStart + piece * pieceLength, pieceLength, clockRate, killing,
                            filter, lineage * K);
                    bridge.transitionRow(filter, lineage * K, filter, lineage * K);
                    normalize(filter, lineage * K, "filter of node " + lineage + " near height " + stepStart);
                }
            }
        }
    }

    private StateHistory drawBranch(int node, int[] states) {
        int K = stateCount;
        int parent = parentOf[node];
        int state = states[parent];
        int first = branchCellStart[node];
        int count = branchCellCount[node];

        branchPath.clear();
        if (chunkStarts.length < count + 1) {
            chunkStarts = new int[Math.max(count + 1, 2 * chunkStarts.length)];
        }
        // Cells are drawn from the top of the branch down; each cell's jumps are in increasing time.
        for (int k = count - 1; k >= 0; k--) {
            int cell = first + k;
            buildGenerator(cellEpoch[cell] * K * K, cellClockRate[cell], cellKilling, cell * K);
            bridge.setGenerator(generator, cellLength[cell]);
            chunkStarts[count - 1 - k] = branchPath.size();
            state = bridge.sampleStartAndPath(cellFilterStart, cell * K, state,
                    cellStartTime[cell], cellEndTime[cell], branchPath);
        }
        chunkStarts[count] = branchPath.size();
        states[node] = state;

        StateHistory history = new StateHistory(nodeHeights[node], state, K);
        int current = state;
        for (int chunk = count - 1; chunk >= 0; chunk--) {
            for (int j = chunkStarts[chunk]; j < chunkStarts[chunk + 1]; j++) {
                current = branchPath.getState(j);
                history.addChange(new StateChange(branchPath.getTime(j), current));
            }
        }
        if (current != states[parent]) {
            throw new IllegalStateException("branch above node " + node + " ends in state " + current +
                    " but its parent was drawn in state " + states[parent]);
        }
        history.addEndingState(new StateChange(nodeHeights[parent], states[parent]));
        return history;
    }

    /** generator = clockRate M_e - diag(killing), with M_e at {@code mOffset}. */
    private void buildGenerator(int mOffset, double clockRate, double[] killingRates, int killingOffset) {
        int K = stateCount;
        for (int s = 0; s < K; s++) {
            for (int t = 0; t < K; t++) {
                generator[s * K + t] = clockRate * migration[mOffset + s * K + t];
            }
            generator[s * K + s] -= killingRates[killingOffset + s];
        }
    }

    private void addRawCell(int lineage, int epoch, double startTime, double length, double clockRate,
                            double[] killingRates, double[] filter, int filterOffset) {
        int K = stateCount;
        if (rawCount == rawBranch.length) {
            int capacity = Math.max(64, 2 * rawCount);
            rawBranch = Arrays.copyOf(rawBranch, capacity);
            rawEpoch = Arrays.copyOf(rawEpoch, capacity);
            rawStartTime = Arrays.copyOf(rawStartTime, capacity);
            rawLength = Arrays.copyOf(rawLength, capacity);
            rawClockRate = Arrays.copyOf(rawClockRate, capacity);
            rawKilling = Arrays.copyOf(rawKilling, capacity * K);
            rawFilterStart = Arrays.copyOf(rawFilterStart, capacity * K);
        }
        rawBranch[rawCount] = lineage;
        rawEpoch[rawCount] = epoch;
        rawStartTime[rawCount] = startTime;
        rawLength[rawCount] = length;
        rawClockRate[rawCount] = clockRate;
        System.arraycopy(killingRates, 0, rawKilling, rawCount * K, K);
        System.arraycopy(filter, filterOffset, rawFilterStart, rawCount * K, K);
        rawCount++;
    }

    /**
     * Stable counting sort of the tape-ordered cells by branch, then snaps each branch's cell boundaries to the
     * tree: the lowest cell starts at the node's height, each cell ends where the next starts, and the highest ends
     * at the parent's height. Cell lengths used for the transition matrices are left as the forward pass had them.
     */
    private void groupCellsByBranch() {
        int K = stateCount;
        cellCount = rawCount;
        branchCellStart = ensureInt(branchCellStart, nodeCount);
        branchCellCount = ensureInt(branchCellCount, nodeCount);
        Arrays.fill(branchCellCount, 0, nodeCount, 0);
        for (int c = 0; c < rawCount; c++) {
            branchCellCount[rawBranch[c]]++;
        }
        int next = 0;
        for (int b = 0; b < nodeCount; b++) {
            branchCellStart[b] = next;
            next += branchCellCount[b];
        }

        cellEpoch = ensureInt(cellEpoch, cellCount);
        cellStartTime = ensure(cellStartTime, cellCount);
        cellEndTime = ensure(cellEndTime, cellCount);
        cellLength = ensure(cellLength, cellCount);
        cellClockRate = ensure(cellClockRate, cellCount);
        cellKilling = ensure(cellKilling, cellCount * K);
        cellFilterStart = ensure(cellFilterStart, cellCount * K);

        int[] fill = new int[nodeCount];
        for (int c = 0; c < rawCount; c++) {
            int b = rawBranch[c];
            int target = branchCellStart[b] + fill[b]++;
            cellEpoch[target] = rawEpoch[c];
            cellStartTime[target] = rawStartTime[c];
            cellLength[target] = rawLength[c];
            cellClockRate[target] = rawClockRate[c];
            System.arraycopy(rawKilling, c * K, cellKilling, target * K, K);
            System.arraycopy(rawFilterStart, c * K, cellFilterStart, target * K, K);
        }

        for (int b = 0; b < nodeCount; b++) {
            if (b == rootNode) {
                if (branchCellCount[b] != 0) {
                    throw new IllegalStateException("root node " + b + " has cells above it");
                }
                continue;
            }
            int parent = parentOf[b];
            if (parent < 0) {
                throw new IllegalStateException("node " + b + " has no parent in the tape replay");
            }
            double bottom = nodeHeights[b];
            double top = nodeHeights[parent];
            int count = branchCellCount[b];
            if (count == 0) {
                requireClose(bottom, top, "zero-cell branch above node " + b);
                continue;
            }
            int first = branchCellStart[b];
            int last = first + count - 1;
            requireClose(cellStartTime[first], bottom, "lowest cell of the branch above node " + b);
            requireClose(cellStartTime[last] + cellLength[last], top, "highest cell of the branch above node " + b);
            cellStartTime[first] = bottom;
            for (int c = first; c < last; c++) {
                requireClose(cellStartTime[c] + cellLength[c], cellStartTime[c + 1],
                        "consecutive cells of the branch above node " + b);
                cellEndTime[c] = cellStartTime[c + 1];
            }
            cellEndTime[last] = top;
            for (int c = first; c <= last; c++) {
                if (!(cellEndTime[c] > cellStartTime[c])) {
                    throw new IllegalStateException("empty cell in the branch above node " + b);
                }
            }
        }
    }

    private int drawCategorical(double[] weights, int offset) {
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

    private void normalize(double[] values, int offset, String what) {
        double sum = 0.0;
        for (int s = 0; s < stateCount; s++) {
            sum += values[offset + s];
        }
        if (!(sum > 0.0) || !Double.isFinite(sum)) {
            throw new IllegalStateException("Pi_h filter has total mass " + sum + " at the " + what);
        }
        for (int s = 0; s < stateCount; s++) {
            values[offset + s] /= sum;
        }
    }

    private static void requireClose(double a, double b, String what) {
        if (Math.abs(a - b) > SNAP_TOLERANCE * (1.0 + Math.abs(b))) {
            throw new IllegalStateException("time mismatch at the " + what + ": " + a + " vs " + b);
        }
    }

    private static double[] ensure(double[] array, int size) {
        return array.length >= size ? array : new double[size];
    }

    private static int[] ensureInt(int[] array, int size) {
        return array.length >= size ? array : new int[size];
    }
}
