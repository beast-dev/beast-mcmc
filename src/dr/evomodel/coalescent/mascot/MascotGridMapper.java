/*
 * MascotGridMapper.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

/**
 * Exact draws from Pi_h, the frozen-environment surrogate discretized on MASCOT's own RK4 grid: on each step (cell)
 * lineage i's killing is constant, lambda_i(s) = invN_e(s) sum_{j != i} ybar_j(s) with ybar_j the RK4-weighted stage
 * average (y0 + 2 y2 + 2 y3 + y4) / 6. The normalized post-order filter of Pi_h is MASCOT's lineage probabilities up to
 * O(h^2), and Pi_h converges to the continuous surrogate at O(h^2) as maxStep shrinks.
 * <p>
 * Within a cell the generator is constant, so a draw is exact for Pi_h: a forward filter over cell boundaries
 * (post-order), the root drawn from its filter, then for each branch from the top down a joint draw of the state at
 * each cell's lower boundary and the path through the cell ({@link KilledUniformizedBridge}). Cells with
 * mu h > {@link KilledUniformizedBridge#MAX_RECOMMENDED_RATE_TIMES_LENGTH} are split into equal pieces.
 */
final class MascotGridMapper extends MascotStochasticMapper {

    // Cell payload: killing rates [0, K), Pi_h filter at the cell's lower boundary [K, 2K).
    private final double[] payload;
    private final double[] sums;
    private double[] averages = new double[0];
    private double[] filter = new double[0];
    private double[] nodeFilter = new double[0];

    MascotGridMapper(int stateCount) {
        super(stateCount, 2 * stateCount);
        this.payload = new double[2 * stateCount];
        this.sums = new double[stateCount];
    }

    @Override
    protected StochasticMapSampler sampler() {
        return StochasticMapSampler.GRID;
    }

    @Override
    protected void beginReplay() {
        filter = ensure(filter, nodeCount * stateCount);
        nodeFilter = ensure(nodeFilter, nodeCount * stateCount);
    }

    @Override
    protected void onSample(int lineage, MascotPreparedInput input) {
        int K = stateCount;
        input.tipData.writeTipPartials(lineage, filter, lineage * K);
        normalize(filter, lineage * K, "tip partial of node " + lineage);
        System.arraycopy(filter, lineage * K, nodeFilter, lineage * K, K);
    }

    @Override
    protected void onInterval(MascotAdjointTape.Interval interval, int[] slotLineage) {
        int K = stateCount;
        averages = ensure(averages, interval.activeCount * K);
        double h = interval.h;
        for (int step = 0; step < interval.steps; step++) {
            stageAverages(interval, step, averages, sums);
            double stepStart = interval.startTime + step * h;
            for (int slot = 0; slot < interval.activeCount; slot++) {
                int lineage = slotLineage[slot];
                double clockRate = interval.clockRates == null ? 1.0 : interval.clockRates[slot];
                stepKilling(interval, slot, sums, averages, payload, 0);
                double muh = maxExitRate(interval.epoch, clockRate, payload, 0) * h;
                maxRateTimesLength = Math.max(maxRateTimesLength, muh);
                int pieces = muh > KilledUniformizedBridge.MAX_RECOMMENDED_RATE_TIMES_LENGTH ?
                        (int) Math.ceil(muh / KilledUniformizedBridge.MAX_RECOMMENDED_RATE_TIMES_LENGTH) : 1;
                double pieceLength = h / pieces;
                buildGenerator(interval.epoch, clockRate, payload, 0);
                bridge.setGenerator(generator, pieceLength);
                for (int piece = 0; piece < pieces; piece++) {
                    System.arraycopy(filter, lineage * K, payload, K, K);
                    cells.add(lineage, interval.epoch, stepStart + piece * pieceLength, pieceLength, clockRate,
                            payload, 0);
                    bridge.transitionRow(filter, lineage * K, filter, lineage * K);
                    normalize(filter, lineage * K, "filter of node " + lineage + " near height " + stepStart);
                }
            }
        }
    }

    @Override
    protected void onCoalescent(int parent, int child1, int child2, MascotAdjointTape.Coalescent coalescent) {
        int K = stateCount;
        int qOffset = coalescent.epoch * K;
        for (int s = 0; s < K; s++) {
            filter[parent * K + s] = filter[child1 * K + s] * filter[child2 * K + s] * inversePopulation[qOffset + s];
        }
        normalize(filter, parent * K, "coalescence at node " + parent);
        System.arraycopy(filter, parent * K, nodeFilter, parent * K, K);
    }

    @Override
    protected int drawRootState() {
        return drawCategorical(nodeFilter, rootNode * stateCount);
    }

    @Override
    protected int drawBranch(int node, int topState, KilledUniformizedBridge.Path changes) {
        int K = stateCount;
        double[] cellPayload = cells.payload();
        int first = cells.first(node);
        int state = topState;
        for (int cell = first + cells.countOf(node) - 1; cell >= first; cell--) {
            int offset = cells.payloadOffset(cell);
            buildGenerator(cells.epoch(cell), cells.clock(cell), cellPayload, offset);
            bridge.setGenerator(generator, cells.length(cell));
            state = bridge.sampleStartAndPath(cellPayload, offset + K, state, cells.start(cell), cells.end(cell),
                    changes);
        }
        return state;
    }
}
