/*
 * MascotStochasticMap.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.inference.markovjumps.StateHistory;

/**
 * One joint draw of every node's location and every branch's complete location history from Pi_h, MASCOT's
 * frozen-environment surrogate on its own RK4 grid (see {@link MascotStochasticMapper}).
 * <p>
 * A branch history runs in increasing height, from the child node (at its height, in its drawn state) to the parent
 * (at its height, in its drawn state). A change to state b at height t after state a is a backward-time move a -> b,
 * i.e. a forward-time migration b -> a; {@code StateHistory.getJumpCounts()[a * K + b]} counts these.
 * <p>
 * Histories are returned as copies because {@link StateHistory} is mutable.
 */
public final class MascotStochasticMap {

    private final int stateCount;
    private final int[] nodeStates;
    private final StateHistory[] histories;
    private final double logLikelihood;
    private final double maxRateTimesLength;
    private final int cellCount;
    private final int totalJumpCount;

    MascotStochasticMap(int stateCount, int[] nodeStates, StateHistory[] histories, double logLikelihood,
                        double maxRateTimesLength, int cellCount, int totalJumpCount) {
        this.stateCount = stateCount;
        this.nodeStates = nodeStates;
        this.histories = histories;
        this.logLikelihood = logLikelihood;
        this.maxRateTimesLength = maxRateTimesLength;
        this.cellCount = cellCount;
        this.totalJumpCount = totalJumpCount;
    }

    public int getStateCount() {
        return stateCount;
    }

    public int getNodeCount() {
        return nodeStates.length;
    }

    /** The drawn location of a node (root, internal node or tip) by node number. */
    public int getNodeState(int nodeNumber) {
        return nodeStates[nodeNumber];
    }

    public int[] getNodeStates() {
        return nodeStates.clone();
    }

    /** A copy of the complete history of the branch above a node, or null for the root. */
    public StateHistory getBranchHistory(int nodeNumber) {
        StateHistory history = histories[nodeNumber];
        return history == null ? null : new StateHistory(history);
    }

    /** MASCOT log likelihood of the forward pass this draw was conditioned on. */
    public double getLogLikelihood() {
        return logLikelihood;
    }

    /**
     * Largest h * mu over all cells, where mu is the largest total exit rate (migration plus killing). Values well
     * above one mean the RK4 steps are coarse relative to the dynamics, which affects the likelihood as well.
     */
    public double getMaxRateTimesLength() {
        return maxRateTimesLength;
    }

    public int getCellCount() {
        return cellCount;
    }

    public int getTotalJumpCount() {
        return totalJumpCount;
    }
}
