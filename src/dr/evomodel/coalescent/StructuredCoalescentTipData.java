/*
 * StructuredCoalescentTipData.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent;

import java.util.Arrays;

/**
 * Engine-neutral fixed tip payload for structured-coalescent schedules.
 *
 * Sample timing and parent/child order live in {@link StructuredCoalescentSchedule};
 * this class only owns the state information attached to sampled nodes.
 */
public final class StructuredCoalescentTipData {

    public final int nodeCount;
    public final int stateCount;

    private final double[][] tipPartials;
    private final int[] observedStates;

    public StructuredCoalescentTipData(int nodeCount, int stateCount,
                                       double[][] tipPartials, int[] observedStates) {
        if (nodeCount <= 0) {
            throw new IllegalArgumentException("nodeCount must be positive");
        }
        if (stateCount <= 0) {
            throw new IllegalArgumentException("stateCount must be positive");
        }
        if (tipPartials == null && observedStates == null) {
            throw new IllegalArgumentException("at least one tip data source is required");
        }
        if (tipPartials != null && tipPartials.length > nodeCount) {
            throw new IllegalArgumentException("tip partials length " + tipPartials.length +
                    " exceeds nodeCount " + nodeCount);
        }
        if (observedStates != null && observedStates.length != nodeCount) {
            throw new IllegalArgumentException("observedStates length " + observedStates.length +
                    " does not match nodeCount " + nodeCount);
        }

        this.nodeCount = nodeCount;
        this.stateCount = stateCount;
        this.tipPartials = cloneNormalizedPartials(tipPartials, nodeCount, stateCount);
        this.observedStates = observedStates == null ? filledObservedStates(nodeCount) : observedStates.clone();

        for (int node = 0; node < this.observedStates.length; node++) {
            int state = this.observedStates[node];
            if (state != -1 && (state < 0 || state >= stateCount)) {
                throw new IllegalArgumentException("observed state for node " + node + " out of range: " +
                        state + " (stateCount=" + stateCount + ")");
            }
        }
    }

    public static StructuredCoalescentTipData fromPartials(int nodeCount, int stateCount, double[][] tipPartials) {
        return new StructuredCoalescentTipData(nodeCount, stateCount, tipPartials, null);
    }

    public void validateCompatibleWith(StructuredCoalescentSchedule schedule) {
        if (schedule == null) {
            throw new IllegalArgumentException("schedule must not be null");
        }
        if (schedule.nodeCount != nodeCount) {
            throw new IllegalArgumentException("tip data nodeCount " + nodeCount +
                    " does not match schedule nodeCount " + schedule.nodeCount);
        }
        validateSampleCovered(schedule.initialSampleNode);
        for (int i = 0; i < schedule.getIntervalCount(); i++) {
            if (schedule.eventTypes[i] == StructuredCoalescentSchedule.SAMPLE) {
                validateSampleCovered(schedule.sampleNodes[i]);
            }
        }
    }

    public boolean hasTipPartials(int node) {
        validateNode(node);
        return tipPartials[node] != null;
    }

    public boolean hasObservedState(int node) {
        validateNode(node);
        return observedStates[node] >= 0;
    }

    public double[] getTipPartialsCopy(int node) {
        validateSampleCovered(node);
        double[] copy = new double[stateCount];
        writeTipPartials(node, copy, 0);
        return copy;
    }

    public void writeTipPartials(int node, double[] out, int offset) {
        validateNode(node);
        if (out == null || offset < 0 || offset + stateCount > out.length) {
            throw new IllegalArgumentException("invalid output array or offset");
        }
        double[] partials = tipPartials[node];
        if (partials != null) {
            System.arraycopy(partials, 0, out, offset, stateCount);
            return;
        }
        int observedState = observedStates[node];
        if (observedState >= 0) {
            Arrays.fill(out, offset, offset + stateCount, 0.0);
            out[offset + observedState] = 1.0;
            return;
        }
        throw new IllegalArgumentException("missing tip state or partials for sample node " + node);
    }

    private void validateSampleCovered(int node) {
        validateNode(node);
        if (tipPartials[node] == null && observedStates[node] < 0) {
            throw new IllegalArgumentException("missing tip state or partials for sample node " + node);
        }
    }

    private void validateNode(int node) {
        if (node < 0 || node >= nodeCount) {
            throw new IllegalArgumentException("node out of range: " + node + " (nodeCount=" + nodeCount + ")");
        }
    }

    private static double[][] cloneNormalizedPartials(double[][] source, int nodeCount, int stateCount) {
        double[][] copy = new double[nodeCount][];
        if (source == null) {
            return copy;
        }
        for (int node = 0; node < source.length; node++) {
            double[] partials = source[node];
            if (partials == null) {
                continue;
            }
            if (partials.length != stateCount) {
                throw new IllegalArgumentException("tip partial dimension mismatch for node " + node +
                        ": " + partials.length + " != " + stateCount);
            }
            double sum = 0.0;
            for (int state = 0; state < stateCount; state++) {
                double value = partials[state];
                if (value < 0.0 || !Double.isFinite(value)) {
                    throw new IllegalArgumentException("invalid tip partial for node " + node +
                            ", state " + state + ": " + value);
                }
                sum += value;
            }
            if (!(sum > 0.0) || !Double.isFinite(sum)) {
                throw new IllegalArgumentException("tip partials must have positive finite sum for node " + node);
            }
            copy[node] = new double[stateCount];
            for (int state = 0; state < stateCount; state++) {
                copy[node][state] = partials[state] / sum;
            }
        }
        return copy;
    }

    private static int[] filledObservedStates(int nodeCount) {
        int[] states = new int[nodeCount];
        Arrays.fill(states, -1);
        return states;
    }
}
