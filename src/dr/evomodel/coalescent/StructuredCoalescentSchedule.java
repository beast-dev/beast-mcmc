/*
 * StructuredCoalescentSchedule.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;

import java.util.Arrays;

/**
 * Primitive, engine-neutral tree/event schedule for structured-coalescent
 * likelihoods. Each entry represents one interval of elapsed backward time,
 * followed by the biological event at the top of that interval. The most
 * recent sample is stored separately and is active at {@link #initialTime}.
 *
 * No BASTA buffers, transition-matrix identifiers, BEAGLE layout, or MASCOT
 * active-lineage slots live here. Engine-specific prepared objects should lower
 * these arrays into their own hot-path representation.
 */
public final class StructuredCoalescentSchedule {

    public static final int NO_NODE = -1;
    public static final int SAMPLE = 1;
    public static final int COALESCENT = 2;

    public final double initialTime;
    public final int initialSampleNode;
    public final int nodeCount;
    public final int maxLineageId;

    public final double[] intervalLengths;
    public final int[] eventTypes;
    public final int[] sampleNodes;
    public final int[] coalescentChild1;
    public final int[] coalescentChild2;
    public final int[] coalescentParent;

    public StructuredCoalescentSchedule(double initialTime,
                                        int initialSampleNode,
                                        int nodeCount,
                                        double[] intervalLengths,
                                        int[] eventTypes,
                                        int[] sampleNodes,
                                        int[] coalescentChild1,
                                        int[] coalescentChild2,
                                        int[] coalescentParent) {
        if (nodeCount <= 0) {
            throw new IllegalArgumentException("nodeCount must be positive");
        }
        if (initialTime < 0.0 || !Double.isFinite(initialTime)) {
            throw new IllegalArgumentException("invalid initial time: " + initialTime);
        }
        validateNode(initialSampleNode, nodeCount, "initial sample");
        int intervalCount = requireSameLength(intervalLengths, eventTypes, sampleNodes,
                coalescentChild1, coalescentChild2, coalescentParent);

        this.initialTime = initialTime;
        this.initialSampleNode = initialSampleNode;
        this.nodeCount = nodeCount;
        this.maxLineageId = nodeCount - 1;
        this.intervalLengths = intervalLengths.clone();
        this.eventTypes = eventTypes.clone();
        this.sampleNodes = sampleNodes.clone();
        this.coalescentChild1 = coalescentChild1.clone();
        this.coalescentChild2 = coalescentChild2.clone();
        this.coalescentParent = coalescentParent.clone();

        for (int i = 0; i < intervalCount; i++) {
            if (this.intervalLengths[i] < 0.0 || !Double.isFinite(this.intervalLengths[i])) {
                throw new IllegalArgumentException("invalid interval length at " + i + ": " +
                        this.intervalLengths[i]);
            }
            if (this.eventTypes[i] == SAMPLE) {
                validateNode(this.sampleNodes[i], nodeCount, "sample event " + i);
                requireNoNode(this.coalescentChild1[i], "coalescentChild1", i);
                requireNoNode(this.coalescentChild2[i], "coalescentChild2", i);
                requireNoNode(this.coalescentParent[i], "coalescentParent", i);
            } else if (this.eventTypes[i] == COALESCENT) {
                requireNoNode(this.sampleNodes[i], "sampleNode", i);
                validateNode(this.coalescentChild1[i], nodeCount, "coalescent child1 event " + i);
                validateNode(this.coalescentChild2[i], nodeCount, "coalescent child2 event " + i);
                validateNode(this.coalescentParent[i], nodeCount, "coalescent parent event " + i);
                if (this.coalescentChild1[i] == this.coalescentChild2[i]) {
                    throw new IllegalArgumentException("coalescent children must be distinct at event " + i);
                }
            } else {
                throw new IllegalArgumentException("unknown event type at " + i + ": " + this.eventTypes[i]);
            }
        }
        StructuredCoalescentScheduleWalker.validate(this);
    }

    public int getIntervalCount() {
        return intervalLengths.length;
    }

    public static StructuredCoalescentSchedule fromTreeIntervals(Tree tree,
                                                                 BigFastTreeIntervals treeIntervals,
                                                                 boolean checkBinaryTree,
                                                                 boolean checkForZeroLengthCoalescentIntervals) {
        int intervalCount = treeIntervals.getIntervalCount();
        NodeRef initialSample = treeIntervals.getSamplingNode(-1);
        int nodeCount = tree.getNodeCount();

        double[] lengths = new double[intervalCount];
        int[] types = new int[intervalCount];
        int[] samples = filledNodes(intervalCount);
        int[] child1 = filledNodes(intervalCount);
        int[] child2 = filledNodes(intervalCount);
        int[] parents = filledNodes(intervalCount);

        for (int interval = 0; interval < intervalCount; interval++) {
            lengths[interval] = treeIntervals.getInterval(interval);
            IntervalType type = treeIntervals.getIntervalType(interval);
            if (type == IntervalType.SAMPLE) {
                types[interval] = SAMPLE;
                samples[interval] = treeIntervals.getSamplingNode(interval).getNumber();
            } else if (type == IntervalType.COALESCENT) {
                NodeRef parent = treeIntervals.getCoalescentNode(interval);
                if (checkBinaryTree && tree.getChildCount(parent) != 2) {
                    throw new IllegalArgumentException("The structured coalescent requires binary trees; node " +
                            parent.getNumber() + " has " + tree.getChildCount(parent) + " children");
                }
                if (checkForZeroLengthCoalescentIntervals && lengths[interval] <= 0.0) {
                    throw new IllegalArgumentException("Cannot coalesce in <= 0.0 time");
                }
                types[interval] = COALESCENT;
                child1[interval] = tree.getChild(parent, 0).getNumber();
                child2[interval] = tree.getChild(parent, 1).getNumber();
                parents[interval] = parent.getNumber();
            } else {
                throw new IllegalArgumentException("Unknown interval type: " + type);
            }
        }

        if (intervalCount == 0 || types[intervalCount - 1] != COALESCENT) {
            throw new IllegalArgumentException("Last structured-coalescent event must be a coalescence");
        }

        return new StructuredCoalescentSchedule(tree.getNodeHeight(initialSample), initialSample.getNumber(),
                nodeCount, lengths, types, samples, child1, child2, parents);
    }

    private static int[] filledNodes(int length) {
        int[] nodes = new int[length];
        Arrays.fill(nodes, NO_NODE);
        return nodes;
    }

    private static int requireSameLength(double[] intervalLengths, int[] eventTypes, int[] sampleNodes,
                                         int[] coalescentChild1, int[] coalescentChild2, int[] coalescentParent) {
        if (intervalLengths == null || eventTypes == null || sampleNodes == null ||
                coalescentChild1 == null || coalescentChild2 == null || coalescentParent == null) {
            throw new IllegalArgumentException("schedule arrays must not be null");
        }
        int length = intervalLengths.length;
        if (eventTypes.length != length || sampleNodes.length != length ||
                coalescentChild1.length != length || coalescentChild2.length != length ||
                coalescentParent.length != length) {
            throw new IllegalArgumentException("schedule arrays have inconsistent lengths");
        }
        if (length == 0) {
            throw new IllegalArgumentException("schedule must contain at least one interval");
        }
        return length;
    }

    private static void validateNode(int node, int nodeCount, String label) {
        if (node < 0 || node >= nodeCount) {
            throw new IllegalArgumentException(label + " node out of range: " + node +
                    " (nodeCount=" + nodeCount + ")");
        }
    }

    private static void requireNoNode(int node, String label, int event) {
        if (node != NO_NODE) {
            throw new IllegalArgumentException(label + " must be NO_NODE at event " + event + ": " + node);
        }
    }

}
