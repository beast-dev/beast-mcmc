/*
 * MascotBranchCells.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import java.util.Arrays;

/**
 * Time cells (one per lineage per RK4 step, possibly split) recorded in tape order, then grouped by branch and
 * snapped to the tree: a branch's lowest cell starts at its node's height, each cell ends where the next one starts,
 * and the highest ends at the parent's height. The lengths the forward pass used are kept for the cells' dynamics;
 * the snapped boundaries are used only to place jumps. Each cell carries {@code payloadWidth} doubles chosen by the
 * sampler.
 */
final class MascotBranchCells {

    private static final double SNAP_TOLERANCE = 1.0e-9;

    private final int payloadWidth;

    private int rawCount;
    private int[] rawBranch = new int[0];
    private int[] rawEpoch = new int[0];
    private double[] rawStart = new double[0];
    private double[] rawLength = new double[0];
    private double[] rawClock = new double[0];
    private double[] rawPayload = new double[0];

    private int count;
    private int[] branchFirst = new int[0];
    private int[] branchCount = new int[0];
    private int[] epoch = new int[0];
    private double[] start = new double[0];
    private double[] end = new double[0];
    private double[] length = new double[0];
    private double[] clock = new double[0];
    private double[] payload = new double[0];

    MascotBranchCells(int payloadWidth) {
        this.payloadWidth = payloadWidth;
    }

    void reset() {
        rawCount = 0;
        count = 0;
    }

    /** Appends a cell; {@code payloadWidth} doubles are copied from {@code source} at {@code offset}. */
    void add(int branch, int cellEpoch, double cellStart, double cellLength, double clockRate,
             double[] source, int offset) {
        if (rawCount == rawBranch.length) {
            int capacity = Math.max(64, 2 * rawCount);
            rawBranch = Arrays.copyOf(rawBranch, capacity);
            rawEpoch = Arrays.copyOf(rawEpoch, capacity);
            rawStart = Arrays.copyOf(rawStart, capacity);
            rawLength = Arrays.copyOf(rawLength, capacity);
            rawClock = Arrays.copyOf(rawClock, capacity);
            rawPayload = Arrays.copyOf(rawPayload, capacity * payloadWidth);
        }
        rawBranch[rawCount] = branch;
        rawEpoch[rawCount] = cellEpoch;
        rawStart[rawCount] = cellStart;
        rawLength[rawCount] = cellLength;
        rawClock[rawCount] = clockRate;
        System.arraycopy(source, offset, rawPayload, rawCount * payloadWidth, payloadWidth);
        rawCount++;
    }

    /** Stable counting sort by branch, then snapping and consistency checks against the tree. */
    void group(int nodeCount, double[] nodeHeights, int[] parentOf, int rootNode) {
        count = rawCount;
        branchFirst = ensureInt(branchFirst, nodeCount);
        branchCount = ensureInt(branchCount, nodeCount);
        Arrays.fill(branchCount, 0, nodeCount, 0);
        for (int c = 0; c < rawCount; c++) {
            branchCount[rawBranch[c]]++;
        }
        int next = 0;
        for (int b = 0; b < nodeCount; b++) {
            branchFirst[b] = next;
            next += branchCount[b];
        }

        epoch = ensureInt(epoch, count);
        start = ensure(start, count);
        end = ensure(end, count);
        length = ensure(length, count);
        clock = ensure(clock, count);
        payload = ensure(payload, count * payloadWidth);
        int[] fill = new int[nodeCount];
        for (int c = 0; c < rawCount; c++) {
            int b = rawBranch[c];
            int target = branchFirst[b] + fill[b]++;
            epoch[target] = rawEpoch[c];
            start[target] = rawStart[c];
            length[target] = rawLength[c];
            clock[target] = rawClock[c];
            System.arraycopy(rawPayload, c * payloadWidth, payload, target * payloadWidth, payloadWidth);
        }

        for (int b = 0; b < nodeCount; b++) {
            if (b == rootNode) {
                if (branchCount[b] != 0) {
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
            if (branchCount[b] == 0) {
                requireClose(bottom, top, "zero-cell branch above node " + b);
                continue;
            }
            int first = branchFirst[b];
            int last = first + branchCount[b] - 1;
            requireClose(start[first], bottom, "lowest cell of the branch above node " + b);
            requireClose(start[last] + length[last], top, "highest cell of the branch above node " + b);
            start[first] = bottom;
            for (int c = first; c < last; c++) {
                requireClose(start[c] + length[c], start[c + 1], "consecutive cells of the branch above node " + b);
                end[c] = start[c + 1];
            }
            end[last] = top;
            for (int c = first; c <= last; c++) {
                if (!(end[c] > start[c])) {
                    throw new IllegalStateException("empty cell in the branch above node " + b);
                }
            }
        }
    }

    int count() {
        return count;
    }

    int first(int branch) {
        return branchFirst[branch];
    }

    int countOf(int branch) {
        return branchCount[branch];
    }

    int epoch(int cell) {
        return epoch[cell];
    }

    double start(int cell) {
        return start[cell];
    }

    double end(int cell) {
        return end[cell];
    }

    /** Length the forward pass used for this cell. */
    double length(int cell) {
        return length[cell];
    }

    double clock(int cell) {
        return clock[cell];
    }

    double[] payload() {
        return payload;
    }

    int payloadOffset(int cell) {
        return cell * payloadWidth;
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
