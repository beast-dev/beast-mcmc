/*
 * StructuredCoalescentActiveLineages.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent;

import java.util.Arrays;

/**
 * Primitive active-lineage set for structured-coalescent schedules.
 *
 * This class tracks only biological node ids and their active/inactive state.
 * It deliberately contains no BASTA buffer indices, BEAGLE matrix ids, MASCOT
 * probability slots, or likelihood state.
 */
public final class StructuredCoalescentActiveLineages {

    private final int nodeCount;
    private final int[] activeNodes;
    private final int[] activeIndexByNode;
    private final boolean preserveOrder;
    private int activeCount;

    public StructuredCoalescentActiveLineages(int nodeCount) {
        this(nodeCount, true);
    }

    public static StructuredCoalescentActiveLineages unordered(int nodeCount) {
        return new StructuredCoalescentActiveLineages(nodeCount, false);
    }

    private StructuredCoalescentActiveLineages(int nodeCount, boolean preserveOrder) {
        if (nodeCount <= 0) {
            throw new IllegalArgumentException("nodeCount must be positive");
        }
        this.nodeCount = nodeCount;
        this.activeNodes = new int[nodeCount];
        this.activeIndexByNode = new int[nodeCount];
        this.preserveOrder = preserveOrder;
        Arrays.fill(activeIndexByNode, -1);
    }

    public void reset(int initialSampleNode) {
        clear();
        addSample(initialSampleNode);
    }

    public int getActiveCount() {
        return activeCount;
    }

    public int getNodeCount() {
        return nodeCount;
    }

    public int getActiveNode(int index) {
        if (index < 0 || index >= activeCount) {
            throw new IllegalArgumentException("active lineage index out of range: " + index);
        }
        return activeNodes[index];
    }

    public boolean isActive(int node) {
        validateNode(node, "node");
        return activeIndexByNode[node] >= 0;
    }

    public void addSample(int node) {
        validateNode(node, "sample");
        if (activeIndexByNode[node] >= 0) {
            throw new IllegalArgumentException("sample lineage is already active: " + node);
        }
        activeNodes[activeCount] = node;
        activeIndexByNode[node] = activeCount;
        activeCount++;
    }

    public void replaceChildrenWithParent(int leftChild, int rightChild, int parent) {
        validateNode(leftChild, "left child");
        validateNode(rightChild, "right child");
        validateNode(parent, "parent");
        if (leftChild == rightChild) {
            throw new IllegalArgumentException("coalescent children must be distinct");
        }
        if (activeIndexByNode[leftChild] < 0 || activeIndexByNode[rightChild] < 0) {
            throw new IllegalArgumentException("coalescent event uses inactive children: " +
                    leftChild + ", " + rightChild);
        }
        if (activeIndexByNode[parent] >= 0) {
            throw new IllegalArgumentException("coalescent parent is already active: " + parent);
        }

        if (preserveOrder) {
            removePreservingOrder(leftChild);
            removePreservingOrder(rightChild);
        } else {
            removeSwapWithLast(leftChild);
            removeSwapWithLast(rightChild);
        }
        addSample(parent);
    }

    public void requireSingleActiveLineage() {
        if (activeCount != 1) {
            throw new IllegalArgumentException("structured-coalescent schedule ends with " + activeCount +
                    " active lineages");
        }
    }

    private void clear() {
        for (int i = 0; i < activeCount; i++) {
            activeIndexByNode[activeNodes[i]] = -1;
        }
        activeCount = 0;
    }

    private void removePreservingOrder(int node) {
        int index = activeIndexByNode[node];
        if (index < 0) {
            throw new IllegalArgumentException("lineage is not active: " + node);
        }
        for (int i = index + 1; i < activeCount; i++) {
            int movedNode = activeNodes[i];
            activeNodes[i - 1] = movedNode;
            activeIndexByNode[movedNode] = i - 1;
        }
        activeCount--;
        activeIndexByNode[node] = -1;
    }

    private void removeSwapWithLast(int node) {
        int index = activeIndexByNode[node];
        if (index < 0) {
            throw new IllegalArgumentException("lineage is not active: " + node);
        }
        int lastIndex = activeCount - 1;
        int lastNode = activeNodes[lastIndex];
        if (index != lastIndex) {
            activeNodes[index] = lastNode;
            activeIndexByNode[lastNode] = index;
        }
        activeCount--;
        activeIndexByNode[node] = -1;
    }

    private void validateNode(int node, String label) {
        if (node < 0 || node >= nodeCount) {
            throw new IllegalArgumentException(label + " node out of range: " + node +
                    " (nodeCount=" + nodeCount + ")");
        }
    }
}
