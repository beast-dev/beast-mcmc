/*
 * AugmentedNodeRegistry.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodel.treedatalikelihood;

import java.util.Arrays;

/**
 * Persistent bookkeeping for the degree-2 nodes that an epoch model adds along the branches of a tree.
 * <p>
 * The nodes of the tree are numbered 0 ... nodeCount - 1. Augmented nodes take the numbers
 * nodeCount ... nodeCount + capacity - 1, so that they can index partial and matrix buffers beside the original nodes.
 * A branch (identified by its lower node) that crosses boundaries f, f + 1, ..., f + k - 1 of the epoch transition
 * times owns a chain of k augmented nodes, ordered from the bottom (boundary f) to the top (boundary f + k - 1).
 * <p>
 * Chains are allocated from a free list and persist between likelihood evaluations, so a branch that has not
 * changed keeps its buffers and its partials remain valid. All state is stored and restored together with the
 * buffer helpers of the delegate (copy-on-write, so proposals that do not change any chain cost nothing).
 * <p>
 * The matrix owned by a branch id (an original node or an augmented node) is the transition matrix of the segment
 * immediately above it. Its epoch is
 * <ul>
 * <li>f for the original node, the epoch that contains the height of the node, and</li>
 * <li>f + j + 1 for the augmented node at position j of the chain.</li>
 * </ul>
 * Epoch e lies between boundaries e - 1 and e (epoch 0 starts at the most recent time).
 *
 * @author Marc A Suchard
 */
public final class AugmentedNodeRegistry {

    private static final int NONE = -1;

    private final int nodeCount;
    private final int boundaryCount;
    private final int capacity;

    private final State current;
    private final State stored;
    private boolean storePending = true; // stored is not yet a snapshot of current

    private int highWaterMark; // like a capacity, never stored or restored

    /**
     * @param nodeCount     number of nodes in the tree
     * @param boundaryCount number of epoch transition times
     * @param capacity      maximum number of augmented nodes that can be in use at once; at most one lineage
     *                      per taxon can cross each boundary, so taxonCount * boundaryCount suffices when chains
     *                      are released before others are taken (as EpochLikelihoodTraversal does)
     */
    public AugmentedNodeRegistry(int nodeCount, int boundaryCount, int capacity) {
        if (nodeCount < 1 || boundaryCount < 0 || capacity < 0) {
            throw new IllegalArgumentException("Invalid registry dimensions");
        }
        this.nodeCount = nodeCount;
        this.boundaryCount = boundaryCount;
        this.capacity = capacity;

        current = new State(nodeCount, boundaryCount, capacity);
        stored = new State(nodeCount, boundaryCount, capacity);
    }

    public int getNodeCount() {
        return nodeCount;
    }

    public int getBoundaryCount() {
        return boundaryCount;
    }

    public int getCapacity() {
        return capacity;
    }

    /**
     * @return the number of buffers to allocate: original plus augmented nodes
     */
    public int getTotalNodeCount() {
        return nodeCount + capacity;
    }

    public int getFreeCount() {
        return current.freeTop;
    }

    /**
     * Ids are taken lowest first and released ids are reused before new ones, so augmented nodes
     * nodeCount ... nodeCount + getHighWaterMark() - 1 are the only ones ever used.
     *
     * @return the largest number of augmented nodes in use at once so far, including in states that were restored
     */
    public int getHighWaterMark() {
        return highWaterMark;
    }

    public boolean isAugmented(int id) {
        return id >= nodeCount;
    }

    /**
     * @return the number of augmented nodes on the branch above node
     */
    public int getChainLength(int branch) {
        return current.chainCount[branch];
    }

    /**
     * @return the highest augmented node on the branch above node, or node itself if there is none
     */
    public int getChainTop(int branch) {
        final int tail = current.chainTail[branch];
        return tail == NONE ? branch : tail;
    }

    /**
     * @return the index of the first (lowest) boundary crossed by the branch above node
     */
    public int getFirstBoundary(int branch) {
        return current.chainFirst[branch];
    }

    /**
     * Copies the augmented nodes on the branch above node from bottom to top.
     *
     * @param out array with room for at least getChainLength(branch) entries
     * @return the number of augmented nodes
     */
    public int copyChain(int branch, int[] out) {
        int count = current.chainCount[branch];
        int id = current.chainHead[branch];
        for (int j = 0; j < count; ++j) {
            out[j] = id;
            id = current.next[id - nodeCount];
        }
        return count;
    }

    /**
     * @return the epoch of the transition matrix owned by an original or augmented node
     */
    public int getMatrixEpoch(int id) {
        return current.epochOfMatrix[id];
    }

    /**
     * @return the epoch transition time that the chains were assigned with
     */
    public double getBoundary(int index) {
        return current.boundaries[index];
    }

    /**
     * @return true if the chains were assigned with these epoch transition times
     */
    public boolean hasBoundaries(double[] boundaries) {
        for (int i = 0; i < boundaryCount; ++i) {
            if (boundaries[i] != current.boundaries[i]) {
                return false;
            }
        }
        return true;
    }

    /**
     * Sets the epoch transition times.
     *
     * @return true if the times changed, in which case every chain must be assigned again
     */
    public boolean setBoundaries(double[] boundaries) {
        if (boundaries.length != boundaryCount) {
            throw new IllegalArgumentException("Expected " + boundaryCount + " epoch transition times but found " +
                    boundaries.length);
        }

        boolean changed = false;
        for (int i = 0; i < boundaryCount; ++i) {
            if (boundaries[i] != current.boundaries[i]) { // NaN is never equal, so the first call always changes
                changed = true;
                break;
            }
        }

        if (changed) {
            beforeMutation();
            System.arraycopy(boundaries, 0, current.boundaries, 0, boundaryCount);
        }
        return changed;
    }

    /**
     * Sets the chain on the branch above node to k = count augmented nodes at boundaries firstBoundary, ...,
     * firstBoundary + count - 1. Nodes already in the chain are kept (from the bottom) and only the difference
     * is released to or taken from the free list.
     *
     * @throws IllegalStateException if there are too few free augmented nodes
     */
    public void assign(int branch, int firstBoundary, int count) {

        if (branch < 0 || branch >= nodeCount) {
            throw new IllegalArgumentException("Invalid branch " + branch);
        }
        if (firstBoundary < 0 || count < 0 || firstBoundary + count > boundaryCount) {
            throw new IllegalArgumentException("Invalid chain boundaries [" + firstBoundary + ", " +
                    (firstBoundary + count) + ") of " + boundaryCount);
        }

        final State s = current;
        final int old = s.chainCount[branch];

        if (old == count && s.chainFirst[branch] == firstBoundary) {
            return;
        }

        // not reached when chains are released before others are taken; BEAGLE buffers follow getHighWaterMark()
        if (count - old > s.freeTop) {
            throw new IllegalStateException("Out of augmented nodes: capacity " + capacity +
                    " is too small for the epoch transition times and tree");
        }

        beforeMutation();

        final int keep = Math.min(old, count);

        // Walk to the end of the part of the chain that is kept
        int last = NONE;
        int id = s.chainHead[branch];
        for (int j = 0; j < keep; ++j) {
            last = id;
            id = s.next[id - nodeCount];
        }

        // Release the rest of the old chain
        for (int j = keep; j < old; ++j) {
            final int following = s.next[id - nodeCount];
            s.freeStack[s.freeTop++] = id;
            id = following;
        }

        // Take what is missing
        for (int j = keep; j < count; ++j) {
            final int taken = s.freeStack[--s.freeTop];
            if (taken - nodeCount >= highWaterMark) {
                highWaterMark = taken - nodeCount + 1;
            }
            s.next[taken - nodeCount] = NONE;
            if (last == NONE) {
                s.chainHead[branch] = taken;
            } else {
                s.next[last - nodeCount] = taken;
            }
            last = taken;
        }

        if (count == 0) {
            s.chainHead[branch] = NONE;
        } else {
            s.next[last - nodeCount] = NONE;
        }
        s.chainTail[branch] = last; // NONE when the chain is empty

        s.chainCount[branch] = count;
        s.chainFirst[branch] = firstBoundary;

        s.epochOfMatrix[branch] = firstBoundary;
        id = s.chainHead[branch];
        for (int j = 0; j < count; ++j) {
            s.epochOfMatrix[id] = firstBoundary + j + 1;
            id = s.next[id - nodeCount];
        }
    }

    public void storeState() {
        storePending = true;
    }

    public void restoreState() {
        if (!storePending) { // otherwise nothing changed since storeState()
            current.copyFrom(stored);
            storePending = true;
        }
    }

    private void beforeMutation() {
        if (storePending) {
            stored.copyFrom(current);
            storePending = false;
        }
    }

    private static final class State {

        final int[] chainHead;
        final int[] chainTail;
        final int[] chainCount;
        final int[] chainFirst;
        final int[] next; // next node up the chain, indexed by (id - nodeCount)
        final int[] epochOfMatrix;
        final int[] freeStack;
        int freeTop;
        final double[] boundaries;

        State(int nodeCount, int boundaryCount, int capacity) {
            chainHead = new int[nodeCount];
            chainTail = new int[nodeCount];
            chainCount = new int[nodeCount];
            chainFirst = new int[nodeCount];
            next = new int[capacity];
            epochOfMatrix = new int[nodeCount + capacity];
            freeStack = new int[capacity];
            boundaries = new double[boundaryCount];

            Arrays.fill(chainHead, NONE);
            Arrays.fill(chainTail, NONE);
            Arrays.fill(next, NONE);
            Arrays.fill(boundaries, Double.NaN);

            for (int i = 0; i < capacity; ++i) { // lowest number is taken first
                freeStack[i] = nodeCount + capacity - 1 - i;
            }
            freeTop = capacity;
        }

        void copyFrom(State source) {
            System.arraycopy(source.chainHead, 0, chainHead, 0, chainHead.length);
            System.arraycopy(source.chainTail, 0, chainTail, 0, chainTail.length);
            System.arraycopy(source.chainCount, 0, chainCount, 0, chainCount.length);
            System.arraycopy(source.chainFirst, 0, chainFirst, 0, chainFirst.length);
            System.arraycopy(source.next, 0, next, 0, next.length);
            System.arraycopy(source.epochOfMatrix, 0, epochOfMatrix, 0, epochOfMatrix.length);
            System.arraycopy(source.freeStack, 0, freeStack, 0, freeStack.length);
            System.arraycopy(source.boundaries, 0, boundaries, 0, boundaries.length);
            freeTop = source.freeTop;
        }
    }
}
