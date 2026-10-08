/*
 * AugmentedBufferIndexHelper.java
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
 * Double-buffered indices for the original nodes of a tree and the augmented (degree-2) nodes of an
 * AugmentedNodeRegistry, laid out so that allocating more augmented nodes only appends buffers.
 * <p>
 * Ids below min have one buffer each, at the id. Original ids min ... nodeCount - 1 have two copies, at id and
 * id + nodeCount - min, as in BufferIndexHelper(nodeCount, min). Augmented id nodeCount + j has its two copies at
 * blockBase + stride * j and blockBase + stride * j + 1; a stride above 2 leaves the other slots of the block to the
 * caller (getSlotIndex). Buffers below blockBase not used by the original ids are also the caller's.
 * <p>
 * The offsets of every augmented id exist from construction, so no index changes when more augmented nodes are
 * allocated; only getBufferCount() depends on the number allocated.
 *
 * @author Marc A Suchard
 */
public final class AugmentedBufferIndexHelper extends BufferIndexHelper {

    private final int min;
    private final int nodeCount;
    private final int blockBase;
    private final int stride;

    private final int[] flipSum; // the offset of copy 0 plus the offset of copy 1, by id - min
    private int[] offsets;       // index = id + offsets[id - min]
    private int[] storedOffsets;
    private final boolean[] flipped;

    private int allocated;

    /**
     * @param nodeCount    number of original nodes
     * @param min          ids below min have a single buffer
     * @param maxAugmented the most augmented nodes that can ever be allocated; initially all are
     * @param blockBase    first index of the block of augmented nodes
     * @param stride       buffers per augmented node, at least 2
     */
    public AugmentedBufferIndexHelper(int nodeCount, int min, int maxAugmented, int blockBase, int stride) {
        super(min, min); // empty arrays; every public method is overridden
        assert stride >= 2 && blockBase >= 2 * nodeCount - min : "augmented block overlaps the original nodes";

        this.min = min;
        this.nodeCount = nodeCount;
        this.blockBase = blockBase;
        this.stride = stride;
        this.allocated = maxAugmented;

        final int n = nodeCount - min + maxAugmented;
        offsets = new int[n];
        flipSum = new int[n];
        flipped = new boolean[n];
        for (int k = 0; k < n; ++k) {
            final int id = k + min;
            if (id < nodeCount) {
                flipSum[k] = nodeCount - min;
            } else {
                final int offset = blockBase + stride * (id - nodeCount) - id;
                offsets[k] = offset;
                flipSum[k] = 2 * offset + 1;
            }
        }
        storedOffsets = offsets.clone();
    }

    @Override
    protected int computeOffset(int bufferSetNumber) {
        return 0; // runs inside the superclass constructor
    }

    /**
     * @return the number of buffers for the augmented nodes allocated now
     */
    @Override
    public int getBufferCount() {
        return getBufferCount(allocated);
    }

    /**
     * @return the number of buffers when a augmented nodes are allocated
     */
    public int getBufferCount(int a) {
        return blockBase + stride * a;
    }

    @Override
    public int getOffsetIndex(int i) {
        final int index = (i < min) ? i : i + offsets[i - min];
        assert index < getBufferCount();
        return index;
    }

    @Override
    public void flipOffset(int i) {
        final int k = i - min;
        if (!flipped[k]) { // only flip once before reject / accept
            offsets[k] = flipSum[k] - offsets[k];
            flipped[k] = true;
        }
    }

    @Override
    public boolean isSafeUpdate(int i) {
        return offsets[i - min] != storedOffsets[i - min];
    }

    @Override
    public void storeState() {
        Arrays.fill(flipped, false);
        System.arraycopy(offsets, 0, storedOffsets, 0, offsets.length);
    }

    @Override
    public void restoreState() {
        final int[] tmp = storedOffsets;
        storedOffsets = offsets;
        offsets = tmp;
        Arrays.fill(flipped, false);
    }

    /**
     * @param id   an augmented id
     * @param slot 0 ... stride - 1; slots 0 and 1 are the two copies
     * @return the index of a buffer in the block of the augmented node
     */
    public int getSlotIndex(int id, int slot) {
        return blockBase + stride * (id - nodeCount) + slot;
    }

    /**
     * Sets the number of augmented nodes allocated, which getBufferCount() reports. Indices do not change.
     */
    public void setAllocated(int a) {
        assert a >= 0 && a <= offsets.length - (nodeCount - min);
        allocated = a;
    }
}
