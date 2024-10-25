/*
 * BufferIndexHelper.java
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

import java.io.Serializable;
import java.util.*;

/**
 * BufferIndexHelper - helper for double buffering of intermediate computation at nodes.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 */
public class BufferIndexHelper implements Serializable {

    /**
     * @param maxIndexValue the number of possible input values for the index
     * @param minIndexValue the minimum index value to have the mirrored buffers
     */
    public BufferIndexHelper(int maxIndexValue, int minIndexValue) {
        this(maxIndexValue, minIndexValue, 0);
    }

    /**
     * @param maxIndexValue the number of possible input values for the index
     * @param minIndexValue the minimum index value to have the mirrored buffers
     * @param bufferSetNumber provides a total offset of bufferSetNumber * bufferCount
     */
    public BufferIndexHelper(int maxIndexValue, int minIndexValue, int bufferSetNumber) {
        this.minIndexValue = minIndexValue;

        doubleBufferCount = maxIndexValue - minIndexValue;
        indexOffsets = new int[doubleBufferCount];
        storedIndexOffsets = new int[doubleBufferCount];
        indexOffsetsFlipped = new boolean[doubleBufferCount];

        this.constantOffset = computeOffset(bufferSetNumber);
    }

    protected int computeOffset(int bufferSetNumber) { return  bufferSetNumber * getBufferCount(); }

    public int getBufferCount() {
        return 2 * doubleBufferCount + minIndexValue;
    }

    public void flipOffset(int i) {
        assert(i >= minIndexValue) : "shouldn't be trying to flip the first 'static' indices";

        if (!indexOffsetsFlipped[i - minIndexValue]){ // only flip once before reject / accept
            indexOffsets[i - minIndexValue] = doubleBufferCount - indexOffsets[i - minIndexValue];
            indexOffsetsFlipped[i - minIndexValue] = true;
        }
    }

    public int getOffsetIndex(int i) {
        if (i < minIndexValue) {
            return i + constantOffset;
        }
        return indexOffsets[i - minIndexValue] + i + constantOffset;
    }

    private int getStoredOffsetIndex(int i) {
        assert (i >= minIndexValue);
        return storedIndexOffsets[i - minIndexValue] + i + constantOffset;
    }

    public boolean isSafeUpdate(int i) {
        return getStoredOffsetIndex(i) != getOffsetIndex(i);
    }

    public void storeState() {
        Arrays.fill(indexOffsetsFlipped, false);
        System.arraycopy(indexOffsets, 0, storedIndexOffsets, 0, indexOffsets.length);
    }

    public void restoreState() {
        int[] tmp = storedIndexOffsets;
        storedIndexOffsets = indexOffsets;
        indexOffsets = tmp;
        Arrays.fill(indexOffsetsFlipped, false);
    }

    private final int minIndexValue;
    private final int constantOffset;
    private final int doubleBufferCount;
    private final boolean[] indexOffsetsFlipped;

    private int[] indexOffsets;
    private int[] storedIndexOffsets;

}//END: class
