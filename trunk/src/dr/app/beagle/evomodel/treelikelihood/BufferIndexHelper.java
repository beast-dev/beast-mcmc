/*
 * BufferIndexHelper.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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
 */

package dr.app.beagle.evomodel.treelikelihood;

import java.io.Serializable;

public class BufferIndexHelper implements Serializable {
	
    /**
     * @param maxIndexValue the number of possible input values for the index
     * @param minIndexValue the minimum index value to have the mirrored buffers
     */
    public BufferIndexHelper(int maxIndexValue, int minIndexValue) {
        this.maxIndexValue = maxIndexValue;
        this.minIndexValue = minIndexValue;

        offsetCount = maxIndexValue - minIndexValue;
        indexOffsets = new int[offsetCount];
        storedIndexOffsets = new int[offsetCount];
    }

    public int getBufferCount() {
        return 2 * offsetCount + minIndexValue;
    }

    public void flipOffset(int i) {
        if (i >= minIndexValue) {
            indexOffsets[i - minIndexValue] = offsetCount - indexOffsets[i - minIndexValue];
        } // else do nothing
    }

    public int getOffsetIndex(int i) {
        if (i < minIndexValue) {
            return i;
        }
        return indexOffsets[i - minIndexValue] + i;
    }

    public void getIndices(int[] outIndices) {
        for (int i = 0; i < maxIndexValue; i++) {
            outIndices[i] = getOffsetIndex(i);
        }
    }

    public void storeState() {
        System.arraycopy(indexOffsets, 0, storedIndexOffsets, 0, indexOffsets.length);

    }

    public void restoreState() {
        int[] tmp = storedIndexOffsets;
        storedIndexOffsets = indexOffsets;
        indexOffsets = tmp;
    }

    private final int maxIndexValue;
    private final int minIndexValue;
    private final int offsetCount;

    private int[] indexOffsets;
    private int[] storedIndexOffsets;

}//END: class
