package dr.app.beagle.evomodel.treelikelihood;

public class BufferIndexHelper {
	
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
