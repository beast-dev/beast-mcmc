package dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.continuous.ou.blockdiagonal.math.BlockDiagonalExpSolver;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;

import java.util.Arrays;

/**
 * Immutable block topology for a block-diagonal OU selection matrix: the
 * dimension and block boundaries, plus the counts derived from them (number
 * of 2x2 blocks, compressed-D dimension, tridiagonal-D dimension) that were
 * previously recomputed independently in several classes across this
 * package.
 */
public final class BlockDiagonalStructure {

    private final int dimension;
    private final int[] blockStarts;
    private final int[] blockSizes;
    private final int numTwoByTwoBlocks;
    private final int compressedDimension;
    private final int tridiagonalDimension;
    private final BlockDiagonalExpSolver.BlockStructure mathStructure;

    public BlockDiagonalStructure(final int dimension,
                                  final int[] blockStarts,
                                  final int[] blockSizes) {
        this.mathStructure = new BlockDiagonalExpSolver.BlockStructure(dimension, blockStarts, blockSizes);
        this.dimension = dimension;
        this.blockStarts = blockStarts.clone();
        this.blockSizes = blockSizes.clone();
        this.numTwoByTwoBlocks = numTwoByTwoBlocks(this.blockSizes);
        this.compressedDimension = dimension + 2 * numTwoByTwoBlocks;
        this.tridiagonalDimension = 3 * dimension - 2;
    }

    public static BlockDiagonalStructure from(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter) {
        return new BlockDiagonalStructure(
                blockParameter.getRowDimension(),
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes());
    }

    public static int numTwoByTwoBlocks(final int[] blockSizes) {
        int count = 0;
        for (final int size : blockSizes) {
            if (size == 2) {
                count++;
            }
        }
        return count;
    }

    public int getDimension() {
        return dimension;
    }

    public int[] getBlockStarts() {
        return blockStarts.clone();
    }

    public int[] getBlockSizes() {
        return blockSizes.clone();
    }

    public int getNumTwoByTwoBlocks() {
        return numTwoByTwoBlocks;
    }

    public int getCompressedDimension() {
        return compressedDimension;
    }

    public int getTridiagonalDimension() {
        return tridiagonalDimension;
    }

    /** The narrower topology type consumed directly by {@link BlockDiagonalExpSolver}. */
    public BlockDiagonalExpSolver.BlockStructure getMathStructure() {
        return mathStructure;
    }

    /**
     * Fail-fast structural-compatibility check: throws if {@code other} does not
     * describe the same block topology as this structure. Intended for callers
     * that combine a builder/helper (constructed once, own fixed topology) with
     * a separately-supplied snapshot/parameter object on every call.
     */
    public void checkMatches(final BlockDiagonalStructure other) {
        if (dimension != other.dimension
                || numTwoByTwoBlocks != other.numTwoByTwoBlocks
                || !Arrays.equals(blockStarts, other.blockStarts)
                || !Arrays.equals(blockSizes, other.blockSizes)) {
            throw new IllegalArgumentException(
                    "Block structure mismatch: expected dimension=" + dimension
                            + " blockStarts=" + Arrays.toString(blockStarts)
                            + " blockSizes=" + Arrays.toString(blockSizes)
                            + " but got dimension=" + other.dimension
                            + " blockStarts=" + Arrays.toString(other.blockStarts)
                            + " blockSizes=" + Arrays.toString(other.blockSizes));
        }
    }
}
