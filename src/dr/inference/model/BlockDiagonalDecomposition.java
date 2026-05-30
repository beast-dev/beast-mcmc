package dr.inference.model;

import java.io.Serializable;
import java.util.Arrays;

/**
 * Data carrier for a block-diagonal similarity decomposition
 *
 *     A = R D R^{-1}.
 *
 * The middle matrix D is stored in the raw block layout used by the OU block
 * helpers: [diag(0..dim-1), upper(0..dim-2), lower(0..dim-2)].
 */
public class BlockDiagonalDecomposition implements Serializable {

    public BlockDiagonalDecomposition(final int dimension,
                                      final int[] blockStarts,
                                      final int[] blockSizes) {
        this(new double[matrixSize(dimension)],
                new double[matrixSize(dimension)],
                new double[blockDiagonalSize(dimension)],
                blockStarts,
                blockSizes);
    }

    public BlockDiagonalDecomposition(final double[] r,
                                      final double[] rInverse,
                                      final double[] blockDiagonal,
                                      final int[] blockStarts,
                                      final int[] blockSizes) {
        this.dimension = inferDimension(r);
        validateSquareArray(r, dimension, "R");
        validateSquareArray(rInverse, dimension, "R inverse");
        validateBlockDiagonalArray(blockDiagonal, dimension);
        validateBlockStructure(dimension, blockStarts, blockSizes);
        this.r = r;
        this.rInverse = rInverse;
        this.blockDiagonal = blockDiagonal;
        this.blockStarts = blockStarts.clone();
        this.blockSizes = blockSizes.clone();
    }

    public BlockDiagonalDecomposition copy() {
        return new BlockDiagonalDecomposition(
                r.clone(),
                rInverse.clone(),
                blockDiagonal.clone(),
                blockStarts,
                blockSizes);
    }

    public int getDimension() {
        return dimension;
    }

    public int getBlockDiagonalDimension() {
        return blockDiagonal.length;
    }

    public int getNumBlocks() {
        return blockStarts.length;
    }

    public boolean isOrthogonal() {
        return false;
    }

    public double[] getR() {
        return r;
    }

    public double[] getRInverse() {
        return rInverse;
    }

    public double[] getBlockDiagonal() {
        return blockDiagonal;
    }

    public int[] getBlockStarts() {
        return blockStarts;
    }

    public int[] getBlockSizes() {
        return blockSizes;
    }

    public void fillR(final double[] out) {
        checkSquareOutput(out, "R output");
        System.arraycopy(r, 0, out, 0, r.length);
    }

    public void fillRInverse(final double[] out) {
        checkSquareOutput(out, "R inverse output");
        System.arraycopy(rInverse, 0, out, 0, rInverse.length);
    }

    public void fillBlockDiagonal(final double[] out) {
        if (out == null || out.length != blockDiagonal.length) {
            throw new IllegalArgumentException(
                    "Block diagonal output must have length " + blockDiagonal.length);
        }
        System.arraycopy(blockDiagonal, 0, out, 0, blockDiagonal.length);
    }

    public void fillBlockDiagonalMatrix(final double[] outRowMajor) {
        checkSquareOutput(outRowMajor, "Block diagonal matrix output");
        Arrays.fill(outRowMajor, 0.0);

        final int upperOffset = dimension;
        final int lowerOffset = dimension + (dimension - 1);
        for (int b = 0; b < blockStarts.length; b++) {
            final int start = blockStarts[b];
            final int size = blockSizes[b];
            if (size == 1) {
                outRowMajor[start * dimension + start] = blockDiagonal[start];
            } else {
                outRowMajor[start * dimension + start] = blockDiagonal[start];
                outRowMajor[start * dimension + start + 1] = blockDiagonal[upperOffset + start];
                outRowMajor[(start + 1) * dimension + start] = blockDiagonal[lowerOffset + start];
                outRowMajor[(start + 1) * dimension + start + 1] = blockDiagonal[start + 1];
            }
        }
    }

    protected static void transpose(final double[] matrix,
                                    final double[] transpose,
                                    final int dimension) {
        validateSquareArray(matrix, dimension, "matrix");
        validateSquareArray(transpose, dimension, "transpose");
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                transpose[i * dimension + j] = matrix[j * dimension + i];
            }
        }
    }

    private void checkSquareOutput(final double[] out,
                                   final String name) {
        validateSquareArray(out, dimension, name);
    }

    private static int inferDimension(final double[] matrix) {
        if (matrix == null) {
            throw new IllegalArgumentException("R must not be null");
        }
        final int dimension = (int) Math.round(Math.sqrt(matrix.length));
        if (dimension * dimension != matrix.length) {
            throw new IllegalArgumentException("R must be a square row-major matrix");
        }
        return dimension;
    }

    private static void validateSquareArray(final double[] matrix,
                                            final int dimension,
                                            final String name) {
        if (matrix == null || matrix.length != dimension * dimension) {
            throw new IllegalArgumentException(
                    name + " must have length " + (dimension * dimension));
        }
    }

    private static void validateBlockDiagonalArray(final double[] blockDiagonal,
                                                   final int dimension) {
        final int expected = blockDiagonalSize(dimension);
        if (blockDiagonal == null || blockDiagonal.length != expected) {
            throw new IllegalArgumentException(
                    "Block diagonal data must have length " + expected);
        }
    }

    private static int matrixSize(final int dimension) {
        validateDimension(dimension);
        return dimension * dimension;
    }

    private static int blockDiagonalSize(final int dimension) {
        validateDimension(dimension);
        return 3 * dimension - 2;
    }

    private static void validateDimension(final int dimension) {
        if (dimension <= 0) {
            throw new IllegalArgumentException("dimension must be positive");
        }
    }

    private static void validateBlockStructure(final int dimension,
                                               final int[] blockStarts,
                                               final int[] blockSizes) {
        if (blockStarts == null || blockSizes == null) {
            throw new IllegalArgumentException("Block starts and sizes must not be null");
        }
        if (blockStarts.length != blockSizes.length) {
            throw new IllegalArgumentException("Block starts and sizes must have the same length");
        }
        int covered = 0;
        for (int b = 0; b < blockStarts.length; b++) {
            final int start = blockStarts[b];
            final int size = blockSizes[b];
            if (size != 1 && size != 2) {
                throw new IllegalArgumentException("Only block sizes 1 and 2 are supported");
            }
            if (start != covered) {
                throw new IllegalArgumentException("Block structure must cover the matrix contiguously");
            }
            covered += size;
        }
        if (covered != dimension) {
            throw new IllegalArgumentException("Block structure must cover dimension " + dimension);
        }
    }

    protected final int dimension;
    protected final double[] r;
    protected final double[] rInverse;
    protected final double[] blockDiagonal;
    protected final int[] blockStarts;
    protected final int[] blockSizes;
}
