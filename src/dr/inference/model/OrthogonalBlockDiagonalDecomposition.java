package dr.inference.model;

/**
 * Block-diagonal decomposition for the special case
 *
 *     A = R D R^T,
 *
 * where R is orthogonal.
 */
public final class OrthogonalBlockDiagonalDecomposition extends BlockDiagonalDecomposition {

    public OrthogonalBlockDiagonalDecomposition(final int dimension,
                                                final int[] blockStarts,
                                                final int[] blockSizes) {
        super(dimension, blockStarts, blockSizes);
    }

    public OrthogonalBlockDiagonalDecomposition(final double[] r,
                                                final double[] rTranspose,
                                                final double[] blockDiagonal,
                                                final int[] blockStarts,
                                                final int[] blockSizes) {
        super(r, rTranspose, blockDiagonal, blockStarts, blockSizes);
    }

    @Override
    public OrthogonalBlockDiagonalDecomposition copy() {
        return new OrthogonalBlockDiagonalDecomposition(
                r.clone(),
                rInverse.clone(),
                blockDiagonal.clone(),
                blockStarts,
                blockSizes);
    }

    @Override
    public boolean isOrthogonal() {
        return true;
    }

    public double[] getRTranspose() {
        return getRInverse();
    }

    public void updateInverseFromRTranspose() {
        transpose(r, rInverse, dimension);
    }
}
