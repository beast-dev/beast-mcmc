package dr.evomodel.continuous.ou.blockdiagonal;

import dr.inference.model.BlockDiagonalDecomposition;

/**
 * Reusable branch-local basis prepared for one orthogonal block-diagonal OU transition.
 */
public final class OrthogonalBlockPreparedBranchBasis extends BlockDiagonalPreparedBranchBasis {

    OrthogonalBlockPreparedBranchBasis(final int dimension,
                                       final int compressedBlockDimension,
                                       final BlockDiagonalDecomposition decomposition) {
        super(dimension, compressedBlockDimension, decomposition);
    }
}
