package dr.evomodel.continuous.ou.blockdiagonal;

import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;

/**
 * Stationary-mean pullback for orthogonal block OU transitions.
 *
 * <p>The transition-matrix adjoint (dLogL_dF -&gt; B/R gradient) is handled by the
 * inherited {@link BlockDiagonalSelectionMatrixParameterization} implementation, not
 * here: it consumes only {@code rMatrix}/{@code rinvMatrix} as dense arrays, and
 * {@code rinvMatrix} is already exactly {@code R^T} for an orthogonal rotation (see
 * {@code AbstractBlockDiagonalTwoByTwoMatrixParameter.fillRAndRinv}), so a separate
 * orthogonal-specific sandwich computation would be numerically identical. Only the
 * mean-gradient path here is genuinely orthogonal-specific, since it composes with
 * {@link OrthogonalMatrixProvider}'s vector-apply operations rather than a dense
 * rotation multiply.
 */
final class OrthogonalBlockSelectionAdjoint {

    private final int[] blockStarts;
    private final int[] blockSizes;
    private final OrthogonalMatrixProvider orthogonalRotation;

    OrthogonalBlockSelectionAdjoint(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                    final OrthogonalMatrixProvider orthogonalRotation) {
        this.blockStarts = blockParameter.getBlockStarts();
        this.blockSizes = blockParameter.getBlockSizes();
        this.orthogonalRotation = orthogonalRotation;
    }

    void accumulateMeanGradientPreparedDBasis(final OrthogonalBlockPreparedBranchBasis prepared,
                                             final double[] dLogL_df,
                                             final double[] dBasisGradientAccumulator,
                                             final BlockDiagonalBranchGradientWorkspace workspace) {
        multiplyRowMajorVector(prepared.rinvMatrix.data, dLogL_df, workspace.tempVector1, prepared.dimension);
        BlockDiagonalMatrixOps.applyBlockDiagonalTranspose(
                prepared.expD, workspace.tempVector1, workspace.tempVector2,
                prepared.dimension, blockStarts, blockSizes);
        for (int i = 0; i < prepared.dimension; ++i) {
            dBasisGradientAccumulator[i] += workspace.tempVector1[i] - workspace.tempVector2[i];
        }
    }

    void finishMeanGradientFromDBasis(final double[] dBasisGradientAccumulator,
                                      final double[] gradientAccumulator,
                                      final BlockDiagonalBranchGradientWorkspace workspace) {
        if (gradientAccumulator.length == 1) {
            orthogonalRotation.applyOrthogonal(dBasisGradientAccumulator, workspace.tempVector1);
            double sum = 0.0;
            for (int i = 0; i < workspace.tempVector1.length; ++i) {
                sum += workspace.tempVector1[i];
            }
            gradientAccumulator[0] += sum;
            return;
        }
        if (gradientAccumulator.length != dBasisGradientAccumulator.length) {
            throw new IllegalArgumentException(
                    "Stationary-mean gradient length must be 1 or "
                            + dBasisGradientAccumulator.length + ", found " + gradientAccumulator.length);
        }
        orthogonalRotation.applyOrthogonal(dBasisGradientAccumulator, workspace.tempVector1);
        for (int i = 0; i < dBasisGradientAccumulator.length; ++i) {
            gradientAccumulator[i] += workspace.tempVector1[i];
        }
    }

    private static void multiplyRowMajorVector(final double[] matrix,
                                               final double[] vector,
                                               final double[] out,
                                               final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            double sum = 0.0;
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                sum += matrix[rowOffset + j] * vector[j];
            }
            out[i] = sum;
        }
    }
}
