package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.continuous.ou.MatrixExponentialUtils;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalFrechetHelper;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Pullback for orthogonal block OU transition matrix and offset adjoints.
 */
final class OrthogonalBlockSelectionAdjoint {

    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final OrthogonalMatrixProvider orthogonalRotation;
    private final BlockDiagonalFrechetHelper frechetHelper;
    private final OrthogonalBlockDenseFallbackPolicy denseFallbackPolicy;
    private final DenseMatrix64F upstreamF;
    private final DenseMatrix64F upstreamFD;
    private final DenseMatrix64F gradD;
    private final DenseMatrix64F gradR;
    private final DenseMatrix64F temp1;
    private final DenseMatrix64F temp2;
    private final DenseMatrix64F temp3;
    private final double[] scaledNegativeBlockDScratch;
    private final double[] denseAdjointScratch;
    private final double[] tempVector1;
    private final double[] tempVector2;

    OrthogonalBlockSelectionAdjoint(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                    final OrthogonalMatrixProvider orthogonalRotation,
                                    final BlockDiagonalFrechetHelper frechetHelper,
                                    final OrthogonalBlockDenseFallbackPolicy denseFallbackPolicy) {
        this.blockParameter = blockParameter;
        this.orthogonalRotation = orthogonalRotation;
        this.frechetHelper = frechetHelper;
        this.denseFallbackPolicy = denseFallbackPolicy;
        final int d = blockParameter.getRowDimension();
        this.upstreamF = new DenseMatrix64F(d, d);
        this.upstreamFD = new DenseMatrix64F(d, d);
        this.gradD = new DenseMatrix64F(d, d);
        this.gradR = new DenseMatrix64F(d, d);
        this.temp1 = new DenseMatrix64F(d, d);
        this.temp2 = new DenseMatrix64F(d, d);
        this.temp3 = new DenseMatrix64F(d, d);
        this.scaledNegativeBlockDScratch = new double[d * d];
        this.denseAdjointScratch = new double[d * d];
        this.tempVector1 = new double[d];
        this.tempVector2 = new double[d];
    }

    void accumulateCurrentFlat(final OrthogonalBlockBasisCache basis,
                               final double dt,
                               final double[] stationaryMean,
                               final double[] dLogL_dF,
                               final double[] dLogL_df,
                               final double[] compressedDAccumulator,
                               final double[] rotationAccumulator) {
        fillTotalUpstreamOnTransition(stationaryMean, dLogL_dF, dLogL_df, upstreamF);
        accumulateTransitionPullback(
                basis.rMatrix,
                basis.rtMatrix,
                basis.expD,
                basis.blockDParams,
                dt,
                upstreamF,
                upstreamFD,
                frechetHelper,
                temp1,
                temp2,
                temp3,
                gradD,
                gradR,
                scaledNegativeBlockDScratch,
                denseAdjointScratch,
                compressedDAccumulator,
                rotationAccumulator);
    }

    void accumulatePreparedFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                final double[] dLogL_dF,
                                final double[] dLogL_df,
                                final OrthogonalBlockBranchGradientWorkspace workspace,
                                final double[] compressedDAccumulator,
                                final double[] rotationAccumulator) {
        fillTotalUpstreamOnTransition(prepared.stationaryMean, dLogL_dF, dLogL_df, workspace.upstreamF);
        accumulateTransitionPullback(
                prepared.rMatrix,
                prepared.rtMatrix,
                prepared.expD,
                prepared.blockDParams,
                prepared.dt,
                workspace.upstreamF,
                workspace.upstreamFD,
                workspace.frechetHelper,
                workspace.temp1,
                workspace.temp2,
                workspace.temp3,
                workspace.gradD,
                workspace.gradR,
                workspace.scaledNegativeBlockDScratch,
                workspace.denseAdjointScratch,
                compressedDAccumulator,
                rotationAccumulator);
    }

    void accumulateMeanGradient(final OrthogonalBlockBasisCache basis,
                                final double[] dLogL_df,
                                final double[] gradientAccumulator) {
        orthogonalRotation.applyOrthogonalTranspose(dLogL_df, tempVector1);
        applyBlockDiagonalTransposeExp(basis.expD.data, tempVector1, tempVector2, basis.expD.numRows);
        orthogonalRotation.applyOrthogonal(tempVector2, tempVector1);
        for (int i = 0; i < basis.expD.numRows; ++i) {
            gradientAccumulator[i] += dLogL_df[i] - tempVector1[i];
        }
    }

    double accumulateScalarMeanGradient(final OrthogonalBlockBasisCache basis,
                                        final double[] dLogL_df) {
        orthogonalRotation.applyOrthogonalTranspose(dLogL_df, tempVector1);
        applyBlockDiagonalTransposeExp(basis.expD.data, tempVector1, tempVector2, basis.expD.numRows);
        orthogonalRotation.applyOrthogonal(tempVector2, tempVector1);
        double sum = 0.0;
        for (int i = 0; i < basis.expD.numRows; ++i) {
            sum += dLogL_df[i] - tempVector1[i];
        }
        return sum;
    }

    void accumulateMeanGradientPrepared(final OrthogonalBlockPreparedBranchBasis prepared,
                                        final double[] dLogL_df,
                                        final double[] gradientAccumulator,
                                        final OrthogonalBlockBranchGradientWorkspace workspace) {
        multiplyRowMajorVector(prepared.rtMatrix.data, dLogL_df, workspace.tempVector1, prepared.dimension);
        applyBlockDiagonalTransposeExp(
                prepared.expD.data, workspace.tempVector1, workspace.tempVector2, prepared.dimension);
        multiplyRowMajorVector(prepared.rMatrix.data, workspace.tempVector2, workspace.tempVector1, prepared.dimension);
        if (gradientAccumulator.length == 1) {
            double sum = 0.0;
            for (int i = 0; i < prepared.dimension; ++i) {
                sum += dLogL_df[i] - workspace.tempVector1[i];
            }
            gradientAccumulator[0] += sum;
            return;
        }
        if (gradientAccumulator.length != prepared.dimension) {
            throw new IllegalArgumentException(
                    "Stationary-mean gradient length must be 1 or "
                            + prepared.dimension + ", found " + gradientAccumulator.length);
        }
        for (int i = 0; i < prepared.dimension; ++i) {
            gradientAccumulator[i] += dLogL_df[i] - workspace.tempVector1[i];
        }
    }

    void accumulateMeanGradientPreparedDBasis(final OrthogonalBlockPreparedBranchBasis prepared,
                                             final double[] dLogL_df,
                                             final double[] dBasisGradientAccumulator,
                                             final OrthogonalBlockBranchGradientWorkspace workspace) {
        multiplyRowMajorVector(prepared.rtMatrix.data, dLogL_df, workspace.tempVector1, prepared.dimension);
        applyBlockDiagonalTransposeExp(
                prepared.expD.data, workspace.tempVector1, workspace.tempVector2, prepared.dimension);
        for (int i = 0; i < prepared.dimension; ++i) {
            dBasisGradientAccumulator[i] += workspace.tempVector1[i] - workspace.tempVector2[i];
        }
    }

    void finishMeanGradientFromDBasis(final double[] dBasisGradientAccumulator,
                                      final double[] gradientAccumulator,
                                      final OrthogonalBlockBranchGradientWorkspace workspace) {
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

    private void accumulateTransitionPullback(final DenseMatrix64F rMatrix,
                                              final DenseMatrix64F rtMatrix,
                                              final DenseMatrix64F expD,
                                              final double[] blockDParams,
                                              final double dt,
                                              final DenseMatrix64F upstreamF,
                                              final DenseMatrix64F upstreamFD,
                                              final BlockDiagonalFrechetHelper frechetHelper,
                                              final DenseMatrix64F temp1,
                                              final DenseMatrix64F temp2,
                                              final DenseMatrix64F temp3,
                                              final DenseMatrix64F gradD,
                                              final DenseMatrix64F gradR,
                                              final double[] scaledNegativeBlockDScratch,
                                              final double[] denseAdjointScratch,
                                              final double[] compressedDAccumulator,
                                              final double[] rotationAccumulator) {
        CommonOps.mult(rtMatrix, upstreamF, temp1);
        CommonOps.mult(temp1, rMatrix, upstreamFD);
        fillTransitionDGradient(
                blockDParams,
                dt,
                upstreamFD,
                frechetHelper,
                temp3,
                gradD,
                scaledNegativeBlockDScratch,
                denseAdjointScratch);
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        CommonOps.multTransB(rMatrix, expD, temp1);
        CommonOps.mult(upstreamF, temp1, gradR);
        CommonOps.multTransA(upstreamF, rMatrix, temp1);
        CommonOps.mult(temp1, expD, temp2);
        CommonOps.addEquals(gradR, temp2);
        addDenseMatrixToFlatArray(gradR, rotationAccumulator);
    }

    private void fillTransitionDGradient(final double[] blockDParams,
                                         final double dt,
                                         final DenseMatrix64F upstreamFD,
                                         final BlockDiagonalFrechetHelper frechetHelper,
                                         final DenseMatrix64F temp,
                                         final DenseMatrix64F gradD,
                                         final double[] scaledNegativeBlockDScratch,
                                         final double[] denseAdjointScratch) {
        if (denseFallbackPolicy.forceDenseAdjointExp()) {
            fillDenseAdjointExpGradient(blockDParams, dt, upstreamFD, scaledNegativeBlockDScratch,
                    denseAdjointScratch, gradD);
            return;
        }
        final DenseMatrix64F frechetInputTransition;
        if (denseFallbackPolicy.transposeNativeFrechetInput()) {
            CommonOps.transpose(upstreamFD, temp);
            frechetInputTransition = temp;
        } else {
            frechetInputTransition = upstreamFD;
        }
        frechetHelper.frechetAdjointExpInDBasis(blockDParams, frechetInputTransition, dt, gradD);
        transposeInPlace(gradD);
        if (!isFinite(gradD.data)) {
            fillDenseAdjointExpGradient(blockDParams, dt, upstreamFD, scaledNegativeBlockDScratch,
                    denseAdjointScratch, gradD);
        }
    }

    private static void fillDenseAdjointExpGradient(final double[] blockDParams,
                                                    final double dt,
                                                    final DenseMatrix64F upstreamFD,
                                                    final double[] scaledNegativeBlockDScratch,
                                                    final double[] denseAdjointScratch,
                                                    final DenseMatrix64F gradD) {
        fillScaledNegativeBlockDMatrix(blockDParams, dt, scaledNegativeBlockDScratch, upstreamFD.numRows);
        MatrixOps.toFlat(upstreamFD, denseAdjointScratch, upstreamFD.numRows);
        MatrixExponentialUtils.adjointExpFlat(
                scaledNegativeBlockDScratch,
                denseAdjointScratch,
                denseAdjointScratch,
                upstreamFD.numRows);
        fillDenseMatrix(denseAdjointScratch, gradD);
        transposeInPlace(gradD);
        CommonOps.scale(-dt, gradD);
    }

    private void fillTotalUpstreamOnTransition(final double[] stationaryMean,
                                               final double[] dLogL_dF,
                                               final double[] dLogL_df,
                                               final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] outData = out.data;
        for (int row = 0; row < dimension; ++row) {
            final int rowOff = row * dimension;
            for (int col = 0; col < dimension; ++col) {
                outData[rowOff + col] =
                        dLogL_dF[rowOff + col] - dLogL_df[row] * stationaryMean[col];
            }
        }
    }

    private void accumulateCompressedGradient(final DenseMatrix64F denseGradient,
                                              final double[] compressedAccumulator) {
        final double[] data = denseGradient.data;
        final int dimension = blockParameter.getRowDimension();
        final int upperBase = dimension;
        final int lowerBase = dimension + blockParameter.getNum2x2Blocks();
        int blockIndex = 0;

        for (int i = 0; i < dimension; ++i) {
            compressedAccumulator[i] += data[i * dimension + i];
        }
        for (int b = 0; b < blockParameter.getNumBlocks(); ++b) {
            if (blockParameter.getBlockSizes()[b] != 2) {
                continue;
            }
            final int start = blockParameter.getBlockStarts()[b];
            compressedAccumulator[upperBase + blockIndex] += data[start * dimension + (start + 1)];
            compressedAccumulator[lowerBase + blockIndex] += data[(start + 1) * dimension + start];
            blockIndex++;
        }
    }

    private void applyBlockDiagonalTransposeExp(final double[] expData,
                                                final double[] in,
                                                final double[] out,
                                                final int dimension) {
        for (int b = 0; b < blockParameter.getNumBlocks(); ++b) {
            final int start = blockParameter.getBlockStarts()[b];
            final int size = blockParameter.getBlockSizes()[b];
            if (size == 1) {
                out[start] = expData[start * dimension + start] * in[start];
            } else {
                final int row0 = start * dimension;
                final int row1 = (start + 1) * dimension;
                final double e00 = expData[row0 + start];
                final double e01 = expData[row0 + start + 1];
                final double e10 = expData[row1 + start];
                final double e11 = expData[row1 + start + 1];
                out[start] = e00 * in[start] + e10 * in[start + 1];
                out[start + 1] = e01 * in[start] + e11 * in[start + 1];
            }
        }
    }

    private static void addDenseMatrixToFlatArray(final DenseMatrix64F src, final double[] dest) {
        final int length = src.numRows * src.numCols;
        final double[] data = src.data;
        for (int i = 0; i < length; ++i) {
            dest[i] += data[i];
        }
    }

    private static void fillScaledNegativeBlockDMatrix(final double[] blockDParams,
                                                       final double dt,
                                                       final double[] out,
                                                       final int dimension) {
        java.util.Arrays.fill(out, 0, dimension * dimension, 0.0);
        final int upperOffset = dimension;
        final int lowerOffset = dimension + (dimension - 1);
        for (int i = 0; i < dimension; ++i) {
            out[i * dimension + i] = -dt * blockDParams[i];
            if (i < dimension - 1) {
                out[i * dimension + i + 1] = -dt * blockDParams[upperOffset + i];
                out[(i + 1) * dimension + i] = -dt * blockDParams[lowerOffset + i];
            }
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

    private static void fillDenseMatrix(final double[] source,
                                        final DenseMatrix64F out) {
        MatrixOps.fromFlat(source, out, out.numRows);
    }

    private static void transposeInPlace(final DenseMatrix64F matrix) {
        final int dimension = matrix.numRows;
        final double[] data = matrix.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final int ij = i * dimension + j;
                final int ji = j * dimension + i;
                final double tmp = data[ij];
                data[ij] = data[ji];
                data[ji] = tmp;
            }
        }
    }

    private static boolean isFinite(final double[] values) {
        for (double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }
}
