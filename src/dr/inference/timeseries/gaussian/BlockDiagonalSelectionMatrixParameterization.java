package dr.inference.timeseries.gaussian;

import dr.evomodel.continuous.ou.DenseSelectionMatrixParameterization;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalFrechetHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovAdjointHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovSolver;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Block-diagonal selection-matrix parametrization with a specialized forward matrix
 * exponential for {@code A = R D R^{-1}}.
 */
public final class BlockDiagonalSelectionMatrixParameterization extends DenseSelectionMatrixParameterization {

    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final BlockDiagonalExpSolver expSolver;
    private final BlockDiagonalFrechetHelper frechetHelper;
    private final BlockDiagonalLyapunovSolver lyapunovSolver;
    private final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper;
    private final double[] rData;
    private final double[] rinvData;
    private final double[] blockDParams;
    private final double[] expD;
    private final DenseMatrix64F rMatrix;
    private final DenseMatrix64F rinvMatrix;
    private final DenseMatrix64F transitionMatrix;
    private final DenseMatrix64F qMatrix;
    private final DenseMatrix64F qDBasis;
    private final DenseMatrix64F stationaryCovDBasis;
    private final DenseMatrix64F transitionCovDBasis;
    private final DenseMatrix64F upstreamF;
    private final DenseMatrix64F upstreamFD;
    private final DenseMatrix64F gradD;
    private final DenseMatrix64F gradR;
    private final DenseMatrix64F gV;
    private final DenseMatrix64F hDBasis;
    private final DenseMatrix64F gS;
    private final DenseMatrix64F yAdjoint;
    private final DenseMatrix64F gECov;
    private final DenseMatrix64F tmpMatrix1;
    private final DenseMatrix64F tmpMatrix2;
    private final DenseMatrix64F tmpMatrix3;
    private final double[][] workMatrix;

    public BlockDiagonalSelectionMatrixParameterization(
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter) {
        super(blockParameter);
        this.blockParameter = blockParameter;
        this.expSolver = new BlockDiagonalExpSolver(
                new BlockDiagonalExpSolver.BlockStructure(
                        blockParameter.getRowDimension(),
                        blockParameter.getBlockStarts(),
                        blockParameter.getBlockSizes()));
        this.frechetHelper = new BlockDiagonalFrechetHelper(expSolver.getStructure());
        this.lyapunovSolver = new BlockDiagonalLyapunovSolver(
                getDimension(),
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes());
        this.lyapunovAdjointHelper = new BlockDiagonalLyapunovAdjointHelper(getDimension(), lyapunovSolver);

        final int dimension = getDimension();
        this.rData = new double[dimension * dimension];
        this.rinvData = new double[dimension * dimension];
        this.blockDParams = new double[blockParameter.getTridiagonalDDimension()];
        this.expD = new double[blockParameter.getCompressedDDimension()];
        this.rMatrix = new DenseMatrix64F(dimension, dimension);
        this.rinvMatrix = new DenseMatrix64F(dimension, dimension);
        this.transitionMatrix = new DenseMatrix64F(dimension, dimension);
        this.qMatrix = new DenseMatrix64F(dimension, dimension);
        this.qDBasis = new DenseMatrix64F(dimension, dimension);
        this.stationaryCovDBasis = new DenseMatrix64F(dimension, dimension);
        this.transitionCovDBasis = new DenseMatrix64F(dimension, dimension);
        this.upstreamF = new DenseMatrix64F(dimension, dimension);
        this.upstreamFD = new DenseMatrix64F(dimension, dimension);
        this.gradD = new DenseMatrix64F(dimension, dimension);
        this.gradR = new DenseMatrix64F(dimension, dimension);
        this.gV = new DenseMatrix64F(dimension, dimension);
        this.hDBasis = new DenseMatrix64F(dimension, dimension);
        this.gS = new DenseMatrix64F(dimension, dimension);
        this.yAdjoint = new DenseMatrix64F(dimension, dimension);
        this.gECov = new DenseMatrix64F(dimension, dimension);
        this.tmpMatrix1 = new DenseMatrix64F(dimension, dimension);
        this.tmpMatrix2 = new DenseMatrix64F(dimension, dimension);
        this.tmpMatrix3 = new DenseMatrix64F(dimension, dimension);
        this.workMatrix = new double[dimension][dimension];
    }

    @Override
    public void fillTransitionMatrix(final double dt, final double[][] out) {
        refreshBasisCaches(dt);
        copyDenseMatrixToArray(transitionMatrix, out);
    }

    public void fillTransitionCovariance(final MatrixParameterInterface diffusionMatrix,
                                         final double dt,
                                         final double[][] out) {
        refreshBasisCaches(dt);
        fillDenseMatrix(diffusionMatrix, qMatrix);
        CommonOps.mult(rinvMatrix, qMatrix, tmpMatrix1);
        CommonOps.multTransB(tmpMatrix1, rinvMatrix, qDBasis);
        lyapunovSolver.solve(blockDParams, qDBasis, stationaryCovDBasis);
        multiplyBlockDiagonalLeft(expD, stationaryCovDBasis, tmpMatrix1, false);
        multiplyRightBlockDiagonal(tmpMatrix1, expD, tmpMatrix2, true);
        transitionCovDBasis.set(stationaryCovDBasis);
        CommonOps.subtractEquals(transitionCovDBasis, tmpMatrix2);
        CommonOps.mult(rMatrix, transitionCovDBasis, tmpMatrix1);
        CommonOps.multTransB(tmpMatrix1, rMatrix, tmpMatrix2);
        symmetrize(tmpMatrix2);
        copyDenseMatrixToArray(tmpMatrix2, out);
    }

    public void accumulateNativeGradientFromTransitionFlat(final double dt,
                                                           final double[] stationaryMean,
                                                           final double[] dLogL_dF,
                                                           final double[] dLogL_df,
                                                           final double[] compressedDAccumulator,
                                                           final double[] rotationAccumulator) {
        refreshBasisCaches(dt);
        fillTotalUpstreamOnTransition(stationaryMean, dLogL_dF, dLogL_df, upstreamF);

        CommonOps.multTransA(rMatrix, upstreamF, tmpMatrix1);
        CommonOps.multTransB(tmpMatrix1, rinvMatrix, upstreamFD);
        frechetHelper.frechetAdjointExpInDBasis(blockDParams, upstreamFD, dt, gradD);
        transposeInPlace(gradD);
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        CommonOps.multTransB(upstreamF, transitionMatrix, tmpMatrix1);
        CommonOps.multTransA(transitionMatrix, upstreamF, tmpMatrix2);
        CommonOps.subtract(tmpMatrix1, tmpMatrix2, tmpMatrix3);
        CommonOps.multTransB(tmpMatrix3, rinvMatrix, gradR);
        addDenseMatrixToFlatArray(gradR, rotationAccumulator);
    }

    public void accumulateNativeGradientFromCovarianceStationaryFlat(final MatrixParameterInterface diffusionMatrix,
                                                                     final double dt,
                                                                     final double[] dLogL_dV,
                                                                     final double[] compressedDAccumulator,
                                                                     final double[] rotationAccumulator) {
        refreshBasisCaches(dt);
        fillDenseMatrix(diffusionMatrix, qMatrix);
        CommonOps.mult(rinvMatrix, qMatrix, tmpMatrix1);
        CommonOps.multTransB(tmpMatrix1, rinvMatrix, qDBasis);
        lyapunovSolver.solve(blockDParams, qDBasis, stationaryCovDBasis);
        multiplyBlockDiagonalLeft(expD, stationaryCovDBasis, tmpMatrix1, false);
        multiplyRightBlockDiagonal(tmpMatrix1, expD, transitionCovDBasis, true);
        CommonOps.subtract(stationaryCovDBasis, transitionCovDBasis, transitionCovDBasis);
        symmetrize(transitionCovDBasis);

        fillSymmetricDenseMatrixFlat(dLogL_dV, gV);

        CommonOps.multTransA(rMatrix, gV, tmpMatrix1);
        CommonOps.mult(tmpMatrix1, rMatrix, hDBasis);
        symmetrize(hDBasis);

        multiplyBlockDiagonalLeft(expD, hDBasis, tmpMatrix1, true);
        multiplyRightBlockDiagonal(tmpMatrix1, expD, gS, false);
        CommonOps.changeSign(gS);
        CommonOps.addEquals(gS, hDBasis);
        lyapunovAdjointHelper.accumulateLyapunovContributionInDBasis(
                stationaryCovDBasis, gS, blockDParams, gradD);
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        lyapunovAdjointHelper.solveAdjointInDBasis(gS, blockDParams, yAdjoint);
        CommonOps.mult(yAdjoint, qDBasis, tmpMatrix1);
        CommonOps.multTransA(yAdjoint, qDBasis, tmpMatrix2);
        CommonOps.add(tmpMatrix1, tmpMatrix2, tmpMatrix3);
        CommonOps.multTransA(rinvMatrix, tmpMatrix3, gradR);
        addDenseMatrixToFlatArray(gradR, rotationAccumulator);

        multiplyRightBlockDiagonal(hDBasis, expD, tmpMatrix1, false);
        CommonOps.mult(tmpMatrix1, stationaryCovDBasis, gECov);
        CommonOps.scale(-2.0, gECov);
        frechetHelper.frechetAdjointExpInDBasis(blockDParams, gECov, dt, gradD);
        transposeInPlace(gradD);
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        CommonOps.mult(gV, rMatrix, tmpMatrix1);
        CommonOps.mult(tmpMatrix1, transitionCovDBasis, gradR);
        CommonOps.scale(2.0, gradR);
        addDenseMatrixToFlatArray(gradR, rotationAccumulator);
    }

    private static void multiplyRowMajor(final double[] leftRowMajor,
                                         final double[] rightRowMajor,
                                         final int dimension,
                                         final double[][] out) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += leftRowMajor[i * dimension + k] * rightRowMajor[k * dimension + j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void multiply(final double[][] left,
                                 final double[] rightRowMajor,
                                 final int dimension,
                                 final double[][] out) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += left[i][k] * rightRowMajor[k * dimension + j];
                }
                out[i][j] = sum;
            }
        }
    }

    public AbstractBlockDiagonalTwoByTwoMatrixParameter getBlockParameter() {
        return blockParameter;
    }

    private void refreshBasisCaches(final double dt) {
        final int dimension = getDimension();
        blockParameter.fillRAndRinv(rData, rinvData);
        System.arraycopy(rData, 0, rMatrix.data, 0, dimension * dimension);
        System.arraycopy(rinvData, 0, rinvMatrix.data, 0, dimension * dimension);
        blockParameter.fillBlockDiagonalElements(blockDParams);
        expSolver.computeCompressed(blockDParams, dt, expD);

        multiplyRightBlockDiagonal(rMatrix, expD, tmpMatrix1, false);
        CommonOps.mult(tmpMatrix1, rinvMatrix, transitionMatrix);
    }

    private void multiplyRightBlockDiagonal(final DenseMatrix64F matrix,
                                            final double[] compressedBlockData,
                                            final DenseMatrix64F out,
                                            final boolean transposeBlock) {
        final int dimension = matrix.numRows;
        final double[] matrixData = matrix.data;
        final double[] outData = out.data;
        final int upperBase = dimension;
        final int lowerBase = dimension + blockParameter.getNum2x2Blocks();
        int block2x2Index = 0;
        for (int row = 0; row < dimension; ++row) {
            final int rowOffset = row * dimension;
            for (int b = 0; b < blockParameter.getNumBlocks(); ++b) {
                final int start = blockParameter.getBlockStarts()[b];
                final int size = blockParameter.getBlockSizes()[b];
                if (size == 1) {
                    outData[rowOffset + start] =
                            matrixData[rowOffset + start] * compressedBlockData[start];
                } else {
                    final double e00 = compressedBlockData[start];
                    final double e11 = compressedBlockData[start + 1];
                    final double e01 = compressedBlockData[upperBase + block2x2Index];
                    final double e10 = compressedBlockData[lowerBase + block2x2Index];
                    final double x0 = matrixData[rowOffset + start];
                    final double x1 = matrixData[rowOffset + start + 1];
                    if (transposeBlock) {
                        outData[rowOffset + start] = x0 * e00 + x1 * e01;
                        outData[rowOffset + start + 1] = x0 * e10 + x1 * e11;
                    } else {
                        outData[rowOffset + start] = x0 * e00 + x1 * e10;
                        outData[rowOffset + start + 1] = x0 * e01 + x1 * e11;
                    }
                    block2x2Index++;
                }
            }
            block2x2Index = 0;
        }
    }

    private void multiplyBlockDiagonalLeft(final double[] compressedBlockData,
                                           final DenseMatrix64F matrix,
                                           final DenseMatrix64F out,
                                           final boolean transposeBlock) {
        final int dimension = matrix.numRows;
        final double[] matrixData = matrix.data;
        final double[] outData = out.data;
        final int upperBase = dimension;
        final int lowerBase = dimension + blockParameter.getNum2x2Blocks();
        int block2x2Index = 0;
        for (int b = 0; b < blockParameter.getNumBlocks(); ++b) {
            final int start = blockParameter.getBlockStarts()[b];
            final int size = blockParameter.getBlockSizes()[b];
            if (size == 1) {
                final int rowOffset = start * dimension;
                for (int col = 0; col < dimension; ++col) {
                    outData[rowOffset + col] =
                            compressedBlockData[start] * matrixData[rowOffset + col];
                }
            } else {
                final int row0 = start * dimension;
                final int row1 = (start + 1) * dimension;
                final double e00 = compressedBlockData[start];
                final double e11 = compressedBlockData[start + 1];
                final double e01 = compressedBlockData[upperBase + block2x2Index];
                final double e10 = compressedBlockData[lowerBase + block2x2Index];
                for (int col = 0; col < dimension; ++col) {
                    final double x0 = matrixData[row0 + col];
                    final double x1 = matrixData[row1 + col];
                    if (transposeBlock) {
                        outData[row0 + col] = e00 * x0 + e10 * x1;
                        outData[row1 + col] = e01 * x0 + e11 * x1;
                    } else {
                        outData[row0 + col] = e00 * x0 + e01 * x1;
                        outData[row1 + col] = e10 * x0 + e11 * x1;
                    }
                }
                block2x2Index++;
            }
        }
    }

    private void fillTotalUpstreamOnTransition(final double[] stationaryMean,
                                               final double[] dLogL_dF,
                                               final double[] dLogL_df,
                                               final DenseMatrix64F out) {
        final int dimension = getDimension();
        final double[] outData = out.data;
        for (int row = 0; row < dimension; ++row) {
            final int rowOffset = row * dimension;
            for (int col = 0; col < dimension; ++col) {
                outData[rowOffset + col] =
                        dLogL_dF[rowOffset + col] - dLogL_df[row] * stationaryMean[col];
            }
        }
    }

    private void accumulateCompressedGradient(final DenseMatrix64F denseGradient,
                                              final double[] compressedAccumulator) {
        CommonOps.fill(gradD, 0.0);
        gradD.set(denseGradient);
        final double[] data = gradD.data;
        final int dimension = getDimension();
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

    private static void addDenseMatrixToFlatArray(final DenseMatrix64F src, final double[] dest) {
        final int dimension = src.numRows;
        final double[] data = src.data;
        for (int i = 0; i < dimension * dimension; ++i) {
            dest[i] += data[i];
        }
    }

    private static void fillDenseMatrix(final MatrixParameterInterface parameter, final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                data[i * dimension + j] = parameter.getParameterValue(i, j);
            }
        }
    }

    private static void fillDenseMatrix(final double[][] source, final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                data[i * dimension + j] = source[i][j];
            }
        }
    }

    private static void fillSymmetricDenseMatrixFlat(final double[] source, final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                data[rowOffset + j] = 0.5 * (source[rowOffset + j] + source[j * dimension + i]);
            }
        }
    }

    private static void copyDenseMatrixToArray(final DenseMatrix64F source, final double[][] out) {
        final int dimension = source.numRows;
        final double[] data = source.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i][j] = data[i * dimension + j];
            }
        }
    }

    private static void symmetrize(final DenseMatrix64F matrix) {
        final int dimension = matrix.numRows;
        final double[] data = matrix.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final double value = 0.5 * (data[i * dimension + j] + data[j * dimension + i]);
                data[i * dimension + j] = value;
                data[j * dimension + i] = value;
            }
        }
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
}
