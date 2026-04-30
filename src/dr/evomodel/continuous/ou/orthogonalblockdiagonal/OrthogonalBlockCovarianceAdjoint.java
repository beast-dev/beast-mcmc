package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.continuous.ou.MatrixExponentialUtils;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalFrechetHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovAdjointHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovSolver;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Pullback for orthogonal block OU transition covariance and diffusion adjoints.
 */
final class OrthogonalBlockCovarianceAdjoint {

    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final BlockDiagonalLyapunovSolver lyapunovSolver;
    private final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper;
    private final BlockDiagonalFrechetHelper frechetHelper;
    private final OrthogonalBlockDenseFallbackPolicy denseFallbackPolicy;
    private final DenseMatrix64F qMatrix;
    private final DenseMatrix64F qDBasis;
    private final DenseMatrix64F stationaryCovDBasis;
    private final DenseMatrix64F transitionCovDBasis;
    private final DenseMatrix64F transitionCovariance;
    private final DenseMatrix64F gV;
    private final DenseMatrix64F hDBasis;
    private final DenseMatrix64F gS;
    private final DenseMatrix64F yAdjoint;
    private final DenseMatrix64F gECov;
    private final DenseMatrix64F gradD;
    private final DenseMatrix64F gradR;
    private final DenseMatrix64F temp1;
    private final DenseMatrix64F temp2;
    private final DenseMatrix64F temp3;
    private final double[][] scaledNegativeBlockDScratch;
    private final double[][] denseAdjointScratch;

    OrthogonalBlockCovarianceAdjoint(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                     final BlockDiagonalLyapunovSolver lyapunovSolver,
                                     final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper,
                                     final BlockDiagonalFrechetHelper frechetHelper,
                                     final OrthogonalBlockDenseFallbackPolicy denseFallbackPolicy) {
        this.blockParameter = blockParameter;
        this.lyapunovSolver = lyapunovSolver;
        this.lyapunovAdjointHelper = lyapunovAdjointHelper;
        this.frechetHelper = frechetHelper;
        this.denseFallbackPolicy = denseFallbackPolicy;
        final int d = blockParameter.getRowDimension();
        this.qMatrix = new DenseMatrix64F(d, d);
        this.qDBasis = new DenseMatrix64F(d, d);
        this.stationaryCovDBasis = new DenseMatrix64F(d, d);
        this.transitionCovDBasis = new DenseMatrix64F(d, d);
        this.transitionCovariance = new DenseMatrix64F(d, d);
        this.gV = new DenseMatrix64F(d, d);
        this.hDBasis = new DenseMatrix64F(d, d);
        this.gS = new DenseMatrix64F(d, d);
        this.yAdjoint = new DenseMatrix64F(d, d);
        this.gECov = new DenseMatrix64F(d, d);
        this.gradD = new DenseMatrix64F(d, d);
        this.gradR = new DenseMatrix64F(d, d);
        this.temp1 = new DenseMatrix64F(d, d);
        this.temp2 = new DenseMatrix64F(d, d);
        this.temp3 = new DenseMatrix64F(d, d);
        this.scaledNegativeBlockDScratch = new double[d][d];
        this.denseAdjointScratch = new double[d][d];
    }

    void fillTransitionCovarianceMatrix(final MatrixParameterInterface diffusionMatrix,
                                        final OrthogonalBlockBasisCache basis,
                                        final DenseMatrix64F out) {
        OrthogonalBlockTransitionCovarianceSolver.fillTransitionCovariance(
                diffusionMatrix,
                basis.rMatrix,
                basis.rtMatrix,
                basis.expD,
                basis.blockDParams,
                lyapunovSolver,
                qMatrix,
                qDBasis,
                stationaryCovDBasis,
                transitionCovDBasis,
                temp1,
                out,
                false);
    }

    void fillTransitionCovarianceMatrixPrepared(final OrthogonalBlockPreparedBranchBasis prepared,
                                                final MatrixParameterInterface diffusionMatrix,
                                                final OrthogonalBlockBranchGradientWorkspace workspace,
                                                final DenseMatrix64F out) {
        OrthogonalBlockTransitionCovarianceSolver.fillTransitionCovariance(
                diffusionMatrix,
                prepared.rMatrix,
                prepared.rtMatrix,
                prepared.expD,
                prepared.blockDParams,
                workspace.lyapunovSolver,
                workspace.qMatrix,
                workspace.qDBasis,
                workspace.stationaryCovDBasis,
                workspace.transitionCovDBasis,
                workspace.temp1,
                out,
                true);
    }

    void prepareAndAccumulateCurrent(final MatrixParameterInterface diffusionMatrix,
                                     final OrthogonalBlockBasisCache basis,
                                     final double dt,
                                     final double[] dLogL_dV,
                                     final double[] compressedDAccumulator,
                                     final double[][] rotationAccumulator) {
        fillTransitionCovarianceMatrix(diffusionMatrix, basis, transitionCovariance);
        accumulateCurrentCached(basis, dt, dLogL_dV, compressedDAccumulator, rotationAccumulator);
    }

    void accumulateCurrentCached(final OrthogonalBlockBasisCache basis,
                                 final double dt,
                                 final double[] dLogL_dV,
                                 final double[] compressedDAccumulator,
                                 final double[][] rotationAccumulator) {
        accumulateCovariancePullback(
                basis.rMatrix,
                basis.rtMatrix,
                basis.expD,
                basis.blockDParams,
                dt,
                qMatrix,
                stationaryCovDBasis,
                transitionCovDBasis,
                gV,
                hDBasis,
                gS,
                yAdjoint,
                gECov,
                gradD,
                gradR,
                temp1,
                temp2,
                temp3,
                lyapunovAdjointHelper,
                frechetHelper,
                scaledNegativeBlockDScratch,
                denseAdjointScratch,
                dLogL_dV,
                compressedDAccumulator,
                rotationAccumulator);
    }

    void accumulatePrepared(final OrthogonalBlockPreparedBranchBasis prepared,
                            final double[] dLogL_dV,
                            final OrthogonalBlockBranchGradientWorkspace workspace,
                            final double[] compressedDAccumulator,
                            final double[][] rotationAccumulator) {
        accumulateCovariancePullback(
                prepared.rMatrix,
                prepared.rtMatrix,
                prepared.expD,
                prepared.blockDParams,
                prepared.dt,
                workspace.qMatrix,
                workspace.stationaryCovDBasis,
                workspace.transitionCovDBasis,
                workspace.gV,
                workspace.hDBasis,
                workspace.gS,
                workspace.yAdjoint,
                workspace.gECov,
                workspace.gradD,
                workspace.gradR,
                workspace.temp1,
                workspace.temp2,
                workspace.temp3,
                workspace.lyapunovAdjointHelper,
                workspace.frechetHelper,
                workspace.scaledNegativeBlockDScratch,
                workspace.denseAdjointScratch,
                dLogL_dV,
                compressedDAccumulator,
                rotationAccumulator);
    }

    void accumulatePreparedFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                final double[] dLogL_dV,
                                final OrthogonalBlockBranchGradientWorkspace workspace,
                                final double[] compressedDAccumulator,
                                final double[] rotationAccumulator) {
        accumulateCovariancePullback(
                prepared.rMatrix,
                prepared.rtMatrix,
                prepared.expD,
                prepared.blockDParams,
                prepared.dt,
                workspace.qMatrix,
                workspace.stationaryCovDBasis,
                workspace.transitionCovDBasis,
                workspace.gV,
                workspace.hDBasis,
                workspace.gS,
                workspace.yAdjoint,
                workspace.gECov,
                workspace.gradD,
                workspace.gradR,
                workspace.temp1,
                workspace.temp2,
                workspace.temp3,
                workspace.lyapunovAdjointHelper,
                workspace.frechetHelper,
                workspace.scaledNegativeBlockDScratch,
                workspace.denseAdjointScratch,
                dLogL_dV,
                compressedDAccumulator,
                rotationAccumulator);
    }

    void accumulateDiffusionGradientCurrent(final OrthogonalBlockBasisCache basis,
                                            final double[][] dLogL_dV,
                                            final double[] gradientAccumulator) {
        fillSymmetricDenseMatrix(dLogL_dV, gV);
        accumulateDiffusionGradientSymmetric(
                basis.rMatrix,
                basis.rtMatrix,
                basis.expD,
                basis.blockDParams,
                gV,
                hDBasis,
                gS,
                yAdjoint,
                temp1,
                temp2,
                lyapunovAdjointHelper,
                gradientAccumulator);
    }

    void accumulateDiffusionGradientPrepared(final OrthogonalBlockPreparedBranchBasis prepared,
                                             final double[][] dLogL_dV,
                                             final double[] gradientAccumulator,
                                             final OrthogonalBlockBranchGradientWorkspace workspace) {
        fillSymmetricDenseMatrix(dLogL_dV, workspace.gV);
        accumulateDiffusionGradientPreparedSymmetric(prepared, gradientAccumulator, workspace);
    }

    void accumulateDiffusionGradientPreparedFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                                 final double[] dLogL_dV,
                                                 final double[] gradientAccumulator,
                                                 final OrthogonalBlockBranchGradientWorkspace workspace) {
        fillSymmetricDenseMatrixFlat(dLogL_dV, workspace.gV);
        accumulateDiffusionGradientPreparedSymmetric(prepared, gradientAccumulator, workspace);
    }

    private void accumulateDiffusionGradientPreparedSymmetric(final OrthogonalBlockPreparedBranchBasis prepared,
                                                              final double[] gradientAccumulator,
                                                              final OrthogonalBlockBranchGradientWorkspace workspace) {
        accumulateDiffusionGradientSymmetric(
                prepared.rMatrix,
                prepared.rtMatrix,
                prepared.expD,
                prepared.blockDParams,
                workspace.gV,
                workspace.hDBasis,
                workspace.gS,
                workspace.yAdjoint,
                workspace.temp1,
                workspace.temp2,
                workspace.lyapunovAdjointHelper,
                gradientAccumulator);
    }

    private void accumulateCovariancePullback(final DenseMatrix64F rMatrix,
                                              final DenseMatrix64F rtMatrix,
                                              final DenseMatrix64F expD,
                                              final double[] blockDParams,
                                              final double dt,
                                              final DenseMatrix64F qMatrix,
                                              final DenseMatrix64F stationaryCovDBasis,
                                              final DenseMatrix64F transitionCovDBasis,
                                              final DenseMatrix64F gV,
                                              final DenseMatrix64F hDBasis,
                                              final DenseMatrix64F gS,
                                              final DenseMatrix64F yAdjoint,
                                              final DenseMatrix64F gECov,
                                              final DenseMatrix64F gradD,
                                              final DenseMatrix64F gradR,
                                              final DenseMatrix64F temp1,
                                              final DenseMatrix64F temp2,
                                              final DenseMatrix64F temp3,
                                              final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper,
                                              final BlockDiagonalFrechetHelper frechetHelper,
                                              final double[][] scaledNegativeBlockDScratch,
                                              final double[][] denseAdjointScratch,
                                              final double[] dLogL_dV,
                                              final double[] compressedDAccumulator,
                                              final double[][] rotationAccumulator) {
        fillSymmetricDenseMatrixFlat(dLogL_dV, gV);
        accumulateCovariancePullbackSymmetric(rMatrix, rtMatrix, expD, blockDParams, dt,
                qMatrix, stationaryCovDBasis, transitionCovDBasis, gV, hDBasis, gS,
                yAdjoint, gECov, gradD, gradR, temp1, temp2, temp3,
                lyapunovAdjointHelper, frechetHelper,
                scaledNegativeBlockDScratch, denseAdjointScratch, compressedDAccumulator);
        addDenseMatrixToArray(gradR, rotationAccumulator);
    }

    private void accumulateCovariancePullback(final DenseMatrix64F rMatrix,
                                              final DenseMatrix64F rtMatrix,
                                              final DenseMatrix64F expD,
                                              final double[] blockDParams,
                                              final double dt,
                                              final DenseMatrix64F qMatrix,
                                              final DenseMatrix64F stationaryCovDBasis,
                                              final DenseMatrix64F transitionCovDBasis,
                                              final DenseMatrix64F gV,
                                              final DenseMatrix64F hDBasis,
                                              final DenseMatrix64F gS,
                                              final DenseMatrix64F yAdjoint,
                                              final DenseMatrix64F gECov,
                                              final DenseMatrix64F gradD,
                                              final DenseMatrix64F gradR,
                                              final DenseMatrix64F temp1,
                                              final DenseMatrix64F temp2,
                                              final DenseMatrix64F temp3,
                                              final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper,
                                              final BlockDiagonalFrechetHelper frechetHelper,
                                              final double[][] scaledNegativeBlockDScratch,
                                              final double[][] denseAdjointScratch,
                                              final double[] dLogL_dV,
                                              final double[] compressedDAccumulator,
                                              final double[] rotationAccumulator) {
        fillSymmetricDenseMatrixFlat(dLogL_dV, gV);
        accumulateCovariancePullbackSymmetric(rMatrix, rtMatrix, expD, blockDParams, dt,
                qMatrix, stationaryCovDBasis, transitionCovDBasis, gV, hDBasis, gS,
                yAdjoint, gECov, gradD, gradR, temp1, temp2, temp3,
                lyapunovAdjointHelper, frechetHelper,
                scaledNegativeBlockDScratch, denseAdjointScratch, compressedDAccumulator);
        addDenseMatrixToFlatArray(gradR, rotationAccumulator);
    }

    private void accumulateCovariancePullbackSymmetric(final DenseMatrix64F rMatrix,
                                                       final DenseMatrix64F rtMatrix,
                                                       final DenseMatrix64F expD,
                                                       final double[] blockDParams,
                                                       final double dt,
                                                       final DenseMatrix64F qMatrix,
                                                       final DenseMatrix64F stationaryCovDBasis,
                                                       final DenseMatrix64F transitionCovDBasis,
                                                       final DenseMatrix64F gV,
                                                       final DenseMatrix64F hDBasis,
                                                       final DenseMatrix64F gS,
                                                       final DenseMatrix64F yAdjoint,
                                                       final DenseMatrix64F gECov,
                                                       final DenseMatrix64F gradD,
                                                       final DenseMatrix64F gradR,
                                                       final DenseMatrix64F temp1,
                                                       final DenseMatrix64F temp2,
                                                       final DenseMatrix64F temp3,
                                                       final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper,
                                                       final BlockDiagonalFrechetHelper frechetHelper,
                                                       final double[][] scaledNegativeBlockDScratch,
                                                       final double[][] denseAdjointScratch,
                                                       final double[] compressedDAccumulator) {
        CommonOps.mult(rtMatrix, gV, temp1);
        CommonOps.mult(temp1, rMatrix, hDBasis);
        symmetrize(hDBasis);

        CommonOps.multTransA(expD, hDBasis, temp1);
        CommonOps.mult(temp1, expD, gS);
        CommonOps.changeSign(gS);
        CommonOps.addEquals(gS, hDBasis);

        CommonOps.fill(gradD, 0.0);
        lyapunovAdjointHelper.accumulateLyapunovContributionInDBasis(
                stationaryCovDBasis, gS, blockDParams, gradD);
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        lyapunovAdjointHelper.solveAdjointInDBasis(gS, blockDParams, yAdjoint);

        CommonOps.mult(qMatrix, rMatrix, temp1);
        CommonOps.multTransB(temp1, yAdjoint, gradR);
        CommonOps.multTransA(qMatrix, rMatrix, temp1);
        CommonOps.mult(temp1, yAdjoint, temp2);
        CommonOps.addEquals(gradR, temp2);
        CommonOps.scale(-1.0, gradR);

        CommonOps.mult(hDBasis, expD, temp1);
        CommonOps.mult(temp1, stationaryCovDBasis, gECov);
        CommonOps.scale(-2.0, gECov);
        fillCovarianceExpDGradient(
                blockDParams,
                dt,
                gECov,
                frechetHelper,
                temp3,
                gradD,
                scaledNegativeBlockDScratch,
                denseAdjointScratch);
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        CommonOps.mult(gV, rMatrix, temp1);
        CommonOps.mult(temp1, transitionCovDBasis, temp2);
        CommonOps.multTransA(gV, rMatrix, temp1);
        CommonOps.mult(temp1, transitionCovDBasis, temp3);
        CommonOps.addEquals(temp2, temp3);
        CommonOps.addEquals(gradR, temp2);
    }

    private void fillCovarianceExpDGradient(final double[] blockDParams,
                                            final double dt,
                                            final DenseMatrix64F gECov,
                                            final BlockDiagonalFrechetHelper frechetHelper,
                                            final DenseMatrix64F temp,
                                            final DenseMatrix64F gradD,
                                            final double[][] scaledNegativeBlockDScratch,
                                            final double[][] denseAdjointScratch) {
        if (denseFallbackPolicy.forceDenseAdjointExp()) {
            fillDenseAdjointExpGradient(blockDParams, dt, gECov, scaledNegativeBlockDScratch,
                    denseAdjointScratch, gradD);
            return;
        }
        final DenseMatrix64F frechetInputCov;
        if (denseFallbackPolicy.transposeNativeFrechetInput()) {
            CommonOps.transpose(gECov, temp);
            frechetInputCov = temp;
        } else {
            frechetInputCov = gECov;
        }
        frechetHelper.frechetAdjointExpInDBasis(blockDParams, frechetInputCov, dt, gradD);
        transposeInPlace(gradD);
    }

    private static void accumulateDiffusionGradientSymmetric(final DenseMatrix64F rMatrix,
                                                             final DenseMatrix64F rtMatrix,
                                                             final DenseMatrix64F expD,
                                                             final double[] blockDParams,
                                                             final DenseMatrix64F gV,
                                                             final DenseMatrix64F hDBasis,
                                                             final DenseMatrix64F gS,
                                                             final DenseMatrix64F yAdjoint,
                                                             final DenseMatrix64F temp1,
                                                             final DenseMatrix64F temp2,
                                                             final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper,
                                                             final double[] gradientAccumulator) {
        CommonOps.mult(rtMatrix, gV, temp1);
        CommonOps.mult(temp1, rMatrix, hDBasis);
        symmetrize(hDBasis);

        CommonOps.multTransA(expD, hDBasis, temp1);
        CommonOps.mult(temp1, expD, gS);
        CommonOps.subtract(hDBasis, gS, gS);
        symmetrize(gS);

        yAdjoint.set(gS);
        CommonOps.scale(-1.0, yAdjoint);
        lyapunovAdjointHelper.solveAdjointInDBasis(yAdjoint, blockDParams, yAdjoint);

        CommonOps.mult(rMatrix, yAdjoint, temp1);
        CommonOps.mult(temp1, rtMatrix, temp2);
        symmetrize(temp2);

        final int dimension = temp2.numRows;
        final double[] data = temp2.data;
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                gradientAccumulator[rowOffset + j] += data[rowOffset + j];
            }
        }
    }

    private static void fillDenseAdjointExpGradient(final double[] blockDParams,
                                                    final double dt,
                                                    final DenseMatrix64F upstream,
                                                    final double[][] scaledNegativeBlockDScratch,
                                                    final double[][] denseAdjointScratch,
                                                    final DenseMatrix64F gradD) {
        fillScaledNegativeBlockDMatrix(blockDParams, dt, scaledNegativeBlockDScratch, upstream.numRows);
        copyDenseMatrixToArray(upstream, denseAdjointScratch);
        MatrixExponentialUtils.adjointExp(
                scaledNegativeBlockDScratch,
                denseAdjointScratch,
                denseAdjointScratch);
        fillDenseMatrix(denseAdjointScratch, gradD);
        transposeInPlace(gradD);
        CommonOps.scale(-dt, gradD);
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

    private static void fillDenseMatrix(final MatrixParameterInterface parameter, final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                data[i * dimension + j] = parameter.getParameterValue(i, j);
            }
        }
    }

    private static void fillSymmetricDenseMatrix(final double[][] source, final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                data[i * dimension + j] = 0.5 * (source[i][j] + source[j][i]);
            }
        }
    }

    private static void fillSymmetricDenseMatrixFlat(final double[] source, final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                data[i * dimension + j] = 0.5 * (source[i * dimension + j] + source[j * dimension + i]);
            }
        }
    }

    private static void addDenseMatrixToArray(final DenseMatrix64F src, final double[][] dest) {
        final int dimension = dest.length;
        final double[] data = src.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                dest[i][j] += data[i * dimension + j];
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

    private static void copyDenseMatrixToArray(final DenseMatrix64F source, final double[][] out) {
        final int dimension = source.numRows;
        final double[] data = source.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i][j] = data[i * dimension + j];
            }
        }
    }

    private static void fillScaledNegativeBlockDMatrix(final double[] blockDParams,
                                                       final double dt,
                                                       final double[][] out,
                                                       final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i][j] = 0.0;
            }
        }
        final int upperOffset = dimension;
        final int lowerOffset = dimension + (dimension - 1);
        for (int i = 0; i < dimension; ++i) {
            out[i][i] = -dt * blockDParams[i];
            if (i < dimension - 1) {
                out[i][i + 1] = -dt * blockDParams[upperOffset + i];
                out[i + 1][i] = -dt * blockDParams[lowerOffset + i];
            }
        }
    }

    private static void fillDenseMatrix(final double[][] source,
                                        final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                data[i * dimension + j] = source[i][j];
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
