package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.continuous.ou.MatrixExponentialUtils;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalFrechetHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovAdjointHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovSolver;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Pullback for orthogonal block OU transition covariance and diffusion adjoints.
 */
final class OrthogonalBlockCovarianceAdjoint {

    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final int[] blockStarts;
    private final int[] blockSizes;
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
    private final double[] scaledNegativeBlockDScratch;
    private final double[] denseAdjointScratch;

    OrthogonalBlockCovarianceAdjoint(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                     final BlockDiagonalLyapunovSolver lyapunovSolver,
                                     final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper,
                                     final BlockDiagonalFrechetHelper frechetHelper,
                                     final OrthogonalBlockDenseFallbackPolicy denseFallbackPolicy) {
        this.blockParameter = blockParameter;
        this.blockStarts = blockParameter.getBlockStarts();
        this.blockSizes = blockParameter.getBlockSizes();
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
        this.scaledNegativeBlockDScratch = new double[d * d];
        this.denseAdjointScratch = new double[d * d];
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
                false,
                blockStarts,
                blockSizes);
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
                true,
                blockStarts,
                blockSizes);
    }

    void prepareAndAccumulateCurrentFlat(final MatrixParameterInterface diffusionMatrix,
                                         final OrthogonalBlockBasisCache basis,
                                         final double dt,
                                         final double[] dLogL_dV,
                                         final double[] compressedDAccumulator,
                                         final double[] rotationAccumulator) {
        fillTransitionCovarianceMatrix(diffusionMatrix, basis, transitionCovariance);
        accumulateCurrentCachedFlat(basis, dt, dLogL_dV, compressedDAccumulator, rotationAccumulator);
    }

    void accumulateCurrentCachedFlat(final OrthogonalBlockBasisCache basis,
                                     final double dt,
                                     final double[] dLogL_dV,
                                     final double[] compressedDAccumulator,
                                     final double[] rotationAccumulator) {
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

    void accumulateDiffusionGradientCurrentFlat(final OrthogonalBlockBasisCache basis,
                                                final double[] dLogL_dV,
                                                final boolean transposeAdjoint,
                                                final double[] gradientAccumulator) {
        fillSymmetricDenseMatrixFlat(dLogL_dV, transposeAdjoint, gV);
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
                blockStarts,
                blockSizes,
                gradientAccumulator);
    }

    void accumulateDiffusionGradientPreparedFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                                 final double[] dLogL_dV,
                                                 final double[] gradientAccumulator,
                                                 final OrthogonalBlockBranchGradientWorkspace workspace) {
        fillSymmetricDenseMatrixFlat(dLogL_dV, workspace.gV);
        accumulateDiffusionGradientPreparedSymmetric(prepared, gradientAccumulator, workspace);
    }

    void accumulateDiffusionGradientPreparedDBasisFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                                       final double[] dLogL_dV,
                                                       final double[] dBasisGradientAccumulator,
                                                       final OrthogonalBlockBranchGradientWorkspace workspace) {
        fillSymmetricDenseMatrixFlat(dLogL_dV, workspace.gV);
        fillDiffusionGradientDBasisSymmetric(
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
                blockStarts,
                blockSizes);
        addDenseMatrixToFlatArray(workspace.yAdjoint, dBasisGradientAccumulator);
    }

    void accumulatePreparedCovarianceAndDiffusionGradientFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                                              final double[] dLogL_dV,
                                                              final OrthogonalBlockBranchGradientWorkspace workspace,
                                                              final double[] compressedDAccumulator,
                                                              final double[] rotationAccumulator,
                                                              final boolean delayDiffusionGradientRotation,
                                                              final double[] diffusionGradientAccumulator) {
        fillSymmetricDenseMatrixFlat(dLogL_dV, workspace.gV);
        fillSharedCovarianceAdjointsInDBasis(
                prepared.rMatrix,
                prepared.rtMatrix,
                prepared.expD,
                workspace.gV,
                workspace.hDBasis,
                workspace.gS,
                workspace.temp1,
                blockStarts,
                blockSizes);
        symmetrize(workspace.gS);

        accumulateCovariancePullbackFromSharedSymmetric(
                prepared.rMatrix,
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
                compressedDAccumulator);
        addDenseMatrixToFlatArray(workspace.gradR, rotationAccumulator);

        CommonOps.scale(-1.0, workspace.yAdjoint);
        if (delayDiffusionGradientRotation) {
            addDenseMatrixToFlatArray(workspace.yAdjoint, diffusionGradientAccumulator);
        } else {
            rotateDBasisDiffusionGradientToOriginalBasis(
                    prepared.rMatrix,
                    workspace.yAdjoint,
                    workspace.temp1,
                    workspace.temp2);
            addDenseMatrixToFlatArray(workspace.temp2, diffusionGradientAccumulator);
        }
    }

    void finishDiffusionGradientFromDBasisFlat(final DenseMatrix64F rMatrix,
                                               final DenseMatrix64F rtMatrix,
                                               final double[] dBasisGradientAccumulator,
                                               final double[] gradientAccumulator,
                                               final OrthogonalBlockBranchGradientWorkspace workspace) {
        fillDenseMatrix(dBasisGradientAccumulator, workspace.yAdjoint);
        rotateDBasisDiffusionGradientToOriginalBasis(
                rMatrix,
                workspace.yAdjoint,
                workspace.temp3,
                workspace.gS);
        addDenseMatrixToFlatArray(workspace.gS, gradientAccumulator);
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
                blockStarts,
                blockSizes,
                gradientAccumulator);
    }

    private void accumulateCovariancePullback(final DenseMatrix64F rMatrix,
                                              final DenseMatrix64F rtMatrix,
                                              final double[] expD,
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
                                              final double[] scaledNegativeBlockDScratch,
                                              final double[] denseAdjointScratch,
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
                                                       final double[] expD,
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
                                                       final double[] scaledNegativeBlockDScratch,
                                                       final double[] denseAdjointScratch,
                                                       final double[] compressedDAccumulator) {
        fillSharedCovarianceAdjointsInDBasis(
                rMatrix, rtMatrix, expD, gV, hDBasis, gS, temp1, blockStarts, blockSizes);
        accumulateCovariancePullbackFromSharedSymmetric(rMatrix, expD, blockDParams, dt,
                qMatrix, stationaryCovDBasis, transitionCovDBasis, gV, hDBasis, gS,
                yAdjoint, gECov, gradD, gradR, temp1, temp2, temp3,
                lyapunovAdjointHelper, frechetHelper,
                scaledNegativeBlockDScratch, denseAdjointScratch, compressedDAccumulator);
    }

    private void accumulateCovariancePullbackFromSharedSymmetric(final DenseMatrix64F rMatrix,
                                                                 final double[] expD,
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
                                                                 final double[] scaledNegativeBlockDScratch,
                                                                 final double[] denseAdjointScratch,
                                                                 final double[] compressedDAccumulator) {
        CommonOps.fill(gradD, 0.0);
        lyapunovAdjointHelper.accumulateLyapunovContributionInDBasis(
                stationaryCovDBasis, gS, blockDParams, gradD);
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        lyapunovAdjointHelper.solveAdjointInDBasis(gS, blockDParams, yAdjoint);

        OrthogonalBlockDenseMatrixOps.multSymmetricLeft(qMatrix, rMatrix, temp1);
        OrthogonalBlockDenseMatrixOps.mult(temp1, yAdjoint, gradR);
        CommonOps.scale(-2.0, gradR);

        multiplyBlockDiagonalRight(hDBasis, expD, temp1, blockStarts, blockSizes);
        OrthogonalBlockDenseMatrixOps.multSymmetricRight(temp1, stationaryCovDBasis, gECov);
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

        OrthogonalBlockDenseMatrixOps.multSymmetricLeft(gV, rMatrix, temp1);
        OrthogonalBlockDenseMatrixOps.multSymmetricRight(temp1, transitionCovDBasis, temp2);
        OrthogonalBlockDenseMatrixOps.addEquals(gradR, 2.0, temp2);
    }

    private void fillCovarianceExpDGradient(final double[] blockDParams,
                                            final double dt,
                                            final DenseMatrix64F gECov,
                                            final BlockDiagonalFrechetHelper frechetHelper,
                                            final DenseMatrix64F temp,
                                            final DenseMatrix64F gradD,
                                            final double[] scaledNegativeBlockDScratch,
                                            final double[] denseAdjointScratch) {
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
                                                             final double[] expD,
                                                             final double[] blockDParams,
                                                             final DenseMatrix64F gV,
                                                             final DenseMatrix64F hDBasis,
                                                             final DenseMatrix64F gS,
                                                             final DenseMatrix64F yAdjoint,
                                                             final DenseMatrix64F temp1,
                                                             final DenseMatrix64F temp2,
                                                             final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper,
                                                             final int[] blockStarts,
                                                             final int[] blockSizes,
                                                             final double[] gradientAccumulator) {
        fillDiffusionGradientDBasisSymmetric(
                rMatrix,
                rtMatrix,
                expD,
                blockDParams,
                gV,
                hDBasis,
                gS,
                yAdjoint,
                temp1,
                temp2,
                lyapunovAdjointHelper,
                blockStarts,
                blockSizes);

        rotateDBasisDiffusionGradientToOriginalBasis(rMatrix, yAdjoint, temp1, temp2);

        final int dimension = temp2.numRows;
        final double[] data = temp2.data;
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                gradientAccumulator[rowOffset + j] += data[rowOffset + j];
            }
        }
    }

    private static void fillDiffusionGradientDBasisSymmetric(final DenseMatrix64F rMatrix,
                                                             final DenseMatrix64F rtMatrix,
                                                             final double[] expD,
                                                             final double[] blockDParams,
                                                             final DenseMatrix64F gV,
                                                             final DenseMatrix64F hDBasis,
                                                             final DenseMatrix64F gS,
                                                             final DenseMatrix64F yAdjoint,
                                                             final DenseMatrix64F temp1,
                                                             final DenseMatrix64F temp2,
                                                             final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper,
                                                             final int[] blockStarts,
                                                             final int[] blockSizes) {
        fillSharedCovarianceAdjointsInDBasis(
                rMatrix, rtMatrix, expD, gV, hDBasis, gS, temp1, blockStarts, blockSizes);
        fillDiffusionGradientDBasisFromSharedSymmetric(
                blockDParams,
                gS,
                yAdjoint,
                lyapunovAdjointHelper);
    }

    private static void fillSharedCovarianceAdjointsInDBasis(final DenseMatrix64F rMatrix,
                                                             final DenseMatrix64F rtMatrix,
                                                             final double[] expD,
                                                             final DenseMatrix64F gV,
                                                             final DenseMatrix64F hDBasis,
                                                             final DenseMatrix64F gS,
                                                             final DenseMatrix64F temp,
                                                             final int[] blockStarts,
                                                             final int[] blockSizes) {
        OrthogonalBlockDenseMatrixOps.multSymmetricRight(rtMatrix, gV, temp);
        OrthogonalBlockDenseMatrixOps.mult(temp, rMatrix, hDBasis);
        symmetrize(hDBasis);

        multiplyBlockDiagonalLeftTranspose(expD, hDBasis, temp, blockStarts, blockSizes);
        multiplyBlockDiagonalRight(temp, expD, gS, blockStarts, blockSizes);
        OrthogonalBlockDenseMatrixOps.subtract(hDBasis, gS, gS);
    }

    private static void fillDiffusionGradientDBasisFromSharedSymmetric(final double[] blockDParams,
                                                                       final DenseMatrix64F gS,
                                                                       final DenseMatrix64F yAdjoint,
                                                                       final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper) {
        symmetrize(gS);
        yAdjoint.set(gS);
        CommonOps.scale(-1.0, yAdjoint);
        lyapunovAdjointHelper.solveAdjointInDBasis(yAdjoint, blockDParams, yAdjoint);
    }

    private static void multiplyBlockDiagonalLeftTranspose(final double[] blockDiagonal,
                                                           final DenseMatrix64F matrix,
                                                           final DenseMatrix64F out,
                                                           final int[] blockStarts,
                                                           final int[] blockSizes) {
        OrthogonalBlockMatrixOps.multiplyBlockDiagonalLeftTranspose(
                blockDiagonal, matrix.data, out.data, matrix.numRows, blockStarts, blockSizes);
    }

    private static void multiplyBlockDiagonalRight(final DenseMatrix64F matrix,
                                                   final double[] blockDiagonal,
                                                   final DenseMatrix64F out,
                                                   final int[] blockStarts,
                                                   final int[] blockSizes) {
        OrthogonalBlockMatrixOps.multiplyRightBlockDiagonal(
                matrix.data, blockDiagonal, matrix.numRows, blockStarts, blockSizes, out.data);
    }

    private static void rotateDBasisDiffusionGradientToOriginalBasis(final DenseMatrix64F rMatrix,
                                                                     final DenseMatrix64F dBasisGradient,
                                                                     final DenseMatrix64F temp,
                                                                     final DenseMatrix64F out) {
        OrthogonalBlockDenseMatrixOps.mult(rMatrix, dBasisGradient, temp);
        multiplyRightTransposeSymmetric(temp, rMatrix, out);
    }

    private static void multiplyRightTransposeSymmetric(final DenseMatrix64F left,
                                                        final DenseMatrix64F right,
                                                        final DenseMatrix64F out) {
        final int dimension = left.numRows;
        final double[] leftData = left.data;
        final double[] rightData = right.data;
        final double[] outData = out.data;
        for (int i = 0; i < dimension; ++i) {
            final int leftRowOffset = i * dimension;
            for (int j = i; j < dimension; ++j) {
                final int rightRowOffset = j * dimension;
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += leftData[leftRowOffset + k] * rightData[rightRowOffset + k];
                }
                outData[i * dimension + j] = sum;
                outData[j * dimension + i] = sum;
            }
        }
    }

    private static void fillDenseAdjointExpGradient(final double[] blockDParams,
                                                    final double dt,
                                                    final DenseMatrix64F upstream,
                                                    final double[] scaledNegativeBlockDScratch,
                                                    final double[] denseAdjointScratch,
                                                    final DenseMatrix64F gradD) {
        fillScaledNegativeBlockDMatrix(blockDParams, dt, scaledNegativeBlockDScratch, upstream.numRows);
        MatrixOps.toFlat(upstream, denseAdjointScratch, upstream.numRows);
        MatrixExponentialUtils.adjointExpFlat(
                scaledNegativeBlockDScratch,
                denseAdjointScratch,
                denseAdjointScratch,
                upstream.numRows);
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

    private static void fillSymmetricDenseMatrixFlat(final double[] source, final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                data[i * dimension + j] = 0.5 * (source[i * dimension + j] + source[j * dimension + i]);
            }
        }
    }

    private static void fillSymmetricDenseMatrixFlat(final double[] source,
                                                     final boolean transpose,
                                                     final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                final double ij = transpose ? source[j * dimension + i] : source[i * dimension + j];
                final double ji = transpose ? source[i * dimension + j] : source[j * dimension + i];
                data[i * dimension + j] = 0.5 * (ij + ji);
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

    private static void fillDenseMatrix(final double[] source,
                                        final DenseMatrix64F out) {
        MatrixOps.fromFlat(source, out, out.numRows);
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
