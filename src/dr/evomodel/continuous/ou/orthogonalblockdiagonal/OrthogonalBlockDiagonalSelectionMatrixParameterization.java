package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalFrechetHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovAdjointHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovSolver;
import dr.evomodel.continuous.ou.DenseSelectionMatrixParameterization;
import dr.evomodel.continuous.ou.MatrixExponentialUtils;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianUtils;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

/**
 * Orthogonal block-diagonal parametrization using block-space helpers for the
 * forward path and the orthogonal native backward path.
 */
public final class OrthogonalBlockDiagonalSelectionMatrixParameterization
        extends DenseSelectionMatrixParameterization
        implements OrthogonalBlockCanonicalParameterization {

    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final OrthogonalMatrixProvider orthogonalRotation;
    private final BlockDiagonalExpSolver expSolver;
    private final BlockDiagonalFrechetHelper frechetHelper;
    private final BlockDiagonalLyapunovSolver lyapunovSolver;
    private final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper;
    private final OrthogonalBlockBasisCache basisCache;
    private final OrthogonalBlockTransitionFactory transitionFactory;
    private final OrthogonalBlockDenseFallbackPolicy denseFallbackPolicy;
    private final DenseMatrix64F qMatrix;
    private final DenseMatrix64F qDBasis;
    private final DenseMatrix64F stationaryCovDBasis;
    private final DenseMatrix64F transitionCovDBasis;
    private final DenseMatrix64F transitionCovariance;
    private final double[] transitionMatrixArrayScratch;
    private final double[] transitionCovarianceArrayScratch;
    private final double[] precisionFlat;
    private final double[] nativeBlockGradientScratch;
    private final double[] transitionOffsetScratch;
    private final double[][] scaledNegativeBlockDScratch;
    private final double[][] denseAdjointScratch;
    private final DenseMatrix64F upstreamF;
    private final DenseMatrix64F upstreamFD;
    private final DenseMatrix64F gradD;
    private final DenseMatrix64F gradR;
    private final DenseMatrix64F gV;
    private final DenseMatrix64F hDBasis;
    private final DenseMatrix64F gS;
    private final DenseMatrix64F yAdjoint;
    private final DenseMatrix64F gECov;
    private final DenseMatrix64F temp1;
    private final DenseMatrix64F temp2;
    private final DenseMatrix64F temp3;
    private final double[] choleskyScratch;
    private final double[] lowerInverseScratch;
    private final CanonicalTransitionAdjointUtils.Workspace canonicalAdjointWorkspace;
    private final double[] tempVector1;
    private final double[] tempVector2;

    public OrthogonalBlockDiagonalSelectionMatrixParameterization(
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final OrthogonalMatrixProvider orthogonalRotation) {
        super(blockParameter);
        this.blockParameter = blockParameter;
        this.orthogonalRotation = orthogonalRotation;
        this.expSolver = new BlockDiagonalExpSolver(
                new BlockDiagonalExpSolver.BlockStructure(
                        blockParameter.getRowDimension(),
                        blockParameter.getBlockStarts(),
                        blockParameter.getBlockSizes()));
        this.basisCache = new OrthogonalBlockBasisCache(blockParameter, orthogonalRotation, expSolver);
        this.frechetHelper = new BlockDiagonalFrechetHelper(expSolver.getStructure());
        this.lyapunovSolver = new BlockDiagonalLyapunovSolver(
                getDimension(),
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes());
        this.lyapunovAdjointHelper = new BlockDiagonalLyapunovAdjointHelper(getDimension(), lyapunovSolver);
        this.transitionFactory = new OrthogonalBlockTransitionFactory(getDimension(), lyapunovSolver);
        this.denseFallbackPolicy = new OrthogonalBlockDenseFallbackPolicy();

        final int d = getDimension();
        this.qMatrix = new DenseMatrix64F(d, d);
        this.qDBasis = new DenseMatrix64F(d, d);
        this.stationaryCovDBasis = new DenseMatrix64F(d, d);
        this.transitionCovDBasis = new DenseMatrix64F(d, d);
        this.transitionCovariance = new DenseMatrix64F(d, d);
        this.transitionMatrixArrayScratch = new double[d * d];
        this.transitionCovarianceArrayScratch = new double[d * d];
        this.precisionFlat = new double[d * d];
        this.nativeBlockGradientScratch = new double[blockParameter.getBlockDiagonalNParameters()];
        this.transitionOffsetScratch = new double[d];
        this.scaledNegativeBlockDScratch = new double[d][d];
        this.denseAdjointScratch = new double[d][d];
        this.upstreamF = new DenseMatrix64F(d, d);
        this.upstreamFD = new DenseMatrix64F(d, d);
        this.gradD = new DenseMatrix64F(d, d);
        this.gradR = new DenseMatrix64F(d, d);
        this.gV = new DenseMatrix64F(d, d);
        this.hDBasis = new DenseMatrix64F(d, d);
        this.gS = new DenseMatrix64F(d, d);
        this.yAdjoint = new DenseMatrix64F(d, d);
        this.gECov = new DenseMatrix64F(d, d);
        this.temp1 = new DenseMatrix64F(d, d);
        this.temp2 = new DenseMatrix64F(d, d);
        this.temp3 = new DenseMatrix64F(d, d);
        this.choleskyScratch = new double[d * d];
        this.lowerInverseScratch = new double[d * d];
        this.canonicalAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(d);
        this.tempVector1 = new double[d];
        this.tempVector2 = new double[d];
    }

    @Override
    public void fillTransitionMatrix(final double dt, final double[][] out) {
        refreshBasisCaches(dt);
        copyDenseMatrixToArray(basisCache.transitionMatrix, out);
    }

    @Override
    public void fillTransitionMatrixFlat(final double dt, final double[] out) {
        if (out == null || out.length != getDimension() * getDimension()) {
            throw new IllegalArgumentException(
                    "transition matrix must have length " + (getDimension() * getDimension()));
        }
        refreshBasisCaches(dt);
        copyDenseMatrixToFlat(basisCache.transitionMatrix, out);
    }

    public void fillTransitionCovariance(final MatrixParameterInterface diffusionMatrix,
                                         final double dt,
                                         final double[][] out) {
        refreshBasisCaches(dt);
        transitionFactory.fillTransitionCovariance(diffusionMatrix, basisCache, out);
    }

    public void fillCanonicalTransition(final MatrixParameterInterface diffusionMatrix,
                                        final double[] stationaryMean,
                                        final double dt,
                                        final CanonicalGaussianTransition out) {
        refreshBasisCaches(dt);
        transitionFactory.fillCanonicalTransition(diffusionMatrix, stationaryMean, basisCache, out);
    }

    public OrthogonalBlockPreparedBranchBasis prepareBranchBasis(final double dt,
                                                  final double[] stationaryMean) {
        final OrthogonalBlockPreparedBranchBasis prepared = createPreparedBranchBasis();
        prepareBranchBasis(dt, stationaryMean, prepared);
        return prepared;
    }

    public OrthogonalBlockPreparedBranchBasis createPreparedBranchBasis() {
        final int dimension = getDimension();
        return new OrthogonalBlockPreparedBranchBasis(dimension, blockParameter.getTridiagonalDDimension());
    }

    public void prepareBranchBasis(final double dt,
                                   final double[] stationaryMean,
                                   final OrthogonalBlockPreparedBranchBasis prepared) {
        final int dimension = getDimension();
        if (prepared.dimension != dimension) {
            throw new IllegalArgumentException(
                    "prepared basis dimension must be " + dimension + " but is " + prepared.dimension);
        }
        OrthogonalBlockPreparedBasisBuilder.prepare(
                blockParameter,
                orthogonalRotation,
                expSolver,
                dt,
                stationaryMean,
                prepared);
    }

    public OrthogonalBlockBranchGradientWorkspace createBranchGradientWorkspace() {
        return new OrthogonalBlockBranchGradientWorkspace(
                getDimension(),
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes());
    }

    public void fillCanonicalTransitionPrepared(final OrthogonalBlockPreparedBranchBasis prepared,
                                                final MatrixParameterInterface diffusionMatrix,
                                                final OrthogonalBlockBranchGradientWorkspace workspace,
                                                final CanonicalGaussianTransition out) {
        loadOrFillPreparedCovariance(prepared, diffusionMatrix, workspace);
        fillCanonicalTransitionDirectPrepared(prepared, workspace.transitionCovariance, workspace, out);
    }

    public void prepareBranchCovariance(final OrthogonalBlockPreparedBranchBasis prepared,
                                        final MatrixParameterInterface diffusionMatrix,
                                        final OrthogonalBlockBranchGradientWorkspace workspace) {
        fillTransitionCovarianceMatrixPrepared(prepared, diffusionMatrix, workspace, workspace.transitionCovariance);
        OrthogonalBlockTransitionCovarianceSolver.storePreparedCovariance(prepared, workspace);
    }

    public void accumulateNativeGradientFromAdjointsPrepared(final OrthogonalBlockPreparedBranchBasis prepared,
                                                             final MatrixParameterInterface diffusionMatrix,
                                                             final CanonicalLocalTransitionAdjoints localAdjoints,
                                                             final OrthogonalBlockBranchGradientWorkspace workspace,
                                                             final double[] compressedDAccumulator,
                                                             final double[][] rotationAccumulator) {
        loadOrFillPreparedCovariance(prepared, diffusionMatrix, workspace);

        accumulateNativeGradientFromTransitionPrepared(
                prepared,
                localAdjoints.dLogL_dF,
                localAdjoints.dLogL_df,
                workspace,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException(
                    "Non-finite orthogonal native transition contribution at dt=" + prepared.dt);
        }

        accumulateNativeGradientFromCovarianceStationaryPrepared(
                prepared,
                localAdjoints.dLogL_dOmega,
                workspace,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException(
                    "Non-finite orthogonal native covariance contribution at dt=" + prepared.dt);
        }
    }

    public void accumulateNativeGradientFromAdjointsPreparedFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                                                 final MatrixParameterInterface diffusionMatrix,
                                                                 final CanonicalLocalTransitionAdjoints localAdjoints,
                                                                 final OrthogonalBlockBranchGradientWorkspace workspace,
                                                                 final double[] compressedDAccumulator,
                                                                 final double[] rotationAccumulator) {
        loadOrFillPreparedCovariance(prepared, diffusionMatrix, workspace);

        accumulateNativeGradientFromTransitionPreparedFlat(
                prepared,
                localAdjoints.dLogL_dF,
                localAdjoints.dLogL_df,
                workspace,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException(
                    "Non-finite orthogonal native transition contribution at dt=" + prepared.dt);
        }

        accumulateNativeGradientFromCovarianceStationaryPreparedFlat(
                prepared,
                localAdjoints.dLogL_dOmega,
                workspace,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException(
                    "Non-finite orthogonal native covariance contribution at dt=" + prepared.dt);
        }
    }

    private void loadOrFillPreparedCovariance(final OrthogonalBlockPreparedBranchBasis prepared,
                                              final MatrixParameterInterface diffusionMatrix,
                                              final OrthogonalBlockBranchGradientWorkspace workspace) {
        if (!prepared.covariancePrepared) {
            prepareBranchCovariance(prepared, diffusionMatrix, workspace);
            return;
        }
        OrthogonalBlockTransitionCovarianceSolver.copyPreparedCovariance(prepared, workspace);
    }

    public void accumulateMeanGradientPrepared(final OrthogonalBlockPreparedBranchBasis prepared,
                                               final double[] dLogL_df,
                                               final double[] gradientAccumulator,
                                               final OrthogonalBlockBranchGradientWorkspace workspace) {
        multiplyRowMajorVector(prepared.rtMatrix.data, dLogL_df, workspace.tempVector1, prepared.dimension);
        applyBlockDiagonalTransposeExpPrepared(prepared, workspace.tempVector1, workspace.tempVector2);
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

    public void accumulateDiffusionGradientPrepared(final OrthogonalBlockPreparedBranchBasis prepared,
                                                    final double[][] dLogL_dV,
                                                    final double[] gradientAccumulator,
                                                    final OrthogonalBlockBranchGradientWorkspace workspace) {
        fillSymmetricDenseMatrix(dLogL_dV, workspace.gV);
        accumulateDiffusionGradientPreparedSymmetric(prepared, gradientAccumulator, workspace);
    }

    public void accumulateDiffusionGradientPreparedFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                                        final double[] dLogL_dV,
                                                        final boolean transposeAdjoint,
                                                        final double[] gradientAccumulator,
                                                        final OrthogonalBlockBranchGradientWorkspace workspace) {
        fillSymmetricDenseMatrixFlat(dLogL_dV, workspace.gV);
        accumulateDiffusionGradientPreparedSymmetric(prepared, gradientAccumulator, workspace);
    }

    private void accumulateDiffusionGradientPreparedSymmetric(final OrthogonalBlockPreparedBranchBasis prepared,
                                                              final double[] gradientAccumulator,
                                                              final OrthogonalBlockBranchGradientWorkspace workspace) {
        CommonOps.mult(prepared.rtMatrix, workspace.gV, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.rMatrix, workspace.hDBasis);
        symmetrize(workspace.hDBasis);

        CommonOps.multTransA(prepared.expD, workspace.hDBasis, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.expD, workspace.gS);
        CommonOps.subtract(workspace.hDBasis, workspace.gS, workspace.gS);
        symmetrize(workspace.gS);

        workspace.yAdjoint.set(workspace.gS);
        CommonOps.scale(-1.0, workspace.yAdjoint);
        workspace.lyapunovAdjointHelper.solveAdjointInDBasis(
                workspace.yAdjoint, prepared.blockDParams, workspace.yAdjoint);

        CommonOps.mult(prepared.rMatrix, workspace.yAdjoint, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.rtMatrix, workspace.temp2);
        symmetrize(workspace.temp2);

        final int dimension = prepared.dimension;
        final double[] data = workspace.temp2.data;
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                gradientAccumulator[rowOffset + j] += data[rowOffset + j];
            }
        }
    }

    public void fillCanonicalLocalAdjoints(final MatrixParameterInterface diffusionMatrix,
                                           final double[] stationaryMean,
                                           final double dt,
                                           final CanonicalBranchMessageContribution contribution,
                                           final CanonicalLocalTransitionAdjoints out) {
        refreshBasisCaches(dt);
        transitionFactory.fillCanonicalLocalAdjoints(
                diffusionMatrix, stationaryMean, basisCache, contribution, out);
    }

    private void fillTransitionCovarianceMatrix(final MatrixParameterInterface diffusionMatrix,
                                                final double dt,
                                                final DenseMatrix64F out) {
        refreshBasisCaches(dt);
        OrthogonalBlockTransitionCovarianceSolver.fillTransitionCovariance(
                diffusionMatrix,
                basisCache.rMatrix,
                basisCache.rtMatrix,
                basisCache.expD,
                basisCache.blockDParams,
                lyapunovSolver,
                qMatrix,
                qDBasis,
                stationaryCovDBasis,
                transitionCovDBasis,
                temp1,
                out,
                false);
    }

    public void accumulateNativeGradientFromTransition(final double dt,
                                                       final double[] stationaryMean,
                                                       final double[][] dLogL_dF,
                                                       final double[] dLogL_df,
                                                       final double[] compressedDAccumulator,
                                                       final double[][] rotationAccumulator) {
        copySquareMatrixToFlat(dLogL_dF, transitionMatrixArrayScratch);
        accumulateNativeGradientFromTransition(
                dt, stationaryMean, transitionMatrixArrayScratch, dLogL_df,
                compressedDAccumulator, rotationAccumulator);
    }

    @Override
    public void accumulateGradientFromTransition(final double dt,
                                                 final double[] stationaryMean,
                                                 final double[][] dLogL_dF,
                                                 final double[] dLogL_df,
                                                 final double[] gradientAccumulator) {
        final int nativeBlockDim = blockParameter.getBlockDiagonalNParameters();
        final int angleDim = orthogonalRotation.getOrthogonalParameter().getDimension();
        final int nativeDim = nativeBlockDim + angleDim;
        if (gradientAccumulator.length != nativeDim || nativeDim == getDimension() * getDimension()) {
            super.accumulateGradientFromTransition(
                    dt, stationaryMean, dLogL_dF, dLogL_df, gradientAccumulator);
            return;
        }
        copySquareMatrixToFlat(dLogL_dF, transitionCovarianceArrayScratch);
        accumulateGradientFromTransitionFlat(
                dt, stationaryMean, transitionCovarianceArrayScratch, dLogL_df, gradientAccumulator);
    }

    public void accumulateNativeGradientFromTransition(final double dt,
                                                       final double[] stationaryMean,
                                                       final double[] dLogL_dF,
                                                       final double[] dLogL_df,
                                                       final double[] compressedDAccumulator,
                                                       final double[][] rotationAccumulator) {
        refreshBasisCaches(dt);
        accumulateNativeGradientFromTransitionCached(
                dt, stationaryMean, dLogL_dF, dLogL_df, compressedDAccumulator, rotationAccumulator);
    }

    @Override
    public void accumulateGradientFromTransitionFlat(final double dt,
                                                     final double[] stationaryMean,
                                                     final double[] dLogL_dF,
                                                     final double[] dLogL_df,
                                                     final double[] gradientAccumulator) {
        final int nativeBlockDim = blockParameter.getBlockDiagonalNParameters();
        final int angleDim = orthogonalRotation.getOrthogonalParameter().getDimension();
        final int nativeDim = nativeBlockDim + angleDim;
        if (gradientAccumulator.length != nativeDim || nativeDim == getDimension() * getDimension()) {
            super.accumulateGradientFromTransitionFlat(
                    dt, stationaryMean, dLogL_dF, dLogL_df, gradientAccumulator);
            return;
        }

        final int compressedDim = blockParameter.getCompressedDDimension();
        Arrays.fill(transitionMatrixArrayScratch, 0, compressedDim, 0.0);
        Arrays.fill(precisionFlat, 0, nativeBlockDim, 0.0);
        clearSquare(denseAdjointScratch);

        accumulateNativeGradientFromTransition(
                dt,
                stationaryMean,
                dLogL_dF,
                dLogL_df,
                transitionMatrixArrayScratch,
                denseAdjointScratch);

        OrthogonalBlockGradientPullback.accumulateNativeParameterGradient(
                blockParameter,
                orthogonalRotation,
                transitionMatrixArrayScratch,
                denseAdjointScratch,
                nativeBlockGradientScratch,
                gradientAccumulator);
    }

    public void accumulateNativeGradientFromCovarianceStationary(final MatrixParameterInterface diffusionMatrix,
                                                                 final double dt,
                                                                 final double[][] dLogL_dV,
                                                                 final double[] compressedDAccumulator,
                                                                 final double[][] rotationAccumulator) {
        copySquareMatrixToFlat(dLogL_dV, transitionCovarianceArrayScratch);
        accumulateNativeGradientFromCovarianceStationary(
                diffusionMatrix, dt, transitionCovarianceArrayScratch,
                compressedDAccumulator, rotationAccumulator);
    }

    public void accumulateNativeGradientFromCovarianceStationary(final MatrixParameterInterface diffusionMatrix,
                                                                 final double dt,
                                                                 final double[] dLogL_dV,
                                                                 final double[] compressedDAccumulator,
                                                                 final double[][] rotationAccumulator) {
        refreshBasisCaches(dt);
        fillDenseMatrix(diffusionMatrix, qMatrix);
        CommonOps.mult(basisCache.rtMatrix, qMatrix, temp1);
        CommonOps.mult(temp1, basisCache.rMatrix, qDBasis);
        lyapunovSolver.solve(basisCache.blockDParams, qDBasis, stationaryCovDBasis);
        CommonOps.mult(basisCache.expD, stationaryCovDBasis, temp1);
        CommonOps.multTransB(temp1, basisCache.expD, transitionCovDBasis);
        CommonOps.subtract(stationaryCovDBasis, transitionCovDBasis, transitionCovDBasis);
        symmetrize(transitionCovDBasis);
        accumulateNativeGradientFromCovarianceStationaryCached(
                dt, dLogL_dV, compressedDAccumulator, rotationAccumulator);
    }

    public void accumulateDiffusionGradient(final MatrixParameterInterface diffusionMatrix,
                                            final double dt,
                                            final double[][] dLogL_dV,
                                            final double[] gradientAccumulator) {
        refreshBasisCaches(dt);
        fillSymmetricDenseMatrix(dLogL_dV, gV);

        // Map covariance adjoint to block-D basis: G_V^D = R^T G_V R.
        CommonOps.mult(basisCache.rtMatrix, gV, temp1);
        CommonOps.mult(temp1, basisCache.rMatrix, hDBasis);
        symmetrize(hDBasis);

        // For V = R (S - E S E^T) R^T with D S + S D^T = Q_D, the exact adjoint is:
        //   G_S = G_V^D - E^T G_V^D E
        //   D^T Y + Y D = G_S
        //   G_Q^D = Y
        //   G_Q = R G_Q^D R^T
        CommonOps.multTransA(basisCache.expD, hDBasis, temp1); // E^T G_V^D
        CommonOps.mult(temp1, basisCache.expD, gS);            // E^T G_V^D E
        CommonOps.subtract(hDBasis, gS, gS);        // G_S
        symmetrize(gS);

        // Helper solves D^T Y + Y D = -hBlock, so pass hBlock = -G_S.
        yAdjoint.set(gS);
        CommonOps.scale(-1.0, yAdjoint);
        lyapunovAdjointHelper.solveAdjointInDBasis(yAdjoint, basisCache.blockDParams, yAdjoint);

        CommonOps.mult(basisCache.rMatrix, yAdjoint, temp1);
        CommonOps.mult(temp1, basisCache.rtMatrix, temp2);
        symmetrize(temp2);

        final int d = getDimension();
        final double[] data = temp2.data;
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                gradientAccumulator[rowOffset + j] += data[rowOffset + j];
            }
        }
    }

    public void accumulateNativeGradientFromCanonicalContribution(final MatrixParameterInterface diffusionMatrix,
                                                                  final double[] stationaryMean,
                                                                  final double dt,
                                                                  final CanonicalBranchMessageContribution contribution,
                                                                  final CanonicalLocalTransitionAdjoints localAdjoints,
                                                                  final double[] compressedDAccumulator,
                                                                  final double[][] rotationAccumulator) {
        refreshBasisCaches(dt);
        fillTransitionCovarianceMatrix(diffusionMatrix, dt, transitionCovariance);
        copyDenseMatrixToFlat(basisCache.transitionMatrix, transitionMatrixArrayScratch);
        fillTransitionOffset(stationaryMean, transitionOffsetScratch);

        copyAndInvertPositiveDefiniteFlat(transitionCovariance, transitionCovarianceArrayScratch, precisionFlat);
        CanonicalTransitionAdjointUtils.fillFromMoments(
                precisionFlat,
                transitionCovarianceArrayScratch,
                transitionMatrixArrayScratch,
                transitionOffsetScratch,
                contribution,
                canonicalAdjointWorkspace,
                localAdjoints);

        accumulateNativeGradientFromTransitionCached(
                dt,
                stationaryMean,
                localAdjoints.dLogL_dF,
                localAdjoints.dLogL_df,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException("Non-finite orthogonal native transition contribution at dt=" + dt);
        }
        accumulateNativeGradientFromCovarianceStationaryCached(
                dt,
                localAdjoints.dLogL_dOmega,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException("Non-finite orthogonal native covariance contribution at dt=" + dt);
        }
    }

    public void accumulateNativeGradientFromAdjoints(final MatrixParameterInterface diffusionMatrix,
                                                     final double[] stationaryMean,
                                                     final double dt,
                                                     final CanonicalLocalTransitionAdjoints localAdjoints,
                                                     final double[] compressedDAccumulator,
                                                     final double[][] rotationAccumulator) {
        refreshBasisCaches(dt);
        fillTransitionCovarianceMatrix(diffusionMatrix, dt, transitionCovariance);

        accumulateNativeGradientFromTransitionCached(
                dt,
                stationaryMean,
                localAdjoints.dLogL_dF,
                localAdjoints.dLogL_df,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException("Non-finite orthogonal native transition contribution at dt=" + dt);
        }

        accumulateNativeGradientFromCovarianceStationaryCached(
                dt,
                localAdjoints.dLogL_dOmega,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException("Non-finite orthogonal native covariance contribution at dt=" + dt);
        }
    }

    public void accumulateMeanGradient(final double dt,
                                       final double[] dLogL_df,
                                       final double[] gradientAccumulator) {
        refreshMeanGradientCaches(dt);
        orthogonalRotation.applyOrthogonalTranspose(dLogL_df, tempVector1);
        applyBlockDiagonalTransposeExp(tempVector1, tempVector2);
        orthogonalRotation.applyOrthogonal(tempVector2, tempVector1);
        for (int i = 0; i < getDimension(); ++i) {
            gradientAccumulator[i] += dLogL_df[i] - tempVector1[i];
        }
    }

    public double accumulateScalarMeanGradient(final double dt,
                                               final double[] dLogL_df) {
        refreshMeanGradientCaches(dt);
        orthogonalRotation.applyOrthogonalTranspose(dLogL_df, tempVector1);
        applyBlockDiagonalTransposeExp(tempVector1, tempVector2);
        orthogonalRotation.applyOrthogonal(tempVector2, tempVector1);
        double sum = 0.0;
        for (int i = 0; i < getDimension(); ++i) {
            sum += dLogL_df[i] - tempVector1[i];
        }
        return sum;
    }

    private void accumulateNativeGradientFromTransitionCached(final double dt,
                                                              final double[] stationaryMean,
                                                              final double[] dLogL_dF,
                                                              final double[] dLogL_df,
                                                              final double[] compressedDAccumulator,
                                                              final double[][] rotationAccumulator) {
        fillTotalUpstreamOnTransition(stationaryMean, dLogL_dF, dLogL_df, upstreamF);
        final boolean forceDenseAdjointExp = denseFallbackPolicy.forceDenseAdjointExp();

        CommonOps.mult(basisCache.rtMatrix, upstreamF, temp1);
        CommonOps.mult(temp1, basisCache.rMatrix, upstreamFD);
        if (forceDenseAdjointExp) {
            fillScaledNegativeBlockDMatrix(dt, scaledNegativeBlockDScratch);
            copyDenseMatrixToArray(upstreamFD, denseAdjointScratch);
            MatrixExponentialUtils.adjointExp(
                    scaledNegativeBlockDScratch,
                    denseAdjointScratch,
                    denseAdjointScratch);
            fillDenseMatrix(denseAdjointScratch, gradD);
            transposeInPlace(gradD);
            CommonOps.scale(-dt, gradD);
        } else {
            final DenseMatrix64F frechetInputTransition;
            if (denseFallbackPolicy.transposeNativeFrechetInput()) {
                CommonOps.transpose(upstreamFD, temp3);
                frechetInputTransition = temp3;
            } else {
                frechetInputTransition = upstreamFD;
            }
            frechetHelper.frechetAdjointExpInDBasis(basisCache.blockDParams, frechetInputTransition, dt, gradD);
            transposeInPlace(gradD);
            if (!isFinite(gradD.data)) {
                fillScaledNegativeBlockDMatrix(dt, scaledNegativeBlockDScratch);
                copyDenseMatrixToArray(upstreamFD, denseAdjointScratch);
                MatrixExponentialUtils.adjointExp(
                        scaledNegativeBlockDScratch,
                        denseAdjointScratch,
                        denseAdjointScratch);
                fillDenseMatrix(denseAdjointScratch, gradD);
                transposeInPlace(gradD);
                CommonOps.scale(-dt, gradD);
            }
        }
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        CommonOps.multTransB(basisCache.rMatrix, basisCache.expD, temp1);     // R * E^T
        CommonOps.mult(upstreamF, temp1, gradR);        // U * R * E^T
        CommonOps.multTransA(upstreamF, basisCache.rMatrix, temp1);// U^T * R
        CommonOps.mult(temp1, basisCache.expD, temp2);             // U^T * R * E
        CommonOps.addEquals(gradR, temp2);
        addDenseMatrixToArray(gradR, rotationAccumulator);
    }

    private void fillTransitionCovarianceMatrixPrepared(final OrthogonalBlockPreparedBranchBasis prepared,
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

    private void fillCanonicalTransitionDirectPrepared(final OrthogonalBlockPreparedBranchBasis prepared,
                                                       final DenseMatrix64F transitionCovariance,
                                                       final OrthogonalBlockBranchGradientWorkspace workspace,
                                                       final CanonicalGaussianTransition out) {
        OrthogonalBlockCanonicalTransitionAssembler.fillCanonicalTransition(
                prepared.transitionMatrix,
                transitionCovariance,
                prepared.stationaryMean,
                workspace.transitionOffsetScratch,
                workspace.transitionCovarianceArrayScratch,
                workspace.choleskyScratch,
                workspace.lowerInverseScratch,
                out);
    }

    private void accumulateNativeGradientFromTransitionPrepared(final OrthogonalBlockPreparedBranchBasis prepared,
                                                                final double[] dLogL_dF,
                                                                final double[] dLogL_df,
                                                                final OrthogonalBlockBranchGradientWorkspace workspace,
                                                                final double[] compressedDAccumulator,
                                                                final double[][] rotationAccumulator) {
        fillTotalUpstreamOnTransition(prepared.stationaryMean, dLogL_dF, dLogL_df, workspace.upstreamF);
        final boolean forceDenseAdjointExp = denseFallbackPolicy.forceDenseAdjointExp();

        CommonOps.mult(prepared.rtMatrix, workspace.upstreamF, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.rMatrix, workspace.upstreamFD);
        if (forceDenseAdjointExp) {
            fillScaledNegativeBlockDMatrix(
                    prepared.blockDParams,
                    prepared.dt,
                    workspace.scaledNegativeBlockDScratch,
                    prepared.dimension);
            copyDenseMatrixToArray(workspace.upstreamFD, workspace.denseAdjointScratch);
            MatrixExponentialUtils.adjointExp(
                    workspace.scaledNegativeBlockDScratch,
                    workspace.denseAdjointScratch,
                    workspace.denseAdjointScratch);
            fillDenseMatrix(workspace.denseAdjointScratch, workspace.gradD);
            transposeInPlace(workspace.gradD);
            CommonOps.scale(-prepared.dt, workspace.gradD);
        } else {
            final DenseMatrix64F frechetInputTransition;
            if (denseFallbackPolicy.transposeNativeFrechetInput()) {
                CommonOps.transpose(workspace.upstreamFD, workspace.temp3);
                frechetInputTransition = workspace.temp3;
            } else {
                frechetInputTransition = workspace.upstreamFD;
            }
            workspace.frechetHelper.frechetAdjointExpInDBasis(
                    prepared.blockDParams, frechetInputTransition, prepared.dt, workspace.gradD);
            transposeInPlace(workspace.gradD);
            if (!isFinite(workspace.gradD.data)) {
                fillScaledNegativeBlockDMatrix(
                        prepared.blockDParams,
                        prepared.dt,
                        workspace.scaledNegativeBlockDScratch,
                        prepared.dimension);
                copyDenseMatrixToArray(workspace.upstreamFD, workspace.denseAdjointScratch);
                MatrixExponentialUtils.adjointExp(
                        workspace.scaledNegativeBlockDScratch,
                        workspace.denseAdjointScratch,
                        workspace.denseAdjointScratch);
                fillDenseMatrix(workspace.denseAdjointScratch, workspace.gradD);
                transposeInPlace(workspace.gradD);
                CommonOps.scale(-prepared.dt, workspace.gradD);
            }
        }
        accumulateCompressedGradient(workspace.gradD, compressedDAccumulator);

        CommonOps.multTransB(prepared.rMatrix, prepared.expD, workspace.temp1);
        CommonOps.mult(workspace.upstreamF, workspace.temp1, workspace.gradR);
        CommonOps.multTransA(workspace.upstreamF, prepared.rMatrix, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.expD, workspace.temp2);
        CommonOps.addEquals(workspace.gradR, workspace.temp2);
        addDenseMatrixToArray(workspace.gradR, rotationAccumulator);
    }

    private void accumulateNativeGradientFromTransitionPreparedFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                                                    final double[] dLogL_dF,
                                                                    final double[] dLogL_df,
                                                                    final OrthogonalBlockBranchGradientWorkspace workspace,
                                                                    final double[] compressedDAccumulator,
                                                                    final double[] rotationAccumulator) {
        fillTotalUpstreamOnTransition(prepared.stationaryMean, dLogL_dF, dLogL_df, workspace.upstreamF);
        final boolean forceDenseAdjointExp = denseFallbackPolicy.forceDenseAdjointExp();

        CommonOps.mult(prepared.rtMatrix, workspace.upstreamF, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.rMatrix, workspace.upstreamFD);
        if (forceDenseAdjointExp) {
            fillScaledNegativeBlockDMatrix(
                    prepared.blockDParams,
                    prepared.dt,
                    workspace.scaledNegativeBlockDScratch,
                    prepared.dimension);
            copyDenseMatrixToArray(workspace.upstreamFD, workspace.denseAdjointScratch);
            MatrixExponentialUtils.adjointExp(
                    workspace.scaledNegativeBlockDScratch,
                    workspace.denseAdjointScratch,
                    workspace.denseAdjointScratch);
            fillDenseMatrix(workspace.denseAdjointScratch, workspace.gradD);
            transposeInPlace(workspace.gradD);
            CommonOps.scale(-prepared.dt, workspace.gradD);
        } else {
            final DenseMatrix64F frechetInputTransition;
            if (denseFallbackPolicy.transposeNativeFrechetInput()) {
                CommonOps.transpose(workspace.upstreamFD, workspace.temp3);
                frechetInputTransition = workspace.temp3;
            } else {
                frechetInputTransition = workspace.upstreamFD;
            }
            workspace.frechetHelper.frechetAdjointExpInDBasis(
                    prepared.blockDParams, frechetInputTransition, prepared.dt, workspace.gradD);
            transposeInPlace(workspace.gradD);
            if (!isFinite(workspace.gradD.data)) {
                fillScaledNegativeBlockDMatrix(
                        prepared.blockDParams,
                        prepared.dt,
                        workspace.scaledNegativeBlockDScratch,
                        prepared.dimension);
                copyDenseMatrixToArray(workspace.upstreamFD, workspace.denseAdjointScratch);
                MatrixExponentialUtils.adjointExp(
                        workspace.scaledNegativeBlockDScratch,
                        workspace.denseAdjointScratch,
                        workspace.denseAdjointScratch);
                fillDenseMatrix(workspace.denseAdjointScratch, workspace.gradD);
                transposeInPlace(workspace.gradD);
                CommonOps.scale(-prepared.dt, workspace.gradD);
            }
        }
        accumulateCompressedGradient(workspace.gradD, compressedDAccumulator);

        CommonOps.multTransB(prepared.rMatrix, prepared.expD, workspace.temp1);
        CommonOps.mult(workspace.upstreamF, workspace.temp1, workspace.gradR);
        CommonOps.multTransA(workspace.upstreamF, prepared.rMatrix, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.expD, workspace.temp2);
        CommonOps.addEquals(workspace.gradR, workspace.temp2);
        addDenseMatrixToFlatArray(workspace.gradR, rotationAccumulator);
    }

    private void accumulateNativeGradientFromCovarianceStationaryPrepared(final OrthogonalBlockPreparedBranchBasis prepared,
                                                                          final double[] dLogL_dV,
                                                                          final OrthogonalBlockBranchGradientWorkspace workspace,
                                                                          final double[] compressedDAccumulator,
                                                                          final double[][] rotationAccumulator) {
        final boolean forceDenseAdjointExp = denseFallbackPolicy.forceDenseAdjointExp();
        fillSymmetricDenseMatrixFlat(dLogL_dV, workspace.gV);

        CommonOps.mult(prepared.rtMatrix, workspace.gV, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.rMatrix, workspace.hDBasis);
        symmetrize(workspace.hDBasis);

        CommonOps.multTransA(prepared.expD, workspace.hDBasis, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.expD, workspace.gS);
        CommonOps.changeSign(workspace.gS);
        CommonOps.addEquals(workspace.gS, workspace.hDBasis);

        CommonOps.fill(workspace.gradD, 0.0);
        workspace.lyapunovAdjointHelper.accumulateLyapunovContributionInDBasis(
                workspace.stationaryCovDBasis, workspace.gS, prepared.blockDParams, workspace.gradD);
        accumulateCompressedGradient(workspace.gradD, compressedDAccumulator);

        workspace.lyapunovAdjointHelper.solveAdjointInDBasis(
                workspace.gS, prepared.blockDParams, workspace.yAdjoint);

        CommonOps.mult(workspace.qMatrix, prepared.rMatrix, workspace.temp1);
        CommonOps.multTransB(workspace.temp1, workspace.yAdjoint, workspace.gradR);
        CommonOps.multTransA(workspace.qMatrix, prepared.rMatrix, workspace.temp1);
        CommonOps.mult(workspace.temp1, workspace.yAdjoint, workspace.temp2);
        CommonOps.addEquals(workspace.gradR, workspace.temp2);
        CommonOps.scale(-1.0, workspace.gradR);
        addDenseMatrixToArray(workspace.gradR, rotationAccumulator);

        CommonOps.mult(workspace.hDBasis, prepared.expD, workspace.temp1);
        CommonOps.mult(workspace.temp1, workspace.stationaryCovDBasis, workspace.gECov);
        CommonOps.scale(-2.0, workspace.gECov);
        if (forceDenseAdjointExp) {
            fillScaledNegativeBlockDMatrix(
                    prepared.blockDParams,
                    prepared.dt,
                    workspace.scaledNegativeBlockDScratch,
                    prepared.dimension);
            copyDenseMatrixToArray(workspace.gECov, workspace.denseAdjointScratch);
            MatrixExponentialUtils.adjointExp(
                    workspace.scaledNegativeBlockDScratch,
                    workspace.denseAdjointScratch,
                    workspace.denseAdjointScratch);
            fillDenseMatrix(workspace.denseAdjointScratch, workspace.gradD);
            transposeInPlace(workspace.gradD);
            CommonOps.scale(-prepared.dt, workspace.gradD);
        } else {
            final DenseMatrix64F frechetInputCov;
            if (denseFallbackPolicy.transposeNativeFrechetInput()) {
                CommonOps.transpose(workspace.gECov, workspace.temp3);
                frechetInputCov = workspace.temp3;
            } else {
                frechetInputCov = workspace.gECov;
            }
            workspace.frechetHelper.frechetAdjointExpInDBasis(
                    prepared.blockDParams, frechetInputCov, prepared.dt, workspace.gradD);
            transposeInPlace(workspace.gradD);
        }
        accumulateCompressedGradient(workspace.gradD, compressedDAccumulator);

        CommonOps.mult(workspace.gV, prepared.rMatrix, workspace.temp1);
        CommonOps.mult(workspace.temp1, workspace.transitionCovDBasis, workspace.gradR);
        CommonOps.multTransA(workspace.gV, prepared.rMatrix, workspace.temp2);
        CommonOps.mult(workspace.temp2, workspace.transitionCovDBasis, workspace.temp3);
        CommonOps.addEquals(workspace.gradR, workspace.temp3);
        addDenseMatrixToArray(workspace.gradR, rotationAccumulator);
    }

    private void accumulateNativeGradientFromCovarianceStationaryPreparedFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                                                              final double[] dLogL_dV,
                                                                              final OrthogonalBlockBranchGradientWorkspace workspace,
                                                                              final double[] compressedDAccumulator,
                                                                              final double[] rotationAccumulator) {
        final boolean forceDenseAdjointExp = denseFallbackPolicy.forceDenseAdjointExp();
        fillSymmetricDenseMatrixFlat(dLogL_dV, workspace.gV);

        CommonOps.mult(prepared.rtMatrix, workspace.gV, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.rMatrix, workspace.hDBasis);
        symmetrize(workspace.hDBasis);

        CommonOps.multTransA(prepared.expD, workspace.hDBasis, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.expD, workspace.gS);
        CommonOps.changeSign(workspace.gS);
        CommonOps.addEquals(workspace.gS, workspace.hDBasis);

        CommonOps.fill(workspace.gradD, 0.0);
        workspace.lyapunovAdjointHelper.accumulateLyapunovContributionInDBasis(
                workspace.stationaryCovDBasis, workspace.gS, prepared.blockDParams, workspace.gradD);
        accumulateCompressedGradient(workspace.gradD, compressedDAccumulator);

        workspace.lyapunovAdjointHelper.solveAdjointInDBasis(
                workspace.gS, prepared.blockDParams, workspace.yAdjoint);

        CommonOps.mult(workspace.qMatrix, prepared.rMatrix, workspace.temp1);
        CommonOps.multTransB(workspace.temp1, workspace.yAdjoint, workspace.gradR);
        CommonOps.multTransA(workspace.qMatrix, prepared.rMatrix, workspace.temp1);
        CommonOps.mult(workspace.temp1, workspace.yAdjoint, workspace.temp2);
        CommonOps.addEquals(workspace.gradR, workspace.temp2);
        CommonOps.scale(-1.0, workspace.gradR);
        addDenseMatrixToFlatArray(workspace.gradR, rotationAccumulator);

        CommonOps.mult(workspace.hDBasis, prepared.expD, workspace.temp1);
        CommonOps.mult(workspace.temp1, workspace.stationaryCovDBasis, workspace.gECov);
        CommonOps.scale(-2.0, workspace.gECov);
        if (forceDenseAdjointExp) {
            fillScaledNegativeBlockDMatrix(
                    prepared.blockDParams,
                    prepared.dt,
                    workspace.scaledNegativeBlockDScratch,
                    prepared.dimension);
            copyDenseMatrixToArray(workspace.gECov, workspace.denseAdjointScratch);
            MatrixExponentialUtils.adjointExp(
                    workspace.scaledNegativeBlockDScratch,
                    workspace.denseAdjointScratch,
                    workspace.denseAdjointScratch);
            fillDenseMatrix(workspace.denseAdjointScratch, workspace.gradD);
            transposeInPlace(workspace.gradD);
            CommonOps.scale(-prepared.dt, workspace.gradD);
        } else {
            final DenseMatrix64F frechetInputCov;
            if (denseFallbackPolicy.transposeNativeFrechetInput()) {
                CommonOps.transpose(workspace.gECov, workspace.temp3);
                frechetInputCov = workspace.temp3;
            } else {
                frechetInputCov = workspace.gECov;
            }
            workspace.frechetHelper.frechetAdjointExpInDBasis(
                    prepared.blockDParams, frechetInputCov, prepared.dt, workspace.gradD);
            transposeInPlace(workspace.gradD);
        }
        accumulateCompressedGradient(workspace.gradD, compressedDAccumulator);

        CommonOps.mult(workspace.gV, prepared.rMatrix, workspace.temp1);
        CommonOps.mult(workspace.temp1, workspace.transitionCovDBasis, workspace.gradR);
        CommonOps.multTransA(workspace.gV, prepared.rMatrix, workspace.temp2);
        CommonOps.mult(workspace.temp2, workspace.transitionCovDBasis, workspace.temp3);
        CommonOps.addEquals(workspace.gradR, workspace.temp3);
        addDenseMatrixToFlatArray(workspace.gradR, rotationAccumulator);
    }

    private void accumulateNativeGradientFromCovarianceStationaryCached(final double dt,
                                                                        final double[] dLogL_dV,
                                                                        final double[] compressedDAccumulator,
                                                                        final double[][] rotationAccumulator) {
        final boolean forceDenseAdjointExp = denseFallbackPolicy.forceDenseAdjointExp();
        fillSymmetricDenseMatrixFlat(dLogL_dV, gV);

        CommonOps.mult(basisCache.rtMatrix, gV, temp1);
        CommonOps.mult(temp1, basisCache.rMatrix, hDBasis);
        symmetrize(hDBasis);

        CommonOps.multTransA(basisCache.expD, hDBasis, temp1);     // E^T H
        CommonOps.mult(temp1, basisCache.expD, gS);                // E^T H E
        CommonOps.changeSign(gS);
        CommonOps.addEquals(gS, hDBasis);               // H - E^T H E

        CommonOps.fill(gradD, 0.0);
        lyapunovAdjointHelper.accumulateLyapunovContributionInDBasis(
                stationaryCovDBasis, gS, basisCache.blockDParams, gradD);
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        lyapunovAdjointHelper.solveAdjointInDBasis(gS, basisCache.blockDParams, yAdjoint);

        // Q-basis contribution: helper returns yAdjoint = -Y, where Y solves
        // D^T Y + Y D = G_S.  Since Q_D = R^T Q R, the raw R gradient is
        // Q R Y^T + Q^T R Y = -(Q R yAdjoint^T + Q^T R yAdjoint).
        CommonOps.mult(qMatrix, basisCache.rMatrix, temp1);
        CommonOps.multTransB(temp1, yAdjoint, gradR);   // Q R y^T
        CommonOps.multTransA(qMatrix, basisCache.rMatrix, temp1);  // Q^T R
        CommonOps.mult(temp1, yAdjoint, temp2);         // Q^T R y
        CommonOps.addEquals(gradR, temp2);
        CommonOps.scale(-1.0, gradR);
        addDenseMatrixToArray(gradR, rotationAccumulator);

        // E-inside-V contribution.
        CommonOps.mult(hDBasis, basisCache.expD, temp1);
        CommonOps.mult(temp1, stationaryCovDBasis, gECov);
        CommonOps.scale(-2.0, gECov);
        if (forceDenseAdjointExp) {
            fillScaledNegativeBlockDMatrix(dt, scaledNegativeBlockDScratch);
            copyDenseMatrixToArray(gECov, denseAdjointScratch);
            MatrixExponentialUtils.adjointExp(
                    scaledNegativeBlockDScratch,
                    denseAdjointScratch,
                    denseAdjointScratch);
            fillDenseMatrix(denseAdjointScratch, gradD);
            transposeInPlace(gradD);
            CommonOps.scale(-dt, gradD);
        } else {
            final DenseMatrix64F frechetInputCov;
            if (denseFallbackPolicy.transposeNativeFrechetInput()) {
                CommonOps.transpose(gECov, temp3);
                frechetInputCov = temp3;
            } else {
                frechetInputCov = gECov;
            }
            frechetHelper.frechetAdjointExpInDBasis(basisCache.blockDParams, frechetInputCov, dt, gradD);
            transposeInPlace(gradD);
        }
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        // Outer basis reconstruction V = R V_D R^T.
        CommonOps.mult(gV, basisCache.rMatrix, temp1);
        CommonOps.mult(temp1, transitionCovDBasis, gradR);
        CommonOps.multTransA(gV, basisCache.rMatrix, temp2);
        CommonOps.mult(temp2, transitionCovDBasis, temp3);
        CommonOps.addEquals(gradR, temp3);
        addDenseMatrixToArray(gradR, rotationAccumulator);
    }

    private void fillTransitionOffset(final double[] stationaryMean,
                                      final double[] out) {
        final int dimension = getDimension();
        final double[] transitionData = basisCache.transitionMatrix.data;
        for (int i = 0; i < dimension; ++i) {
            double transformedMean = 0.0;
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                transformedMean += transitionData[rowOffset + j] * stationaryMean[j];
            }
            out[i] = stationaryMean[i] - transformedMean;
        }
    }

    private void refreshBasisCaches(final double dt) {
        basisCache.refresh(dt);
    }

    private void refreshMeanGradientCaches(final double dt) {
        basisCache.refreshExpOnly(dt);
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

    private void fillTotalUpstreamOnTransition(final double[] stationaryMean,
                                               final double[] dLogL_dF,
                                               final double[] dLogL_df,
                                               final DenseMatrix64F out) {
        final int dimension = getDimension();
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

    private void applyBlockDiagonalTransposeExp(final double[] in,
                                                final double[] out) {
        for (int b = 0; b < blockParameter.getNumBlocks(); ++b) {
            final int start = blockParameter.getBlockStarts()[b];
            final int size = blockParameter.getBlockSizes()[b];
            if (size == 1) {
                out[start] = basisCache.expD.data[start * getDimension() + start] * in[start];
            } else {
                final int row0 = start * getDimension();
                final int row1 = (start + 1) * getDimension();
                final double e00 = basisCache.expD.data[row0 + start];
                final double e01 = basisCache.expD.data[row0 + start + 1];
                final double e10 = basisCache.expD.data[row1 + start];
                final double e11 = basisCache.expD.data[row1 + start + 1];
                out[start] = e00 * in[start] + e10 * in[start + 1];
                out[start + 1] = e01 * in[start] + e11 * in[start + 1];
            }
        }
    }

    private void applyBlockDiagonalTransposeExpPrepared(final OrthogonalBlockPreparedBranchBasis prepared,
                                                        final double[] in,
                                                        final double[] out) {
        final double[] expData = prepared.expD.data;
        for (int b = 0; b < blockParameter.getNumBlocks(); ++b) {
            final int start = blockParameter.getBlockStarts()[b];
            final int size = blockParameter.getBlockSizes()[b];
            if (size == 1) {
                out[start] = expData[start * prepared.dimension + start] * in[start];
            } else {
                final int row0 = start * prepared.dimension;
                final int row1 = (start + 1) * prepared.dimension;
                final double e00 = expData[row0 + start];
                final double e01 = expData[row0 + start + 1];
                final double e10 = expData[row1 + start];
                final double e11 = expData[row1 + start + 1];
                out[start] = e00 * in[start] + e10 * in[start + 1];
                out[start + 1] = e01 * in[start] + e11 * in[start + 1];
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

    private static void copyDenseMatrixToArray(final DenseMatrix64F source, final double[][] out) {
        final int dimension = source.numRows;
        final double[] data = source.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i][j] = data[i * dimension + j];
            }
        }
    }

    private static void copyDenseMatrixToFlat(final DenseMatrix64F source, final double[] out) {
        System.arraycopy(source.data, 0, out, 0, source.numRows * source.numCols);
    }

    private static void copySquareMatrixToFlat(final double[][] source, final double[] out) {
        final int dimension = source.length;
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(source[i], 0, out, i * dimension, dimension);
        }
    }

    private static void clearSquare(final double[][] matrix) {
        for (double[] row : matrix) {
            Arrays.fill(row, 0.0);
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

    private void fillScaledNegativeBlockDMatrix(final double dt,
                                                final double[][] out) {
        final int d = getDimension();
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                out[i][j] = 0.0;
            }
        }
        final int upperOffset = d;
        final int lowerOffset = d + (d - 1);
        for (int i = 0; i < d; ++i) {
            out[i][i] = -dt * basisCache.blockDParams[i];
            if (i < d - 1) {
                out[i][i + 1] = -dt * basisCache.blockDParams[upperOffset + i];
                out[i + 1][i] = -dt * basisCache.blockDParams[lowerOffset + i];
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

    private double copyAndInvertPositiveDefiniteFlat(final DenseMatrix64F source,
                                                     final double[] matrixOut,
                                                     final double[] inverseOut) {
        return OrthogonalBlockPositiveDefiniteInverter.copyAndInvertFlat(
                source, matrixOut, inverseOut, choleskyScratch, lowerInverseScratch);
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

    private static boolean isFinite(final double[] values) {
        for (double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    private static boolean isFinite(final double[][] values) {
        for (double[] row : values) {
            if (!isFinite(row)) {
                return false;
            }
        }
        return true;
    }

}
