package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.continuous.ou.DenseSelectionMatrixParameterization;
import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalFrechetHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovAdjointHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovSolver;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.TransposedMatrixParameter;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

/**
 * Block-diagonal parametrization for {@code A = R D R^{-1}} with a reusable
 * forward/canonical transition path.
 */
public class BlockDiagonalSelectionMatrixParameterization
        extends DenseSelectionMatrixParameterization
        implements BlockDiagonalNativeCanonicalParameterization {

    protected final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    protected final BlockDiagonalExpSolver expSolver;
    protected final BlockDiagonalLyapunovSolver lyapunovSolver;
    protected final BlockDiagonalTransitionCache basisCache;
    protected final BlockDiagonalTransitionFactory transitionFactory;
    protected final DenseMatrix64F qMatrix;
    protected final DenseMatrix64F qDBasis;
    protected final DenseMatrix64F stationaryCovDBasis;
    protected final DenseMatrix64F transitionCovDBasis;
    protected final DenseMatrix64F transitionCovariance;
    protected final DenseMatrix64F tempMatrix;
    protected final double[] transitionMatrixArrayScratch;
    protected final double[] transitionCovarianceArrayScratch;
    protected final double[] compressedDGradientScratch;
    protected final double[] precisionFlat;
    protected final double[] transitionOffsetScratch;
    protected final double[] denseAdjointScratch;
    protected final double[] nativeBlockGradientScratch;
    protected final double[] choleskyScratch;
    protected final double[] lowerInverseScratch;
    protected final CanonicalTransitionAdjointUtils.Workspace canonicalAdjointWorkspace;
    private final BlockDiagonalBranchGradientWorkspace currentGradientWorkspace;

    public BlockDiagonalSelectionMatrixParameterization(
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter) {
        super(blockParameter);
        this.blockParameter = blockParameter;
        this.expSolver = new BlockDiagonalExpSolver(
                new BlockDiagonalExpSolver.BlockStructure(
                        blockParameter.getRowDimension(),
                        blockParameter.getBlockStarts(),
                        blockParameter.getBlockSizes()));
        this.basisCache = new BlockDiagonalTransitionCache(blockParameter, expSolver);
        this.lyapunovSolver = new BlockDiagonalLyapunovSolver(
                getDimension(),
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes());
        this.transitionFactory = new BlockDiagonalTransitionFactory(getDimension(), lyapunovSolver);

        final int d = getDimension();
        this.qMatrix = new DenseMatrix64F(d, d);
        this.qDBasis = new DenseMatrix64F(d, d);
        this.stationaryCovDBasis = new DenseMatrix64F(d, d);
        this.transitionCovDBasis = new DenseMatrix64F(d, d);
        this.transitionCovariance = new DenseMatrix64F(d, d);
        this.tempMatrix = new DenseMatrix64F(d, d);
        this.transitionMatrixArrayScratch = new double[d * d];
        this.transitionCovarianceArrayScratch = new double[d * d];
        this.compressedDGradientScratch = new double[blockParameter.getCompressedDDimension()];
        this.precisionFlat = new double[d * d];
        this.transitionOffsetScratch = new double[d];
        this.denseAdjointScratch = new double[d * d];
        this.nativeBlockGradientScratch = new double[blockParameter.getBlockDiagonalNParameters()];
        this.choleskyScratch = new double[d * d];
        this.lowerInverseScratch = new double[d * d];
        this.canonicalAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(d);
        this.currentGradientWorkspace = new BlockDiagonalBranchGradientWorkspace(
                d,
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes());
    }

    public AbstractBlockDiagonalTwoByTwoMatrixParameter getBlockParameter() {
        return blockParameter;
    }

    @Override
    public int getSelectionGradientDimension() {
        return blockParameter.getBlockDiagonalNParameters() + getDimension() * getDimension();
    }

    @Override
    public int getCompressedSelectionGradientDimension() {
        return blockParameter.getCompressedDDimension();
    }

    @Override
    public int getNativeSelectionGradientScratchDimension() {
        return blockParameter.getBlockDiagonalNParameters();
    }

    @Override
    public void finishNativeSelectionGradient(final double[] compressedGradient,
                                             final double[] nativeGradientScratch,
                                             final double[] rotationGradientFlat,
                                             final double[] gradientOut) {
        final int nativeBlockDim = blockParameter.getBlockDiagonalNParameters();
        if (gradientOut.length != getSelectionGradientDimension()) {
            throw new IllegalArgumentException(
                    "selection gradient must have length " + getSelectionGradientDimension()
                            + " but has " + gradientOut.length);
        }
        blockParameter.chainGradient(compressedGradient, nativeGradientScratch);
        System.arraycopy(nativeGradientScratch, 0, gradientOut, 0, nativeBlockDim);
        flattenRotationGradientForParameter(rotationGradientFlat, gradientOut, nativeBlockDim, getDimension());
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
        MatrixOps.toFlat(basisCache.transitionMatrix, out, getDimension());
    }

    @Override
    public void accumulateGradientFromTransition(final double dt,
                                                 final double[] stationaryMean,
                                                 final double[][] dLogL_dF,
                                                 final double[] dLogL_df,
                                                 final double[] gradientAccumulator) {
        final int dimension = getDimension();
        for (int row = 0; row < dimension; ++row) {
            System.arraycopy(dLogL_dF[row], 0, denseAdjointScratch, row * dimension, dimension);
        }
        accumulateGradientFromTransitionFlat(
                dt, stationaryMean, denseAdjointScratch, dLogL_df, gradientAccumulator);
    }

    @Override
    public void accumulateGradientFromTransitionFlat(final double dt,
                                                     final double[] stationaryMean,
                                                     final double[] dLogL_dF,
                                                     final double[] dLogL_df,
                                                     final double[] gradientAccumulator) {
        final int nativeDim = getSelectionGradientDimension();
        if (gradientAccumulator.length != nativeDim || nativeDim == getDimension() * getDimension()) {
            super.accumulateGradientFromTransitionFlat(
                    dt, stationaryMean, dLogL_dF, dLogL_df, gradientAccumulator);
            return;
        }

        Arrays.fill(compressedDGradientScratch, 0.0);
        Arrays.fill(denseAdjointScratch, 0.0);
        Arrays.fill(nativeBlockGradientScratch, 0.0);

        accumulateNativeGradientFromTransitionFlat(
                dt,
                stationaryMean,
                dLogL_dF,
                dLogL_df,
                compressedDGradientScratch,
                denseAdjointScratch);

        blockParameter.chainGradient(compressedDGradientScratch, nativeBlockGradientScratch);
        final int nativeBlockDim = blockParameter.getBlockDiagonalNParameters();
        for (int i = 0; i < nativeBlockDim; ++i) {
            gradientAccumulator[i] += nativeBlockGradientScratch[i];
        }
        accumulateRotationGradientForParameter(
                denseAdjointScratch,
                gradientAccumulator,
                nativeBlockDim,
                getDimension());
    }

    @Override
    public void fillTransitionCovarianceFlat(final MatrixParameterInterface diffusionMatrix,
                                             final double dt,
                                             final double[] out) {
        refreshBasisCaches(dt);
        transitionFactory.fillTransitionCovarianceFlat(diffusionMatrix, basisCache, out);
    }

    @Override
    public void fillCanonicalTransition(final MatrixParameterInterface diffusionMatrix,
                                        final double[] stationaryMean,
                                        final double dt,
                                        final CanonicalGaussianTransition out) {
        refreshBasisCaches(dt);
        transitionFactory.fillCanonicalTransition(diffusionMatrix, stationaryMean, basisCache, out);
    }

    public BlockDiagonalPreparedBranchBasis prepareBranchBasis(final double dt,
                                                               final double[] stationaryMean) {
        final BlockDiagonalPreparedBranchBasis prepared = createPreparedBranchBasis();
        prepareBranchBasis(dt, stationaryMean, prepared);
        return prepared;
    }

    @Override
    public BlockDiagonalPreparedBranchBasis createPreparedBranchBasis() {
        final int dimension = getDimension();
        return new BlockDiagonalPreparedBranchBasis(
                dimension,
                blockParameter.getCompressedDDimension(),
                blockParameter.createBlockDiagonalDecomposition());
    }

    @Override
    public CanonicalPreparedBranchHandle createPreparedBranchHandle() {
        return new BlockDiagonalPreparedBranchHandle(createPreparedBranchBasis());
    }

    @Override
    public CanonicalBranchWorkspace createBranchWorkspace() {
        return createBranchGradientWorkspace();
    }

    @Override
    public BlockDiagonalBranchGradientWorkspace createBranchGradientWorkspace() {
        return new BlockDiagonalBranchGradientWorkspace(
                getDimension(),
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes());
    }

    @Override
    public boolean supportsDelayedDiffusionGradientRotation() {
        return true;
    }

    @Override
    public void prepareBranch(final double dt,
                              final double[] stationaryMean,
                              final CanonicalPreparedBranchHandle prepared) {
        prepareBranchBasis(dt, stationaryMean, asPreparedBasis(prepared));
    }

    @Override
    public void prepareBranchBasis(final double dt,
                                   final double[] stationaryMean,
                                   final BlockDiagonalPreparedBranchBasis prepared) {
        final int dimension = getDimension();
        if (prepared.dimension != dimension) {
            throw new IllegalArgumentException(
                    "prepared basis dimension must be " + dimension + " but is " + prepared.dimension);
        }
        fillPreparedBranchBasis(dt, stationaryMean, prepared);
    }

    @Override
    public void fillCanonicalTransitionPrepared(final BlockDiagonalPreparedBranchBasis prepared,
                                                final MatrixParameterInterface diffusionMatrix,
                                                final BlockDiagonalBranchGradientWorkspace workspace,
                                                final CanonicalGaussianTransition out) {
        loadOrFillPreparedCovariance(prepared, diffusionMatrix, workspace);
        fillCanonicalTransitionDirectPrepared(prepared, workspace.transitionCovariance, workspace, out);
    }

    @Override
    public void fillCanonicalTransitionPrepared(final CanonicalPreparedBranchHandle prepared,
                                                final MatrixParameterInterface diffusionMatrix,
                                                final CanonicalBranchWorkspace workspace,
                                                final CanonicalGaussianTransition out) {
        fillCanonicalTransitionPrepared(
                asPreparedBasis(prepared),
                diffusionMatrix,
                asBranchWorkspace(workspace),
                out);
    }

    @Override
    public void fillTransitionMatrixPreparedFlat(final CanonicalPreparedBranchHandle prepared,
                                                 final CanonicalBranchWorkspace workspace,
                                                 final double[] out) {
        MatrixOps.toFlat(asPreparedBasis(prepared).transitionMatrix, out, getDimension());
    }

    @Override
    public boolean fillTransitionMomentsPreparedFlat(final CanonicalPreparedBranchHandle prepared,
                                                     final MatrixParameterInterface diffusionMatrix,
                                                     final CanonicalBranchWorkspace workspace,
                                                     final double[] transitionMatrixOut,
                                                     final double[] transitionOffsetOut,
                                                     final double[] transitionCovarianceOut) {
        final BlockDiagonalPreparedBranchBasis basis = asPreparedBasis(prepared);
        final BlockDiagonalBranchGradientWorkspace branchWorkspace = asBranchWorkspace(workspace);
        loadOrFillPreparedCovariance(basis, diffusionMatrix, branchWorkspace);
        MatrixOps.toFlat(basis.transitionMatrix, transitionMatrixOut, getDimension());
        fillTransitionOffsetPrepared(basis, transitionOffsetOut);
        MatrixOps.toFlat(basis.transitionCovariance, transitionCovarianceOut, getDimension());
        return true;
    }

    public void prepareBranchCovariance(final BlockDiagonalPreparedBranchBasis prepared,
                                        final MatrixParameterInterface diffusionMatrix,
                                        final BlockDiagonalBranchGradientWorkspace workspace) {
        fillTransitionCovarianceMatrixPrepared(prepared, diffusionMatrix, workspace, workspace.transitionCovariance);
        BlockDiagonalTransitionCovarianceSolver.storePreparedCovariance(prepared, workspace);
    }

    @Override
    public void fillCanonicalLocalAdjoints(final MatrixParameterInterface diffusionMatrix,
                                           final double[] stationaryMean,
                                           final double dt,
                                           final CanonicalBranchMessageContribution contribution,
                                           final CanonicalLocalTransitionAdjoints out) {
        refreshBasisCaches(dt);
        transitionFactory.fillCanonicalLocalAdjoints(
                diffusionMatrix, stationaryMean, basisCache, contribution, out);
    }

    public void accumulateNativeGradientFromAdjointsPreparedFlat(
            final CanonicalPreparedBranchHandle prepared,
            final MatrixParameterInterface diffusionMatrix,
            final CanonicalLocalTransitionAdjoints localAdjoints,
            final CanonicalBranchWorkspace workspace,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator) {
        accumulateNativeGradientFromAdjointsPreparedFlat(
                asPreparedBasis(prepared),
                diffusionMatrix,
                localAdjoints,
                asBranchWorkspace(workspace),
                compressedDAccumulator,
                rotationAccumulator);
    }

    public void accumulateNativeSelectionAndDiffusionGradientFromAdjointsPreparedFlat(
            final CanonicalPreparedBranchHandle prepared,
            final MatrixParameterInterface diffusionMatrix,
            final CanonicalLocalTransitionAdjoints localAdjoints,
            final CanonicalBranchWorkspace workspace,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator,
            final boolean delayDiffusionGradientRotation,
            final double[] diffusionGradientAccumulator) {
        accumulateNativeSelectionAndDiffusionGradientFromAdjointsPreparedFlat(
                asPreparedBasis(prepared),
                diffusionMatrix,
                localAdjoints,
                asBranchWorkspace(workspace),
                compressedDAccumulator,
                rotationAccumulator,
                delayDiffusionGradientRotation,
                diffusionGradientAccumulator);
    }

    public void accumulateMeanGradient(final double dt,
                                       final double[] dLogL_df,
                                       final double[] gradientAccumulator) {
        refreshBasisCaches(dt);
        accumulateMeanGradientFromTransitionData(
                basisCache.transitionMatrix.data, dLogL_df, gradientAccumulator, getDimension());
    }

    public double accumulateScalarMeanGradient(final double dt,
                                               final double[] dLogL_df) {
        refreshBasisCaches(dt);
        return scalarMeanGradientFromTransitionData(
                basisCache.transitionMatrix.data, dLogL_df, getDimension());
    }

    public void accumulateMeanGradientPrepared(final CanonicalPreparedBranchHandle prepared,
                                               final double[] dLogL_df,
                                               final double[] gradientAccumulator,
                                               final CanonicalBranchWorkspace workspace) {
        accumulateMeanGradientPrepared(
                asPreparedBasis(prepared),
                dLogL_df,
                gradientAccumulator,
                asBranchWorkspace(workspace));
    }

    public void accumulateMeanGradientPreparedDBasisFlat(final CanonicalPreparedBranchHandle prepared,
                                                         final double[] dLogL_df,
                                                         final double[] dBasisGradientAccumulator,
                                                         final CanonicalBranchWorkspace workspace) {
        accumulateMeanGradientPreparedDBasisFlat(
                asPreparedBasis(prepared),
                dLogL_df,
                dBasisGradientAccumulator,
                asBranchWorkspace(workspace));
    }

    public void finishMeanGradientFromDBasisFlat(final double[] dBasisGradientAccumulator,
                                                 final double[] gradientAccumulator,
                                                 final CanonicalBranchWorkspace workspace) {
        finishMeanGradientFromDBasisFlat(
                dBasisGradientAccumulator,
                gradientAccumulator,
                asBranchWorkspace(workspace));
    }

    public void accumulateDiffusionGradientPreparedFlat(final CanonicalPreparedBranchHandle prepared,
                                                        final double[] dLogL_dV,
                                                        final boolean transposeAdjoint,
                                                        final double[] gradientAccumulator,
                                                        final CanonicalBranchWorkspace workspace) {
        accumulateDiffusionGradientPreparedFlat(
                asPreparedBasis(prepared),
                dLogL_dV,
                transposeAdjoint,
                gradientAccumulator,
                asBranchWorkspace(workspace));
    }

    public void accumulateDiffusionGradientPreparedDBasisFlat(final CanonicalPreparedBranchHandle prepared,
                                                              final double[] dLogL_dV,
                                                              final boolean transposeAdjoint,
                                                              final double[] dBasisGradientAccumulator,
                                                              final CanonicalBranchWorkspace workspace) {
        accumulateDiffusionGradientPreparedDBasisFlat(
                asPreparedBasis(prepared),
                dLogL_dV,
                transposeAdjoint,
                dBasisGradientAccumulator,
                asBranchWorkspace(workspace));
    }

    public void finishDiffusionGradientFromDBasisFlat(final double[] dBasisGradientAccumulator,
                                                      final double[] gradientAccumulator,
                                                      final CanonicalBranchWorkspace workspace) {
        finishDiffusionGradientFromDBasisFlat(
                dBasisGradientAccumulator,
                gradientAccumulator,
                asBranchWorkspace(workspace));
    }

    public void accumulateNativeGradientFromTransitionFlat(final double dt,
                                                           final double[] stationaryMean,
                                                           final double[] dLogL_dF,
                                                           final double[] dLogL_df,
                                                           final double[] compressedDAccumulator,
                                                           final double[] rotationAccumulator) {
        refreshBasisCaches(dt);
        accumulateNativeTransitionCurrentPullback(
                dt,
                stationaryMean,
                dLogL_dF,
                dLogL_df,
                compressedDAccumulator,
                rotationAccumulator);
    }

    public void accumulateNativeGradientFromCovarianceStationaryFlat(
            final MatrixParameterInterface diffusionMatrix,
            final double dt,
            final double[] dLogL_dV,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator) {
        refreshBasisCaches(dt);
        fillNativeTransitionCovarianceMatrix(diffusionMatrix, transitionCovariance);
        accumulateNativeCovarianceCurrentCachedPullback(
                diffusionMatrix,
                dt,
                dLogL_dV,
                compressedDAccumulator,
                rotationAccumulator);
    }

    public void accumulateDiffusionGradientFlat(final MatrixParameterInterface diffusionMatrix,
                                                final double dt,
                                                final double[] dLogL_dV,
                                                final boolean transposeAdjoint,
                                                final double[] gradientAccumulator) {
        refreshBasisCaches(dt);
        accumulateDiffusionGradientCurrentFlat(
                dLogL_dV, transposeAdjoint, gradientAccumulator);
    }

    public void accumulateNativeGradientFromCanonicalContributionFlat(
            final MatrixParameterInterface diffusionMatrix,
            final double[] stationaryMean,
            final double dt,
            final CanonicalBranchMessageContribution contribution,
            final CanonicalLocalTransitionAdjoints localAdjoints,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator) {
        refreshBasisCaches(dt);
        fillNativeTransitionCovarianceMatrix(diffusionMatrix, transitionCovariance);
        fillTransitionOffset(stationaryMean, transitionOffsetScratch);

        copyAndInvertPositiveDefiniteFlat(transitionCovariance, transitionCovarianceArrayScratch, precisionFlat);
        CanonicalTransitionAdjointUtils.fillFromMoments(
                precisionFlat,
                transitionCovarianceArrayScratch,
                basisCache.transitionMatrix.data,
                transitionOffsetScratch,
                contribution,
                canonicalAdjointWorkspace,
                localAdjoints);

        accumulateCurrentNativeGradientFromLocalAdjoints(
                diffusionMatrix,
                stationaryMean,
                dt,
                localAdjoints,
                compressedDAccumulator,
                rotationAccumulator);
    }

    public void accumulateNativeGradientFromAdjointsFlat(final MatrixParameterInterface diffusionMatrix,
                                                         final double[] stationaryMean,
                                                         final double dt,
                                                         final CanonicalLocalTransitionAdjoints localAdjoints,
                                                         final double[] compressedDAccumulator,
                                                         final double[] rotationAccumulator) {
        refreshBasisCaches(dt);
        fillNativeTransitionCovarianceMatrix(diffusionMatrix, transitionCovariance);
        accumulateCurrentNativeGradientFromLocalAdjoints(
                diffusionMatrix,
                stationaryMean,
                dt,
                localAdjoints,
                compressedDAccumulator,
                rotationAccumulator);
    }

    protected void loadOrFillPreparedCovariance(final BlockDiagonalPreparedBranchBasis prepared,
                                                final MatrixParameterInterface diffusionMatrix,
                                                final BlockDiagonalBranchGradientWorkspace workspace) {
        if (!prepared.covariancePrepared) {
            prepareBranchCovariance(prepared, diffusionMatrix, workspace);
            return;
        }
        BlockDiagonalTransitionCovarianceSolver.copyPreparedCovariance(prepared, workspace);
    }

    protected void accumulateNativeGradientFromAdjointsPreparedFlat(
            final BlockDiagonalPreparedBranchBasis prepared,
            final MatrixParameterInterface diffusionMatrix,
            final CanonicalLocalTransitionAdjoints localAdjoints,
            final BlockDiagonalBranchGradientWorkspace workspace,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator) {
        loadOrFillPreparedCovariance(prepared, diffusionMatrix, workspace);

        accumulateNativeTransitionPreparedPullback(
                prepared,
                localAdjoints.dLogL_dF,
                localAdjoints.dLogL_df,
                workspace,
                compressedDAccumulator,
                rotationAccumulator);
        checkNativeGradientFinite("transition", prepared.dt, compressedDAccumulator, rotationAccumulator);

        accumulateNativeCovariancePreparedPullback(
                prepared,
                diffusionMatrix,
                localAdjoints.dLogL_dOmega,
                workspace,
                compressedDAccumulator,
                rotationAccumulator);
        checkNativeGradientFinite("covariance", prepared.dt, compressedDAccumulator, rotationAccumulator);
    }

    protected void accumulateNativeSelectionAndDiffusionGradientFromAdjointsPreparedFlat(
            final BlockDiagonalPreparedBranchBasis prepared,
            final MatrixParameterInterface diffusionMatrix,
            final CanonicalLocalTransitionAdjoints localAdjoints,
            final BlockDiagonalBranchGradientWorkspace workspace,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator,
            final boolean delayDiffusionGradientRotation,
            final double[] diffusionGradientAccumulator) {
        loadOrFillPreparedCovariance(prepared, diffusionMatrix, workspace);

        accumulateNativeTransitionPreparedPullback(
                prepared,
                localAdjoints.dLogL_dF,
                localAdjoints.dLogL_df,
                workspace,
                compressedDAccumulator,
                rotationAccumulator);
        checkNativeGradientFinite("transition", prepared.dt, compressedDAccumulator, rotationAccumulator);

        accumulatePreparedCovarianceAndDiffusionGradient(
                prepared,
                diffusionMatrix,
                localAdjoints.dLogL_dOmega,
                workspace,
                compressedDAccumulator,
                rotationAccumulator,
                delayDiffusionGradientRotation,
                diffusionGradientAccumulator);
        checkNativeGradientFinite("covariance", prepared.dt, compressedDAccumulator, rotationAccumulator);
    }

    protected void accumulateMeanGradientPrepared(final BlockDiagonalPreparedBranchBasis prepared,
                                                  final double[] dLogL_df,
                                                  final double[] gradientAccumulator,
                                                  final BlockDiagonalBranchGradientWorkspace workspace) {
        accumulateMeanGradientFromTransitionData(
                prepared.transitionMatrix.data, dLogL_df, gradientAccumulator, prepared.dimension);
    }

    protected void accumulateMeanGradientPreparedDBasisFlat(final BlockDiagonalPreparedBranchBasis prepared,
                                                            final double[] dLogL_df,
                                                            final double[] dBasisGradientAccumulator,
                                                            final BlockDiagonalBranchGradientWorkspace workspace) {
        throw unsupportedNativeGradient("delayed stationary-mean gradient rotation");
    }

    protected void finishMeanGradientFromDBasisFlat(final double[] dBasisGradientAccumulator,
                                                    final double[] gradientAccumulator,
                                                    final BlockDiagonalBranchGradientWorkspace workspace) {
        throw unsupportedNativeGradient("delayed stationary-mean gradient rotation");
    }

    protected void accumulateDiffusionGradientPreparedFlat(final BlockDiagonalPreparedBranchBasis prepared,
                                                           final double[] dLogL_dV,
                                                           final boolean transposeAdjoint,
                                                           final double[] gradientAccumulator,
                                                           final BlockDiagonalBranchGradientWorkspace workspace) {
        fillSymmetricSandwichOutputAdjoint(dLogL_dV, transposeAdjoint, workspace.gV);
        fillDiffusionGradientDBasisPrepared(prepared, workspace);
        rotateDBasisDiffusionGradientToOriginalBasis(
                prepared.rinvMatrix,
                workspace.yAdjoint,
                workspace.temp1,
                workspace.temp2);
        addDenseMatrixToFlatArray(workspace.temp2, gradientAccumulator);
    }

    protected void accumulateDiffusionGradientPreparedDBasisFlat(final BlockDiagonalPreparedBranchBasis prepared,
                                                                 final double[] dLogL_dV,
                                                                 final boolean transposeAdjoint,
                                                                 final double[] dBasisGradientAccumulator,
                                                                 final BlockDiagonalBranchGradientWorkspace workspace) {
        fillSymmetricSandwichOutputAdjoint(dLogL_dV, transposeAdjoint, workspace.gV);
        fillDiffusionGradientDBasisPrepared(prepared, workspace);
        addDenseMatrixToFlatArray(workspace.yAdjoint, dBasisGradientAccumulator);
    }

    protected void finishDiffusionGradientFromDBasisFlat(final double[] dBasisGradientAccumulator,
                                                         final double[] gradientAccumulator,
                                                         final BlockDiagonalBranchGradientWorkspace workspace) {
        fillDenseMatrix(dBasisGradientAccumulator, workspace.yAdjoint);
        blockParameter.fillRAndRinv(workspace.temp1.data, workspace.temp2.data);
        rotateDBasisDiffusionGradientToOriginalBasis(
                workspace.temp2,
                workspace.yAdjoint,
                workspace.temp1,
                workspace.temp3);
        addDenseMatrixToFlatArray(workspace.temp3, gradientAccumulator);
    }

    protected void accumulateNativeTransitionPreparedPullback(
            final BlockDiagonalPreparedBranchBasis prepared,
            final double[] dLogL_dF,
            final double[] dLogL_df,
            final BlockDiagonalBranchGradientWorkspace workspace,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator) {
        fillTotalUpstreamOnTransition(
                prepared.stationaryMean, dLogL_dF, dLogL_df, workspace.upstreamF);
        accumulateSpecializedTransitionPullback(
                prepared.rMatrix,
                prepared.rinvMatrix,
                prepared.transitionMatrix,
                prepared.blockDParams,
                prepared.dt,
                workspace.upstreamF,
                workspace.upstreamFD,
                workspace.frechetHelper,
                workspace.gradD,
                workspace.gradR,
                workspace.temp1,
                workspace.temp2,
                workspace.temp3,
                compressedDAccumulator,
                rotationAccumulator);
    }

    protected void accumulateNativeCovariancePreparedPullback(
            final BlockDiagonalPreparedBranchBasis prepared,
            final MatrixParameterInterface diffusionMatrix,
            final double[] dLogL_dV,
            final BlockDiagonalBranchGradientWorkspace workspace,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator) {
        accumulateSpecializedCovariancePullback(
                prepared.rMatrix,
                prepared.rinvMatrix,
                prepared.expD,
                prepared.blockDParams,
                prepared.dt,
                workspace.qDBasis,
                workspace.stationaryCovDBasis,
                workspace.transitionCovDBasis,
                dLogL_dV,
                workspace,
                compressedDAccumulator,
                rotationAccumulator);
    }

    protected void accumulatePreparedCovarianceAndDiffusionGradient(
            final BlockDiagonalPreparedBranchBasis prepared,
            final MatrixParameterInterface diffusionMatrix,
            final double[] dLogL_dV,
            final BlockDiagonalBranchGradientWorkspace workspace,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator,
            final boolean delayDiffusionGradientRotation,
            final double[] diffusionGradientAccumulator) {
        fillSymmetricDenseMatrixFlat(dLogL_dV, false, workspace.gV);
        fillSharedCovarianceAdjointsInDBasis(
                prepared.rMatrix,
                prepared.rinvMatrix,
                prepared.expD,
                workspace.gV,
                workspace.hDBasis,
                workspace.gS,
                workspace.temp1);
        symmetrize(workspace.gS);

        accumulateSpecializedCovariancePullbackFromShared(
                prepared.rMatrix,
                prepared.rinvMatrix,
                prepared.expD,
                prepared.blockDParams,
                prepared.dt,
                workspace.qDBasis,
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
                compressedDAccumulator,
                rotationAccumulator);

        fillSymmetricSandwichOutputAdjoint(dLogL_dV, false, workspace.gV);
        fillDiffusionGradientDBasisPrepared(prepared, workspace);
        if (delayDiffusionGradientRotation) {
            addDenseMatrixToFlatArray(workspace.yAdjoint, diffusionGradientAccumulator);
        } else {
            rotateDBasisDiffusionGradientToOriginalBasis(
                    prepared.rinvMatrix,
                    workspace.yAdjoint,
                    workspace.temp1,
                    workspace.temp2);
            addDenseMatrixToFlatArray(workspace.temp2, diffusionGradientAccumulator);
        }
    }

    protected void fillNativeTransitionCovarianceMatrix(final MatrixParameterInterface diffusionMatrix,
                                                        final DenseMatrix64F out) {
        fillTransitionCovarianceMatrix(diffusionMatrix, out);
    }

    protected void accumulateNativeTransitionCurrentPullback(final double dt,
                                                             final double[] stationaryMean,
                                                             final double[] dLogL_dF,
                                                             final double[] dLogL_df,
                                                             final double[] compressedDAccumulator,
                                                             final double[] rotationAccumulator) {
        final BlockDiagonalBranchGradientWorkspace workspace = currentGradientWorkspace;
        fillTotalUpstreamOnTransition(
                stationaryMean, dLogL_dF, dLogL_df, workspace.upstreamF);
        accumulateSpecializedTransitionPullback(
                basisCache.rMatrix,
                basisCache.rinvMatrix,
                basisCache.transitionMatrix,
                basisCache.blockDParams,
                dt,
                workspace.upstreamF,
                workspace.upstreamFD,
                workspace.frechetHelper,
                workspace.gradD,
                workspace.gradR,
                workspace.temp1,
                workspace.temp2,
                workspace.temp3,
                compressedDAccumulator,
                rotationAccumulator);
    }

    protected void accumulateNativeCovarianceCurrentCachedPullback(
            final MatrixParameterInterface diffusionMatrix,
            final double dt,
            final double[] dLogL_dV,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator) {
        final BlockDiagonalBranchGradientWorkspace workspace = currentGradientWorkspace;
        accumulateSpecializedCovariancePullback(
                basisCache.rMatrix,
                basisCache.rinvMatrix,
                basisCache.expD,
                basisCache.blockDParams,
                dt,
                qDBasis,
                stationaryCovDBasis,
                transitionCovDBasis,
                dLogL_dV,
                workspace,
                compressedDAccumulator,
                rotationAccumulator);
    }

    protected void accumulateDiffusionGradientCurrentFlat(final double[] dLogL_dV,
                                                          final boolean transposeAdjoint,
                                                          final double[] gradientAccumulator) {
        final BlockDiagonalBranchGradientWorkspace workspace = currentGradientWorkspace;
        fillSymmetricSandwichOutputAdjoint(dLogL_dV, transposeAdjoint, workspace.gV);
        fillDiffusionGradientDBasis(
                basisCache.rMatrix,
                basisCache.rinvMatrix,
                basisCache.expD,
                basisCache.blockDParams,
                workspace.gV,
                workspace.hDBasis,
                workspace.gS,
                workspace.yAdjoint,
                workspace.temp1,
                workspace.temp2,
                workspace.lyapunovAdjointHelper);
        rotateDBasisDiffusionGradientToOriginalBasis(
                basisCache.rinvMatrix,
                workspace.yAdjoint,
                workspace.temp1,
                workspace.temp2);
        addDenseMatrixToFlatArray(workspace.temp2, gradientAccumulator);
    }

    protected void fillTransitionCovarianceMatrix(final MatrixParameterInterface diffusionMatrix,
                                                  final double dt,
                                                  final DenseMatrix64F out) {
        refreshBasisCaches(dt);
        fillTransitionCovarianceMatrix(diffusionMatrix, out);
    }

    protected void fillTransitionCovarianceMatrix(final MatrixParameterInterface diffusionMatrix,
                                                  final DenseMatrix64F out) {
        BlockDiagonalTransitionCovarianceSolver.fillTransitionCovariance(
                diffusionMatrix,
                basisCache.rMatrix,
                basisCache.rinvMatrix,
                basisCache.expD,
                basisCache.blockDParams,
                lyapunovSolver,
                qMatrix,
                qDBasis,
                stationaryCovDBasis,
                transitionCovDBasis,
                tempMatrix,
                out,
                basisCache.blockStarts,
                basisCache.blockSizes);
    }

    protected void fillTransitionCovarianceMatrixPrepared(final BlockDiagonalPreparedBranchBasis prepared,
                                                          final MatrixParameterInterface diffusionMatrix,
                                                          final BlockDiagonalBranchGradientWorkspace workspace,
                                                          final DenseMatrix64F out) {
        BlockDiagonalTransitionCovarianceSolver.fillTransitionCovariance(
                diffusionMatrix,
                prepared.rMatrix,
                prepared.rinvMatrix,
                prepared.expD,
                prepared.blockDParams,
                workspace.lyapunovSolver,
                workspace.qMatrix,
                workspace.qDBasis,
                workspace.stationaryCovDBasis,
                workspace.transitionCovDBasis,
                workspace.temp1,
                out,
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes());
    }

    protected void fillCanonicalTransitionDirectPrepared(final BlockDiagonalPreparedBranchBasis prepared,
                                                         final DenseMatrix64F transitionCovariance,
                                                         final BlockDiagonalBranchGradientWorkspace workspace,
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

    protected synchronized void fillPreparedBranchBasis(final double dt,
                                                        final double[] stationaryMean,
                                                        final BlockDiagonalPreparedBranchBasis prepared) {
        final int dimension = prepared.dimension;
        if (stationaryMean.length != dimension) {
            throw new IllegalArgumentException(
                    "stationaryMean must have length " + dimension + " but has " + stationaryMean.length);
        }

        prepared.dt = dt;
        System.arraycopy(stationaryMean, 0, prepared.stationaryMean, 0, dimension);
        blockParameter.fillBlockDiagonalDecomposition(prepared.decomposition);
        expSolver.computeCompressed(prepared.blockDParams, dt, prepared.expD);
        BlockDiagonalMatrixOps.fillTransitionMatrix(
                prepared.rMatrix.data,
                prepared.expD,
                prepared.rinvMatrix.data,
                dimension,
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes(),
                prepared.workMatrix,
                prepared.transitionMatrix.data);
        prepared.covariancePrepared = false;
    }

    protected void fillTransitionOffset(final double[] stationaryMean,
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

    protected void fillTransitionOffsetPrepared(final BlockDiagonalPreparedBranchBasis prepared,
                                                final double[] out) {
        final int dimension = getDimension();
        final double[] transitionData = prepared.transitionMatrix.data;
        final double[] stationaryMean = prepared.stationaryMean;
        for (int i = 0; i < dimension; ++i) {
            double transformedMean = 0.0;
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                transformedMean += transitionData[rowOffset + j] * stationaryMean[j];
            }
            out[i] = stationaryMean[i] - transformedMean;
        }
    }

    protected void refreshBasisCaches(final double dt) {
        basisCache.refresh(dt);
    }

    protected void refreshMeanGradientCaches(final double dt) {
        basisCache.refreshExpOnly(dt);
    }

    protected double copyAndInvertPositiveDefiniteFlat(final DenseMatrix64F source,
                                                       final double[] matrixOut,
                                                       final double[] inverseOut) {
        return OrthogonalBlockPositiveDefiniteInverter.copyAndInvertFlat(
                source, matrixOut, inverseOut, choleskyScratch, lowerInverseScratch);
    }

    protected void accumulateCurrentNativeGradientFromLocalAdjoints(
            final MatrixParameterInterface diffusionMatrix,
            final double[] stationaryMean,
            final double dt,
            final CanonicalLocalTransitionAdjoints localAdjoints,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator) {
        accumulateNativeTransitionCurrentPullback(
                dt,
                stationaryMean,
                localAdjoints.dLogL_dF,
                localAdjoints.dLogL_df,
                compressedDAccumulator,
                rotationAccumulator);
        checkNativeGradientFinite("transition", dt, compressedDAccumulator, rotationAccumulator);

        accumulateNativeCovarianceCurrentCachedPullback(
                diffusionMatrix,
                dt,
                localAdjoints.dLogL_dOmega,
                compressedDAccumulator,
                rotationAccumulator);
        checkNativeGradientFinite("covariance", dt, compressedDAccumulator, rotationAccumulator);
    }

    protected void accumulateMeanGradientFromTransitionData(final double[] transitionData,
                                                           final double[] dLogL_df,
                                                           final double[] gradientAccumulator,
                                                           final int dimension) {
        if (gradientAccumulator.length == 1) {
            gradientAccumulator[0] += scalarMeanGradientFromTransitionData(
                    transitionData, dLogL_df, dimension);
            return;
        }
        if (gradientAccumulator.length != dimension) {
            throw new IllegalArgumentException(
                    "Stationary-mean gradient length must be 1 or "
                            + dimension + ", found " + gradientAccumulator.length);
        }
        for (int j = 0; j < dimension; ++j) {
            double transformedAdjoint = 0.0;
            for (int i = 0; i < dimension; ++i) {
                transformedAdjoint += transitionData[i * dimension + j] * dLogL_df[i];
            }
            gradientAccumulator[j] += dLogL_df[j] - transformedAdjoint;
        }
    }

    protected double scalarMeanGradientFromTransitionData(final double[] transitionData,
                                                         final double[] dLogL_df,
                                                         final int dimension) {
        double sum = 0.0;
        for (int j = 0; j < dimension; ++j) {
            double transformedAdjoint = 0.0;
            for (int i = 0; i < dimension; ++i) {
                transformedAdjoint += transitionData[i * dimension + j] * dLogL_df[i];
            }
            sum += dLogL_df[j] - transformedAdjoint;
        }
        return sum;
    }

    private void accumulateSpecializedTransitionPullback(
            final DenseMatrix64F rMatrix,
            final DenseMatrix64F rinvMatrix,
            final DenseMatrix64F transitionMatrix,
            final double[] blockDParams,
            final double dt,
            final DenseMatrix64F upstreamF,
            final DenseMatrix64F upstreamFD,
            final BlockDiagonalFrechetHelper frechetHelper,
            final DenseMatrix64F gradD,
            final DenseMatrix64F gradR,
            final DenseMatrix64F temp1,
            final DenseMatrix64F temp2,
            final DenseMatrix64F temp3,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator) {
        CommonOps.multTransA(rMatrix, upstreamF, temp1);
        CommonOps.multTransB(temp1, rinvMatrix, upstreamFD);
        frechetHelper.frechetAdjointExpInDBasis(blockDParams, upstreamFD, dt, gradD);
        transposeInPlace(gradD);
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        CommonOps.multTransB(upstreamF, transitionMatrix, temp1);
        CommonOps.multTransA(transitionMatrix, upstreamF, temp2);
        CommonOps.subtract(temp1, temp2, temp3);
        CommonOps.multTransB(temp3, rinvMatrix, gradR);
        addDenseMatrixToFlatArray(gradR, rotationAccumulator);
    }

    private static void fillTotalUpstreamOnTransition(final double[] stationaryMean,
                                                      final double[] dLogL_dF,
                                                      final double[] dLogL_df,
                                                      final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] outData = out.data;
        for (int row = 0; row < dimension; ++row) {
            final int rowOffset = row * dimension;
            for (int col = 0; col < dimension; ++col) {
                outData[rowOffset + col] =
                        dLogL_dF[rowOffset + col] - dLogL_df[row] * stationaryMean[col];
            }
        }
    }

    private void accumulateSpecializedCovariancePullback(
            final DenseMatrix64F rMatrix,
            final DenseMatrix64F rinvMatrix,
            final double[] expD,
            final double[] blockDParams,
            final double dt,
            final DenseMatrix64F qDBasis,
            final DenseMatrix64F stationaryCovDBasis,
            final DenseMatrix64F transitionCovDBasis,
            final double[] dLogL_dV,
            final BlockDiagonalBranchGradientWorkspace workspace,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator) {
        fillSymmetricDenseMatrixFlat(dLogL_dV, false, workspace.gV);
        fillSharedCovarianceAdjointsInDBasis(
                rMatrix,
                rinvMatrix,
                expD,
                workspace.gV,
                workspace.hDBasis,
                workspace.gS,
                workspace.temp1);
        symmetrize(workspace.gS);

        accumulateSpecializedCovariancePullbackFromShared(
                rMatrix,
                rinvMatrix,
                expD,
                blockDParams,
                dt,
                qDBasis,
                stationaryCovDBasis,
                transitionCovDBasis,
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
                compressedDAccumulator,
                rotationAccumulator);
    }

    private void accumulateSpecializedCovariancePullbackFromShared(
            final DenseMatrix64F rMatrix,
            final DenseMatrix64F rinvMatrix,
            final double[] expD,
            final double[] blockDParams,
            final double dt,
            final DenseMatrix64F qDBasis,
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
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator) {
        CommonOps.fill(gradD, 0.0);
        lyapunovAdjointHelper.accumulateLyapunovContributionInDBasis(
                stationaryCovDBasis, gS, blockDParams, gradD);
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        lyapunovAdjointHelper.solveAdjointInDBasis(gS, blockDParams, yAdjoint);
        CommonOps.mult(yAdjoint, qDBasis, temp1);
        CommonOps.multTransA(yAdjoint, qDBasis, temp2);
        CommonOps.add(temp1, temp2, temp3);
        CommonOps.multTransA(rinvMatrix, temp3, gradR);

        BlockDiagonalMatrixOps.multiplyRightBlockDiagonal(
                hDBasis.data,
                expD,
                hDBasis.numRows,
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes(),
                temp1.data);
        CommonOps.mult(temp1, stationaryCovDBasis, gECov);
        CommonOps.scale(-2.0, gECov);
        frechetHelper.frechetAdjointExpInDBasis(blockDParams, gECov, dt, gradD);
        transposeInPlace(gradD);
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        CommonOps.mult(gV, rMatrix, temp1);
        CommonOps.mult(temp1, transitionCovDBasis, temp2);
        CommonOps.scale(2.0, temp2);
        CommonOps.addEquals(gradR, temp2);
        addDenseMatrixToFlatArray(gradR, rotationAccumulator);
    }

    private void fillSharedCovarianceAdjointsInDBasis(final DenseMatrix64F rMatrix,
                                                      final DenseMatrix64F rinvMatrix,
                                                      final double[] expD,
                                                      final DenseMatrix64F gV,
                                                      final DenseMatrix64F hDBasis,
                                                      final DenseMatrix64F gS,
                                                      final DenseMatrix64F temp) {
        CommonOps.multTransA(rMatrix, gV, temp);
        CommonOps.mult(temp, rMatrix, hDBasis);
        symmetrize(hDBasis);

        BlockDiagonalMatrixOps.multiplyBlockDiagonalLeftTranspose(
                expD,
                hDBasis.data,
                temp.data,
                temp.numRows,
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes());
        BlockDiagonalMatrixOps.multiplyRightBlockDiagonal(
                temp.data,
                expD,
                temp.numRows,
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes(),
                gS.data);
        OrthogonalBlockDenseMatrixOps.subtract(hDBasis, gS, gS);
    }

    private void fillDiffusionGradientDBasisPrepared(final BlockDiagonalPreparedBranchBasis prepared,
                                                     final BlockDiagonalBranchGradientWorkspace workspace) {
        fillDiffusionGradientDBasis(
                prepared.rMatrix,
                prepared.rinvMatrix,
                prepared.expD,
                prepared.blockDParams,
                workspace.gV,
                workspace.hDBasis,
                workspace.gS,
                workspace.yAdjoint,
                workspace.temp1,
                workspace.temp2,
                workspace.lyapunovAdjointHelper);
    }

    private void fillDiffusionGradientDBasis(final DenseMatrix64F rMatrix,
                                             final DenseMatrix64F rinvMatrix,
                                             final double[] expD,
                                             final double[] blockDParams,
                                             final DenseMatrix64F gV,
                                             final DenseMatrix64F hDBasis,
                                             final DenseMatrix64F gS,
                                             final DenseMatrix64F yAdjoint,
                                             final DenseMatrix64F temp1,
                                             final DenseMatrix64F temp2,
                                             final dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper) {
        fillSharedCovarianceAdjointsInDBasis(
                rMatrix,
                rinvMatrix,
                expD,
                gV,
                hDBasis,
                gS,
                temp1);
        symmetrize(gS);
        yAdjoint.set(gS);
        CommonOps.scale(-1.0, yAdjoint);
        lyapunovAdjointHelper.solveAdjointInDBasis(yAdjoint, blockDParams, yAdjoint);
    }

    private void rotateDBasisDiffusionGradientToOriginalBasis(final DenseMatrix64F rinvMatrix,
                                                              final DenseMatrix64F dBasisGradient,
                                                              final DenseMatrix64F temp,
                                                              final DenseMatrix64F out) {
        MatrixOps.fillSymmetricSandwichTransposeRightOutputAdjoint(
                dBasisGradient.data, false, out.data, dBasisGradient.numRows);
        MatrixOps.symmetricSandwichTransposeRightMiddleAdjointFromOutputAdjoint(
                rinvMatrix.data, out.data, out.data, temp.data, rinvMatrix.numRows);
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
            compressedAccumulator[upperBase + blockIndex] += data[start * dimension + start + 1];
            compressedAccumulator[lowerBase + blockIndex] += data[(start + 1) * dimension + start];
            blockIndex++;
        }
    }

    private static void fillDenseMatrix(final double[] source, final DenseMatrix64F out) {
        System.arraycopy(source, 0, out.data, 0, out.numRows * out.numCols);
    }

    private static void fillSymmetricDenseMatrixFlat(final double[] source,
                                                     final boolean transposeSource,
                                                     final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                final double ij = transposeSource ? source[j * dimension + i] : source[rowOffset + j];
                final double ji = transposeSource ? source[rowOffset + j] : source[j * dimension + i];
                data[rowOffset + j] = 0.5 * (ij + ji);
            }
        }
    }

    private static void fillSymmetricSandwichOutputAdjoint(final double[] source,
                                                           final boolean transposeSource,
                                                           final DenseMatrix64F out) {
        MatrixOps.fillSymmetricSandwichTransposeRightOutputAdjoint(
                source, transposeSource, out.data, out.numRows);
    }

    private static void addDenseMatrixToFlatArray(final DenseMatrix64F src, final double[] dest) {
        final int length = src.numRows * src.numCols;
        final double[] data = src.data;
        for (int i = 0; i < length; ++i) {
            dest[i] += data[i];
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

    private void flattenRotationGradientForParameter(final double[] rowMajorGradient,
                                                     final double[] out,
                                                     final int offset,
                                                     final int dimension) {
        if (blockParameter.getRotationMatrixParameter() instanceof TransposedMatrixParameter) {
            System.arraycopy(rowMajorGradient, 0, out, offset, dimension * dimension);
            return;
        }
        int index = offset;
        for (int col = 0; col < dimension; ++col) {
            for (int row = 0; row < dimension; ++row) {
                out[index++] = rowMajorGradient[row * dimension + col];
            }
        }
    }

    private void accumulateRotationGradientForParameter(final double[] rowMajorGradient,
                                                        final double[] out,
                                                        final int offset,
                                                        final int dimension) {
        if (blockParameter.getRotationMatrixParameter() instanceof TransposedMatrixParameter) {
            for (int i = 0; i < dimension * dimension; ++i) {
                out[offset + i] += rowMajorGradient[i];
            }
            return;
        }
        int index = offset;
        for (int col = 0; col < dimension; ++col) {
            for (int row = 0; row < dimension; ++row) {
                out[index++] += rowMajorGradient[row * dimension + col];
            }
        }
    }

    protected void checkNativeGradientFinite(final String contribution,
                                             final double dt,
                                             final double[] compressedDAccumulator,
                                             final double[] rotationAccumulator) {
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException(
                    "Non-finite " + nativeGradientPathName()
                            + " native " + contribution + " contribution at dt=" + dt);
        }
    }

    protected String nativeGradientPathName() {
        return "block-diagonal";
    }

    protected UnsupportedOperationException unsupportedNativeGradient(final String operation) {
        return new UnsupportedOperationException(
                operation + " is not supported by this block-diagonal backend.");
    }

    protected BlockDiagonalPreparedBranchBasis asPreparedBasis(
            final CanonicalPreparedBranchHandle prepared) {
        if (!(prepared instanceof BlockDiagonalPreparedBranchHandle)) {
            throw new IllegalArgumentException("Expected a BlockDiagonalPreparedBranchHandle.");
        }
        return ((BlockDiagonalPreparedBranchHandle) prepared).getBasis();
    }

    protected BlockDiagonalBranchGradientWorkspace asBranchWorkspace(
            final CanonicalBranchWorkspace workspace) {
        if (!(workspace instanceof BlockDiagonalBranchGradientWorkspace)) {
            throw new IllegalArgumentException("Expected a BlockDiagonalBranchGradientWorkspace.");
        }
        return (BlockDiagonalBranchGradientWorkspace) workspace;
    }

    protected static boolean isFinite(final double[] values) {
        for (double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    private void copyDenseMatrixToArray(final DenseMatrix64F matrix, final double[][] out) {
        final int dimension = getDimension();
        if (out == null || out.length != dimension) {
            throw new IllegalArgumentException("transition matrix must have " + dimension + " rows");
        }
        final double[] data = matrix.data;
        for (int i = 0; i < dimension; ++i) {
            if (out[i] == null || out[i].length != dimension) {
                throw new IllegalArgumentException("transition matrix row " + i
                        + " must have length " + dimension);
            }
            System.arraycopy(data, i * dimension, out[i], 0, dimension);
        }
    }
}
