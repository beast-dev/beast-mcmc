package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalFrechetHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovAdjointHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovSolver;
import dr.evomodel.continuous.ou.DenseSelectionMatrixParameterization;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import org.ejml.data.DenseMatrix64F;

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
    private final OrthogonalBlockSelectionAdjoint selectionAdjoint;
    private final OrthogonalBlockCovarianceAdjoint covarianceAdjoint;
    private final DenseMatrix64F transitionCovariance;
    private final double[] transitionMatrixArrayScratch;
    private final double[] transitionCovarianceArrayScratch;
    private final double[] precisionFlat;
    private final double[] nativeBlockGradientScratch;
    private final double[] transitionOffsetScratch;
    private final double[] denseAdjointScratch;
    private final double[] choleskyScratch;
    private final double[] lowerInverseScratch;
    private final double[] legacyMatrixScratch;
    private final double[] legacyRotationScratch;
    private final CanonicalTransitionAdjointUtils.Workspace canonicalAdjointWorkspace;

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
        this.selectionAdjoint = new OrthogonalBlockSelectionAdjoint(
                blockParameter, orthogonalRotation, frechetHelper, denseFallbackPolicy);
        this.covarianceAdjoint = new OrthogonalBlockCovarianceAdjoint(
                blockParameter, lyapunovSolver, lyapunovAdjointHelper,
                frechetHelper, denseFallbackPolicy);

        final int d = getDimension();
        this.transitionCovariance = new DenseMatrix64F(d, d);
        this.transitionMatrixArrayScratch = new double[d * d];
        this.transitionCovarianceArrayScratch = new double[d * d];
        this.precisionFlat = new double[d * d];
        this.nativeBlockGradientScratch = new double[blockParameter.getBlockDiagonalNParameters()];
        this.transitionOffsetScratch = new double[d];
        this.denseAdjointScratch = new double[d * d];
        this.choleskyScratch = new double[d * d];
        this.lowerInverseScratch = new double[d * d];
        this.legacyMatrixScratch = new double[d * d];
        this.legacyRotationScratch = new double[d * d];
        this.canonicalAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(d);
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
    public void fillTransitionCovarianceFlat(final MatrixParameterInterface diffusionMatrix,
                                             final double dt,
                                             final double[] out) {
        refreshBasisCaches(dt);
        transitionFactory.fillTransitionCovarianceFlat(diffusionMatrix, basisCache, out);
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
        return new OrthogonalBlockPreparedBranchBasis(
                dimension,
                blockParameter.getTridiagonalDDimension(),
                blockParameter.getCompressedDDimension());
    }

    @Override
    public CanonicalPreparedBranchHandle createPreparedBranchHandle() {
        return new OrthogonalBlockPreparedBranchHandle(createPreparedBranchBasis());
    }

    @Override
    public int getSelectionGradientDimension() {
        return blockParameter.getBlockDiagonalNParameters()
                + orthogonalRotation.getOrthogonalParameter().getDimension();
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
        blockParameter.chainGradient(compressedGradient, nativeGradientScratch);
        System.arraycopy(nativeGradientScratch, 0, gradientOut, 0, nativeBlockDim);
        orthogonalRotation.fillPullBackGradientFlat(
                rotationGradientFlat,
                getDimension(),
                gradientOut,
                nativeBlockDim);
    }

    @Override
    public CanonicalBranchWorkspace createBranchWorkspace() {
        return createBranchGradientWorkspace();
    }

    @Override
    public void prepareBranch(final double dt,
                              final double[] stationaryMean,
                              final CanonicalPreparedBranchHandle prepared) {
        prepareBranchBasis(dt, stationaryMean, asPreparedBasis(prepared));
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
        final OrthogonalBlockPreparedBranchBasis basis = asPreparedBasis(prepared);
        final OrthogonalBlockBranchGradientWorkspace branchWorkspace = asBranchWorkspace(workspace);
        loadOrFillPreparedCovariance(basis, diffusionMatrix, branchWorkspace);
        MatrixOps.toFlat(basis.transitionMatrix, transitionMatrixOut, getDimension());
        fillTransitionOffsetPrepared(basis, transitionOffsetOut);
        MatrixOps.toFlat(basis.transitionCovariance, transitionCovarianceOut, getDimension());
        return true;
    }

    public void prepareBranchCovariance(final OrthogonalBlockPreparedBranchBasis prepared,
                                        final MatrixParameterInterface diffusionMatrix,
                                        final OrthogonalBlockBranchGradientWorkspace workspace) {
        fillTransitionCovarianceMatrixPrepared(prepared, diffusionMatrix, workspace, workspace.transitionCovariance);
        OrthogonalBlockTransitionCovarianceSolver.storePreparedCovariance(prepared, workspace);
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

    public void accumulateNativeSelectionAndDiffusionGradientFromAdjointsPreparedFlat(
            final OrthogonalBlockPreparedBranchBasis prepared,
            final MatrixParameterInterface diffusionMatrix,
            final CanonicalLocalTransitionAdjoints localAdjoints,
            final OrthogonalBlockBranchGradientWorkspace workspace,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator,
            final boolean delayDiffusionGradientRotation,
            final double[] diffusionGradientAccumulator) {
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

        covarianceAdjoint.accumulatePreparedCovarianceAndDiffusionGradientFlat(
                prepared,
                localAdjoints.dLogL_dOmega,
                workspace,
                compressedDAccumulator,
                rotationAccumulator,
                delayDiffusionGradientRotation,
                diffusionGradientAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException(
                    "Non-finite orthogonal native covariance contribution at dt=" + prepared.dt);
        }
    }

    @Override
    public void accumulateNativeGradientFromAdjointsPreparedFlat(final CanonicalPreparedBranchHandle prepared,
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

    @Override
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
        selectionAdjoint.accumulateMeanGradientPrepared(
                prepared, dLogL_df, gradientAccumulator, workspace);
    }

    public void accumulateMeanGradientPreparedDBasisFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                                         final double[] dLogL_df,
                                                         final double[] dBasisGradientAccumulator,
                                                         final OrthogonalBlockBranchGradientWorkspace workspace) {
        selectionAdjoint.accumulateMeanGradientPreparedDBasis(
                prepared, dLogL_df, dBasisGradientAccumulator, workspace);
    }

    public void finishMeanGradientFromDBasisFlat(final double[] dBasisGradientAccumulator,
                                                 final double[] gradientAccumulator,
                                                 final OrthogonalBlockBranchGradientWorkspace workspace) {
        selectionAdjoint.finishMeanGradientFromDBasis(
                dBasisGradientAccumulator, gradientAccumulator, workspace);
    }

    @Override
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

    @Override
    public boolean supportsDelayedMeanGradientRotation() {
        return true;
    }

    @Override
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

    @Override
    public void finishMeanGradientFromDBasisFlat(final double[] dBasisGradientAccumulator,
                                                 final double[] gradientAccumulator,
                                                 final CanonicalBranchWorkspace workspace) {
        finishMeanGradientFromDBasisFlat(
                dBasisGradientAccumulator,
                gradientAccumulator,
                asBranchWorkspace(workspace));
    }

    public void accumulateDiffusionGradientPreparedFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                                        final double[] dLogL_dV,
                                                        final boolean transposeAdjoint,
                                                        final double[] gradientAccumulator,
                                                        final OrthogonalBlockBranchGradientWorkspace workspace) {
        covarianceAdjoint.accumulateDiffusionGradientPreparedFlat(
                prepared, dLogL_dV, gradientAccumulator, workspace);
    }

    public void accumulateDiffusionGradientPreparedDBasisFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                                              final double[] dLogL_dV,
                                                              final boolean transposeAdjoint,
                                                              final double[] dBasisGradientAccumulator,
                                                              final OrthogonalBlockBranchGradientWorkspace workspace) {
        covarianceAdjoint.accumulateDiffusionGradientPreparedDBasisFlat(
                prepared, dLogL_dV, dBasisGradientAccumulator, workspace);
    }

    public void finishDiffusionGradientFromDBasisFlat(final double[] dBasisGradientAccumulator,
                                                      final double[] gradientAccumulator,
                                                      final OrthogonalBlockBranchGradientWorkspace workspace) {
        orthogonalRotation.fillOrthogonalMatrix(workspace.temp1.data);
        orthogonalRotation.fillOrthogonalTranspose(workspace.temp2.data);
        covarianceAdjoint.finishDiffusionGradientFromDBasisFlat(
                workspace.temp1, workspace.temp2, dBasisGradientAccumulator, gradientAccumulator, workspace);
    }

    @Override
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

    @Override
    public boolean supportsDelayedDiffusionGradientRotation() {
        return true;
    }

    @Override
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

    @Override
    public void finishDiffusionGradientFromDBasisFlat(final double[] dBasisGradientAccumulator,
                                                      final double[] gradientAccumulator,
                                                      final CanonicalBranchWorkspace workspace) {
        finishDiffusionGradientFromDBasisFlat(
                dBasisGradientAccumulator,
                gradientAccumulator,
                asBranchWorkspace(workspace));
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

    private static OrthogonalBlockPreparedBranchBasis asPreparedBasis(
            final CanonicalPreparedBranchHandle prepared) {
        if (!(prepared instanceof OrthogonalBlockPreparedBranchHandle)) {
            throw new IllegalArgumentException("Expected an OrthogonalBlockPreparedBranchHandle.");
        }
        return ((OrthogonalBlockPreparedBranchHandle) prepared).getBasis();
    }

    private static OrthogonalBlockBranchGradientWorkspace asBranchWorkspace(
            final CanonicalBranchWorkspace workspace) {
        if (!(workspace instanceof OrthogonalBlockBranchGradientWorkspace)) {
            throw new IllegalArgumentException("Expected an OrthogonalBlockBranchGradientWorkspace.");
        }
        return (OrthogonalBlockBranchGradientWorkspace) workspace;
    }

    private void fillTransitionCovarianceMatrix(final MatrixParameterInterface diffusionMatrix,
                                                final double dt,
                                                final DenseMatrix64F out) {
        refreshBasisCaches(dt);
        covarianceAdjoint.fillTransitionCovarianceMatrix(diffusionMatrix, basisCache, out);
    }

    public void accumulateNativeGradientFromTransitionFlat(final double dt,
                                                           final double[] stationaryMean,
                                                           final double[] dLogL_dF,
                                                           final double[] dLogL_df,
                                                           final double[] compressedDAccumulator,
                                                           final double[] rotationAccumulator) {
        refreshBasisCaches(dt);
        selectionAdjoint.accumulateCurrentFlat(
                basisCache, dt, stationaryMean, dLogL_dF, dLogL_df,
                compressedDAccumulator, rotationAccumulator);
    }

    public void accumulateNativeGradientFromTransition(final double dt,
                                                       final double[] stationaryMean,
                                                       final double[][] dLogL_dF,
                                                       final double[] dLogL_df,
                                                      final double[] compressedDAccumulator,
                                                      final double[][] rotationAccumulator) {
        final int d = getDimension();
        Arrays.fill(legacyRotationScratch, 0.0);
        copyMatrixToFlat(dLogL_dF, legacyMatrixScratch, d);
        accumulateNativeGradientFromTransitionFlat(
                dt, stationaryMean, legacyMatrixScratch, dLogL_df,
                compressedDAccumulator, legacyRotationScratch);
        addFlatToMatrix(legacyRotationScratch, rotationAccumulator, d);
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
        Arrays.fill(denseAdjointScratch, 0.0);

        accumulateNativeGradientFromTransitionFlat(
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
                precisionFlat,
                gradientAccumulator);
    }

    public void accumulateNativeGradientFromCovarianceStationaryFlat(final MatrixParameterInterface diffusionMatrix,
                                                                     final double dt,
                                                                     final double[] dLogL_dV,
                                                                     final double[] compressedDAccumulator,
                                                                     final double[] rotationAccumulator) {
        refreshBasisCaches(dt);
        covarianceAdjoint.prepareAndAccumulateCurrentFlat(
                diffusionMatrix, basisCache, dt, dLogL_dV,
                compressedDAccumulator, rotationAccumulator);
    }

    public void accumulateNativeGradientFromCovarianceStationary(final MatrixParameterInterface diffusionMatrix,
                                                                 final double dt,
                                                                 final double[][] dLogL_dV,
                                                                 final double[] compressedDAccumulator,
                                                                 final double[][] rotationAccumulator) {
        final int d = getDimension();
        Arrays.fill(legacyRotationScratch, 0.0);
        copyMatrixToFlat(dLogL_dV, legacyMatrixScratch, d);
        accumulateNativeGradientFromCovarianceStationaryFlat(
                diffusionMatrix, dt, legacyMatrixScratch,
                compressedDAccumulator, legacyRotationScratch);
        addFlatToMatrix(legacyRotationScratch, rotationAccumulator, d);
    }

    @Override
    public void accumulateDiffusionGradientFlat(final MatrixParameterInterface diffusionMatrix,
                                                final double dt,
                                                final double[] dLogL_dV,
                                                final boolean transposeAdjoint,
                                                final double[] gradientAccumulator) {
        refreshBasisCaches(dt);
        covarianceAdjoint.accumulateDiffusionGradientCurrentFlat(
                basisCache, dLogL_dV, transposeAdjoint, gradientAccumulator);
    }

    public void accumulateDiffusionGradient(final MatrixParameterInterface diffusionMatrix,
                                            final double dt,
                                            final double[][] dLogL_dV,
                                            final double[] gradientAccumulator) {
        final int d = getDimension();
        copyMatrixToFlat(dLogL_dV, legacyMatrixScratch, d);
        accumulateDiffusionGradientFlat(diffusionMatrix, dt, legacyMatrixScratch, false, gradientAccumulator);
    }

    public void accumulateNativeGradientFromCanonicalContributionFlat(final MatrixParameterInterface diffusionMatrix,
                                                                      final double[] stationaryMean,
                                                                      final double dt,
                                                                      final CanonicalBranchMessageContribution contribution,
                                                                      final CanonicalLocalTransitionAdjoints localAdjoints,
                                                                      final double[] compressedDAccumulator,
                                                                      final double[] rotationAccumulator) {
        refreshBasisCaches(dt);
        fillTransitionCovarianceMatrix(diffusionMatrix, dt, transitionCovariance);
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

        selectionAdjoint.accumulateCurrentFlat(
                basisCache,
                dt,
                stationaryMean,
                localAdjoints.dLogL_dF,
                localAdjoints.dLogL_df,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException("Non-finite orthogonal native transition contribution at dt=" + dt);
        }
        covarianceAdjoint.accumulateCurrentCachedFlat(
                basisCache, dt, localAdjoints.dLogL_dOmega,
                compressedDAccumulator, rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException("Non-finite orthogonal native covariance contribution at dt=" + dt);
        }
    }

    public void accumulateNativeGradientFromAdjointsFlat(final MatrixParameterInterface diffusionMatrix,
                                                         final double[] stationaryMean,
                                                         final double dt,
                                                         final CanonicalLocalTransitionAdjoints localAdjoints,
                                                         final double[] compressedDAccumulator,
                                                         final double[] rotationAccumulator) {
        refreshBasisCaches(dt);
        fillTransitionCovarianceMatrix(diffusionMatrix, dt, transitionCovariance);

        selectionAdjoint.accumulateCurrentFlat(
                basisCache,
                dt,
                stationaryMean,
                localAdjoints.dLogL_dF,
                localAdjoints.dLogL_df,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException("Non-finite orthogonal native transition contribution at dt=" + dt);
        }

        covarianceAdjoint.accumulateCurrentCachedFlat(
                basisCache, dt, localAdjoints.dLogL_dOmega,
                compressedDAccumulator, rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException("Non-finite orthogonal native covariance contribution at dt=" + dt);
        }
    }

    public void accumulateMeanGradient(final double dt,
                                       final double[] dLogL_df,
                                       final double[] gradientAccumulator) {
        refreshMeanGradientCaches(dt);
        selectionAdjoint.accumulateMeanGradient(basisCache, dLogL_df, gradientAccumulator);
    }

    public double accumulateScalarMeanGradient(final double dt,
                                               final double[] dLogL_df) {
        refreshMeanGradientCaches(dt);
        return selectionAdjoint.accumulateScalarMeanGradient(basisCache, dLogL_df);
    }

    private void fillTransitionCovarianceMatrixPrepared(final OrthogonalBlockPreparedBranchBasis prepared,
                                                        final MatrixParameterInterface diffusionMatrix,
                                                        final OrthogonalBlockBranchGradientWorkspace workspace,
                                                        final DenseMatrix64F out) {
        covarianceAdjoint.fillTransitionCovarianceMatrixPrepared(
                prepared, diffusionMatrix, workspace, out);
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

    private void accumulateNativeGradientFromTransitionPreparedFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                                                    final double[] dLogL_dF,
                                                                    final double[] dLogL_df,
                                                                    final OrthogonalBlockBranchGradientWorkspace workspace,
                                                                    final double[] compressedDAccumulator,
                                                                    final double[] rotationAccumulator) {
        selectionAdjoint.accumulatePreparedFlat(
                prepared, dLogL_dF, dLogL_df, workspace,
                compressedDAccumulator, rotationAccumulator);
    }

    private void accumulateNativeGradientFromCovarianceStationaryPreparedFlat(final OrthogonalBlockPreparedBranchBasis prepared,
                                                                              final double[] dLogL_dV,
                                                                              final OrthogonalBlockBranchGradientWorkspace workspace,
                                                                              final double[] compressedDAccumulator,
                                                                              final double[] rotationAccumulator) {
        covarianceAdjoint.accumulatePreparedFlat(
                prepared, dLogL_dV, workspace,
                compressedDAccumulator, rotationAccumulator);
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

    private void fillTransitionOffsetPrepared(final OrthogonalBlockPreparedBranchBasis prepared,
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

    private void refreshBasisCaches(final double dt) {
        basisCache.refresh(dt);
    }

    private void refreshMeanGradientCaches(final double dt) {
        basisCache.refreshExpOnly(dt);
    }

    private double copyAndInvertPositiveDefiniteFlat(final DenseMatrix64F source,
                                                     final double[] matrixOut,
                                                     final double[] inverseOut) {
        return OrthogonalBlockPositiveDefiniteInverter.copyAndInvertFlat(
                source, matrixOut, inverseOut, choleskyScratch, lowerInverseScratch);
    }

    private static void copyMatrixToFlat(final double[][] source,
                                         final double[] out,
                                         final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(source[i], 0, out, i * dimension, dimension);
        }
    }

    private static void addFlatToMatrix(final double[] source,
                                        final double[][] out,
                                        final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                out[i][j] += source[rowOffset + j];
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
