package dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalFrechetHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovAdjointHelper;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.OrthogonalMatrixProvider;
import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;

/**
 * Orthogonal specialization of the block-diagonal parametrization. The forward
 * path and native-gradient orchestration are inherited; this class supplies the
 * pullback hooks that require {@code R^{-1} = R^T}.
 */
public final class OrthogonalBlockDiagonalSelectionMatrixParameterization
        extends BlockDiagonalSelectionMatrixParameterization
        implements OrthogonalBlockCanonicalParameterization {

    private final OrthogonalMatrixProvider orthogonalRotation;
    private final OrthogonalBlockSelectionAdjoint selectionAdjoint;
    private final OrthogonalBlockCovarianceAdjoint covarianceAdjoint;
    private final double[] nativeBlockGradientScratch;

    public OrthogonalBlockDiagonalSelectionMatrixParameterization(
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final OrthogonalMatrixProvider orthogonalRotation) {
        super(blockParameter);
        this.orthogonalRotation = orthogonalRotation;

        final BlockDiagonalFrechetHelper frechetHelper =
                new BlockDiagonalFrechetHelper(expSolver.getStructure());
        final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper =
                new BlockDiagonalLyapunovAdjointHelper(getDimension(), lyapunovSolver);
        final BlockDiagonalDenseFallbackPolicy denseFallbackPolicy =
                new BlockDiagonalDenseFallbackPolicy();

        this.selectionAdjoint = new OrthogonalBlockSelectionAdjoint(
                blockParameter, orthogonalRotation, frechetHelper, denseFallbackPolicy);
        this.covarianceAdjoint = new OrthogonalBlockCovarianceAdjoint(
                blockParameter, lyapunovSolver, lyapunovAdjointHelper,
                frechetHelper, denseFallbackPolicy);
        this.nativeBlockGradientScratch = new double[blockParameter.getBlockDiagonalNParameters()];
    }

    @Override
    public OrthogonalBlockPreparedBranchBasis createPreparedBranchBasis() {
        final int dimension = getDimension();
        return new OrthogonalBlockPreparedBranchBasis(
                dimension,
                blockParameter.getCompressedDDimension(),
                blockParameter.createBlockDiagonalDecomposition());
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
    public boolean supportsDelayedMeanGradientRotation() {
        return true;
    }

    @Override
    public boolean supportsDelayedDiffusionGradientRotation() {
        return true;
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

    @Override
    protected void accumulateNativeTransitionPreparedPullback(
            final BlockDiagonalPreparedBranchBasis prepared,
            final double[] dLogL_dF,
            final double[] dLogL_df,
            final BlockDiagonalBranchGradientWorkspace workspace,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator) {
        selectionAdjoint.accumulatePreparedFlat(
                asOrthogonalPreparedBasis(prepared),
                dLogL_dF,
                dLogL_df,
                workspace,
                compressedDAccumulator,
                rotationAccumulator);
    }

    @Override
    protected void accumulateNativeCovariancePreparedPullback(
            final BlockDiagonalPreparedBranchBasis prepared,
            final MatrixParameterInterface diffusionMatrix,
            final double[] dLogL_dV,
            final BlockDiagonalBranchGradientWorkspace workspace,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator) {
        covarianceAdjoint.accumulatePreparedFlat(
                asOrthogonalPreparedBasis(prepared),
                dLogL_dV,
                workspace,
                compressedDAccumulator,
                rotationAccumulator);
    }

    @Override
    protected void accumulatePreparedCovarianceAndDiffusionGradient(
            final BlockDiagonalPreparedBranchBasis prepared,
            final MatrixParameterInterface diffusionMatrix,
            final double[] dLogL_dV,
            final BlockDiagonalBranchGradientWorkspace workspace,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator,
            final boolean delayDiffusionGradientRotation,
            final double[] diffusionGradientAccumulator) {
        covarianceAdjoint.accumulatePreparedCovarianceAndDiffusionGradientFlat(
                asOrthogonalPreparedBasis(prepared),
                dLogL_dV,
                workspace,
                compressedDAccumulator,
                rotationAccumulator,
                delayDiffusionGradientRotation,
                diffusionGradientAccumulator);
    }

    @Override
    protected void accumulateMeanGradientPreparedDBasisFlat(final BlockDiagonalPreparedBranchBasis prepared,
                                                            final double[] dLogL_df,
                                                            final double[] dBasisGradientAccumulator,
                                                            final BlockDiagonalBranchGradientWorkspace workspace) {
        selectionAdjoint.accumulateMeanGradientPreparedDBasis(
                asOrthogonalPreparedBasis(prepared),
                dLogL_df,
                dBasisGradientAccumulator,
                workspace);
    }

    @Override
    protected void finishMeanGradientFromDBasisFlat(final double[] dBasisGradientAccumulator,
                                                    final double[] gradientAccumulator,
                                                    final BlockDiagonalBranchGradientWorkspace workspace) {
        selectionAdjoint.finishMeanGradientFromDBasis(
                dBasisGradientAccumulator,
                gradientAccumulator,
                workspace);
    }

    @Override
    protected void accumulateDiffusionGradientPreparedFlat(final BlockDiagonalPreparedBranchBasis prepared,
                                                           final double[] dLogL_dV,
                                                           final boolean transposeAdjoint,
                                                           final double[] gradientAccumulator,
                                                           final BlockDiagonalBranchGradientWorkspace workspace) {
        covarianceAdjoint.accumulateDiffusionGradientPreparedFlat(
                asOrthogonalPreparedBasis(prepared),
                dLogL_dV,
                gradientAccumulator,
                workspace);
    }

    @Override
    protected void accumulateDiffusionGradientPreparedDBasisFlat(final BlockDiagonalPreparedBranchBasis prepared,
                                                                 final double[] dLogL_dV,
                                                                 final boolean transposeAdjoint,
                                                                 final double[] dBasisGradientAccumulator,
                                                                 final BlockDiagonalBranchGradientWorkspace workspace) {
        covarianceAdjoint.accumulateDiffusionGradientPreparedDBasisFlat(
                asOrthogonalPreparedBasis(prepared),
                dLogL_dV,
                dBasisGradientAccumulator,
                workspace);
    }

    @Override
    protected void finishDiffusionGradientFromDBasisFlat(final double[] dBasisGradientAccumulator,
                                                         final double[] gradientAccumulator,
                                                         final BlockDiagonalBranchGradientWorkspace workspace) {
        orthogonalRotation.fillOrthogonalMatrix(workspace.temp1.data);
        orthogonalRotation.fillOrthogonalTranspose(workspace.temp2.data);
        covarianceAdjoint.finishDiffusionGradientFromDBasisFlat(
                workspace.temp1,
                workspace.temp2,
                dBasisGradientAccumulator,
                gradientAccumulator,
                workspace);
    }

    @Override
    protected void fillNativeTransitionCovarianceMatrix(final MatrixParameterInterface diffusionMatrix,
                                                        final DenseMatrix64F out) {
        covarianceAdjoint.fillTransitionCovarianceMatrix(diffusionMatrix, basisCache, out);
    }

    @Override
    protected void accumulateNativeTransitionCurrentPullback(final double dt,
                                                             final double[] stationaryMean,
                                                             final double[] dLogL_dF,
                                                             final double[] dLogL_df,
                                                             final double[] compressedDAccumulator,
                                                             final double[] rotationAccumulator) {
        selectionAdjoint.accumulateCurrentFlat(
                basisCache, dt, stationaryMean, dLogL_dF, dLogL_df,
                compressedDAccumulator, rotationAccumulator);
    }

    @Override
    protected void accumulateNativeCovarianceCurrentCachedPullback(
            final MatrixParameterInterface diffusionMatrix,
            final double dt,
            final double[] dLogL_dV,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator) {
        covarianceAdjoint.accumulateCurrentCachedFlat(
                basisCache, dt, dLogL_dV, compressedDAccumulator, rotationAccumulator);
    }

    @Override
    protected void accumulateDiffusionGradientCurrentFlat(final double[] dLogL_dV,
                                                          final boolean transposeAdjoint,
                                                          final double[] gradientAccumulator) {
        covarianceAdjoint.accumulateDiffusionGradientCurrentFlat(
                basisCache, dLogL_dV, transposeAdjoint, gradientAccumulator);
    }

    @Override
    protected String nativeGradientPathName() {
        return "orthogonal block-diagonal";
    }

    private static OrthogonalBlockPreparedBranchBasis asOrthogonalPreparedBasis(
            final BlockDiagonalPreparedBranchBasis prepared) {
        if (!(prepared instanceof OrthogonalBlockPreparedBranchBasis)) {
            throw new IllegalArgumentException("Expected an OrthogonalBlockPreparedBranchBasis.");
        }
        return (OrthogonalBlockPreparedBranchBasis) prepared;
    }

}
