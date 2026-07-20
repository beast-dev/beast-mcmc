package dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;

import java.util.Arrays;

/**
 * Orthogonal specialization of the block-diagonal parametrization.
 *
 * <p>All transition/covariance/diffusion sandwich adjoints are inherited unmodified
 * from {@link BlockDiagonalSelectionMatrixParameterization}: they consume only
 * {@code rMatrix}/{@code rinvMatrix} as dense arrays, and
 * {@code AbstractBlockDiagonalTwoByTwoMatrixParameter.fillRAndRinv} already fills
 * {@code rinvMatrix} as exactly {@code R^T} (via
 * {@code InvertibleMatrixParametrization.fillInverse}) for an orthogonal rotation, so
 * the inherited sandwich math is numerically identical to a separate orthogonal-specific
 * implementation — confirmed by finite-difference cross-checks against the general
 * backend in {@code BlockDiagonalGeneralVsOrthogonalGradientCrossCheckTest}.
 *
 * <p>What this class actually supplies:
 * <ul>
 *   <li>Native-parameter bookkeeping ({@link #getSelectionGradientDimension},
 *       {@link #finishNativeSelectionGradient}, {@link #accumulateGradientFromTransitionFlat}):
 *       the orthogonal rotation is parametrized by a small number of Givens angles, not
 *       {@code dimension^2} free dense entries, so gradients must be pulled back through
 *       {@link OrthogonalMatrixProvider#fillPullBackGradientFlat} to respect that
 *       lower-dimensional manifold — this is a genuine difference, not duplication.</li>
 *   <li>The stationary-mean-gradient path ({@link #accumulateMeanGradientPreparedDBasisFlat},
 *       {@link #finishMeanGradientFromDBasisFlat}, via {@link OrthogonalBlockSelectionAdjoint}):
 *       composes with {@link OrthogonalMatrixProvider}'s vector-apply operations rather than
 *       a dense rotation multiply.</li>
 * </ul>
 */
public final class OrthogonalBlockDiagonalSelectionMatrixParameterization
        extends BlockDiagonalSelectionMatrixParameterization
        implements OrthogonalBlockCanonicalParameterization {

    private final OrthogonalMatrixProvider orthogonalRotation;
    private final OrthogonalBlockSelectionAdjoint selectionAdjoint;
    private final double[] nativeBlockGradientScratch;

    public OrthogonalBlockDiagonalSelectionMatrixParameterization(
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final OrthogonalMatrixProvider orthogonalRotation) {
        super(blockParameter);
        this.orthogonalRotation = orthogonalRotation;
        this.selectionAdjoint = new OrthogonalBlockSelectionAdjoint(blockParameter, orthogonalRotation);
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
