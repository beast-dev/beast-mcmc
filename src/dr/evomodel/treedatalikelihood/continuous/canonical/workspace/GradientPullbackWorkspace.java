package dr.evomodel.treedatalikelihood.continuous.canonical.workspace;

import java.util.EnumSet;

public final class GradientPullbackWorkspace {
    public final DenseGradientWorkspace dense;
    public final OrthogonalBlockGradientWorkspace orthogonal;

    public final double[] orthogonalStationaryMeanScratch;
    public final double[] orthogonalCompressedGradientScratch;
    public final double[] orthogonalNativeGradientScratch;
    public final double[] orthogonalRotationGradientFlatScratch;
    public final double[][] orthogonalRotationGradientScratch;
    public final double[][] covariance2;
    public final double[] transitionMatrixFlat;
    public final double[] covarianceAdjointFlat;
    public final double[] matrixProductFlat;
    public final double[] localGradientA;
    public final double[] localGradientQ;
    public final double[] localGradientMuVector;
    public final double[] localGradientMuScalar;

    public GradientPullbackWorkspace(final int dim) {
        this(dim, EnumSet.of(
                WorkspaceCapability.DENSE_GRADIENT,
                WorkspaceCapability.ORTHOGONAL_BLOCK_GRADIENT));
    }

    public GradientPullbackWorkspace(final int dim, final EnumSet<WorkspaceCapability> capabilities) {
        if (dim < 1) throw new IllegalArgumentException("dim must be >= 1");
        this.dense = capabilities.contains(WorkspaceCapability.DENSE_GRADIENT)
                ? new DenseGradientWorkspace(dim)
                : null;
        this.orthogonal = capabilities.contains(WorkspaceCapability.ORTHOGONAL_BLOCK_GRADIENT)
                ? new OrthogonalBlockGradientWorkspace(dim)
                : null;

        this.orthogonalStationaryMeanScratch =
                orthogonal == null ? null : orthogonal.stationaryMeanScratch;
        this.orthogonalCompressedGradientScratch =
                orthogonal == null ? null : orthogonal.compressedGradientScratch;
        this.orthogonalNativeGradientScratch =
                orthogonal == null ? null : orthogonal.nativeGradientScratch;
        this.orthogonalRotationGradientFlatScratch =
                orthogonal == null ? null : orthogonal.rotationGradientFlatScratch;
        this.orthogonalRotationGradientScratch =
                orthogonal == null ? null : orthogonal.rotationGradientScratch;

        this.covariance2 = dense == null ? null : dense.covariance2;
        this.transitionMatrixFlat = dense == null ? null : dense.transitionMatrixFlat;
        this.covarianceAdjointFlat = dense == null ? null : dense.covarianceAdjointFlat;
        this.matrixProductFlat = dense == null ? null : dense.matrixProductFlat;
        this.localGradientA = dense == null ? null : dense.localGradientA;
        this.localGradientQ = dense == null ? null : dense.localGradientQ;
        this.localGradientMuVector = dense == null ? null : dense.localGradientMuVector;
        this.localGradientMuScalar = dense == null ? null : dense.localGradientMuScalar;
    }

    public double[] localGradientMu(final int gradientLength, final int dim) {
        return dense().localGradientMu(gradientLength, dim);
    }

    /**
     * Clears the local gradient accumulators and resets the transition-matrix cache.
     * Specialized (orthogonal-block) fields are NOT cleared here; that is the
     * responsibility of {@link SpecializedCanonicalSelectionGradientPullback#clearWorkerBuffers}.
     */
    public void clearLocalGradientBuffers(final int gradALength, final int gradMuLength, final int dim) {
        dense().clearLocalGradientBuffers(gradALength, gradMuLength, dim);
    }

    public void clearLocalGradientBuffers(final int gradALength,
                                   final int gradMuLength,
                                   final int dim,
                                   final boolean orthogonalSelection,
                                   final int compressedGradientLength) {
        clearLocalGradientBuffers(gradALength, gradMuLength, dim);
        if (orthogonalSelection) {
            orthogonal().clearSpecializedBuffers(compressedGradientLength);
        }
    }

    public void invalidateTransitionMatrixCache() {
        dense().invalidateTransitionMatrixCache();
    }

    public boolean hasTransitionMatrix(final double branchLength) {
        return dense().hasTransitionMatrix(branchLength);
    }

    public void cacheTransitionMatrix(final double branchLength) {
        dense().cacheTransitionMatrix(branchLength);
    }

    private DenseGradientWorkspace dense() {
        if (dense == null) {
            throw new IllegalStateException("Dense gradient workspace capability was not requested");
        }
        return dense;
    }

    public OrthogonalBlockGradientWorkspace orthogonal() {
        if (orthogonal == null) {
            throw new IllegalStateException("Orthogonal-block gradient workspace capability was not requested");
        }
        return orthogonal;
    }
}
