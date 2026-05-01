package dr.evomodel.treedatalikelihood.continuous.canonical;

import java.util.Arrays;

final class DenseGradientWorkspace {

    final double[][] covariance2;
    final double[] transitionMatrixFlat;
    final double[] covarianceAdjointFlat;
    final double[] matrixProductFlat;
    final double[] localGradientA;
    final double[] localGradientQ;
    final double[] localGradientMuVector;
    final double[] localGradientMuScalar;

    private double cachedTransitionMatrixLength;
    private long transitionMatrixCacheEpoch;
    private long activeTransitionMatrixCacheEpoch;
    private boolean cachedTransitionMatrixValid;

    DenseGradientWorkspace(final int dim) {
        if (dim < 1) throw new IllegalArgumentException("dim must be >= 1");
        this.covariance2 = new double[dim][dim];
        this.transitionMatrixFlat = new double[dim * dim];
        this.covarianceAdjointFlat = new double[dim * dim];
        this.matrixProductFlat = new double[dim * dim];
        this.localGradientA = new double[dim * dim];
        this.localGradientQ = new double[dim * dim];
        this.localGradientMuVector = new double[dim];
        this.localGradientMuScalar = new double[1];
        this.cachedTransitionMatrixLength = Double.NaN;
        this.transitionMatrixCacheEpoch = 0L;
        this.activeTransitionMatrixCacheEpoch = -1L;
        this.cachedTransitionMatrixValid = false;
    }

    double[] localGradientMu(final int gradientLength, final int dim) {
        if (gradientLength == 1) {
            return localGradientMuScalar;
        }
        if (gradientLength == dim) {
            return localGradientMuVector;
        }
        throw new IllegalArgumentException(
                "Stationary-mean gradient length must be 1 or " + dim + ", found " + gradientLength);
    }

    void clearLocalGradientBuffers(final int gradALength, final int gradMuLength, final int dim) {
        Arrays.fill(localGradientA, 0, gradALength, 0.0);
        Arrays.fill(localGradientQ, 0.0);
        Arrays.fill(localGradientMu(gradMuLength, dim), 0.0);
        invalidateTransitionMatrixCache();
    }

    void invalidateTransitionMatrixCache() {
        cachedTransitionMatrixValid = false;
        activeTransitionMatrixCacheEpoch = -1L;
        transitionMatrixCacheEpoch++;
    }

    boolean hasTransitionMatrix(final double branchLength) {
        return cachedTransitionMatrixValid
                && activeTransitionMatrixCacheEpoch == transitionMatrixCacheEpoch
                && Double.doubleToLongBits(cachedTransitionMatrixLength)
                == Double.doubleToLongBits(branchLength);
    }

    void cacheTransitionMatrix(final double branchLength) {
        cachedTransitionMatrixLength = branchLength;
        activeTransitionMatrixCacheEpoch = transitionMatrixCacheEpoch;
        cachedTransitionMatrixValid = true;
    }
}
