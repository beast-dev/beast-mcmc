package dr.evomodel.treedatalikelihood.continuous.canonical.adapter;

import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.model.MatrixParameter;
import org.ejml.data.DenseMatrix64F;

/**
 * Maintains the covariance snapshot used by canonical OU process code when
 * the BEAST-facing diffusion model exposes precision.
 */
final class CanonicalOUDiffusionSnapshot {

    private final int dimension;
    private final MultivariateDiffusionModel diffusionModel;
    private final MatrixParameter covarianceParameter;
    private final DenseMatrix64F precision;
    private final DenseMatrix64F covariance;
    private final double[] choleskyForwardScratch;
    private final double[] choleskyBackwardScratch;

    CanonicalOUDiffusionSnapshot(final int dimension,
                                 final MultivariateDiffusionModel diffusionModel,
                                 final MatrixParameter covarianceParameter) {
        this.dimension = dimension;
        this.diffusionModel = diffusionModel;
        this.covarianceParameter = covarianceParameter;
        this.precision = new DenseMatrix64F(dimension, dimension);
        this.covariance = new DenseMatrix64F(dimension, dimension);
        this.choleskyForwardScratch = new double[dimension];
        this.choleskyBackwardScratch = new double[dimension];
    }

    void refresh() {
        final double[][] precisionMatrix = diffusionModel.getPrecisionMatrix();
        MatrixOps.toFlat(precisionMatrix, precision.data, dimension);
        try {
            MatrixOps.invertSPDCompact(
                    precision.data,
                    covariance.data,
                    dimension,
                    choleskyForwardScratch,
                    choleskyBackwardScratch);
        } catch (IllegalArgumentException e) {
            final ArithmeticException arithmeticException =
                    new ArithmeticException("Canonical OU precision matrix is not symmetric positive definite.");
            arithmeticException.initCause(e);
            throw arithmeticException;
        }

        final double[] covarianceData = covariance.data;
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                covarianceParameter.setParameterValueQuietly(i, j, covarianceData[i * dimension + j]);
            }
        }
        covarianceParameter.fireParameterChangedEvent();
    }

    void fillCovarianceFlat(final double[] out) {
        MatrixOps.toFlat(covariance, out, dimension);
    }
}
