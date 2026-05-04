package dr.inference.timeseries.gaussian;

import dr.inference.model.CompoundSymmetricMatrix;
import dr.inference.model.Parameter;

/**
 * Diffusion parametrization for covariance matrices represented by a diagonal
 * variance vector together with an off-diagonal correlation (or Cholesky-correlation)
 * vector through {@link CompoundSymmetricMatrix}.
 */
public final class CompoundSymmetricDiffusionMatrixParameterization implements DiffusionMatrixParameterization {

    private final CompoundSymmetricMatrix matrixParameter;
    private final int dimension;

    public CompoundSymmetricDiffusionMatrixParameterization(final CompoundSymmetricMatrix matrixParameter) {
        if (matrixParameter == null) {
            throw new IllegalArgumentException("matrixParameter must not be null");
        }
        this.matrixParameter = matrixParameter;
        this.dimension = matrixParameter.getRowDimension();
    }

    @Override
    public CompoundSymmetricMatrix getMatrixParameter() {
        return matrixParameter;
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public void fillDiffusionMatrix(final double[][] out) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i][j] = matrixParameter.getParameterValue(i, j);
            }
        }
    }

    @Override
    public boolean supportsParameter(final Parameter parameter) {
        return parameter == matrixParameter
                || parameter == matrixParameter.getDiagonalParameter()
                || parameter == matrixParameter.getOffDiagonalParameter()
                || parameter == matrixParameter.getUntransformedOffDiagonalParameter()
                || parameter == matrixParameter.getUntransformedCompoundParameter();
    }

    @Override
    public double[] pullBackGradient(final Parameter parameter, final double[] denseGradient) {
        final double[] diagonalGradient = matrixParameter.updateGradientDiagonal(denseGradient);
        final double[] offDiagonalGradient = matrixParameter.updateGradientOffDiagonal(denseGradient);

        if (parameter == matrixParameter) {
            final double[] out = new double[diagonalGradient.length + offDiagonalGradient.length];
            System.arraycopy(diagonalGradient, 0, out, 0, diagonalGradient.length);
            System.arraycopy(offDiagonalGradient, 0, out, diagonalGradient.length, offDiagonalGradient.length);
            return out;
        }
        if (parameter == matrixParameter.getDiagonalParameter()) {
            return diagonalGradient;
        }
        if (parameter == matrixParameter.getOffDiagonalParameter()
                || parameter == matrixParameter.getUntransformedOffDiagonalParameter()) {
            return offDiagonalGradient;
        }
        if (parameter == matrixParameter.getUntransformedCompoundParameter()) {
            final double[] out = new double[diagonalGradient.length + offDiagonalGradient.length];
            System.arraycopy(diagonalGradient, 0, out, 0, diagonalGradient.length);
            System.arraycopy(offDiagonalGradient, 0, out, diagonalGradient.length, offDiagonalGradient.length);
            return out;
        }
        throw new IllegalArgumentException("Unsupported compound-symmetric diffusion parameter");
    }
}
