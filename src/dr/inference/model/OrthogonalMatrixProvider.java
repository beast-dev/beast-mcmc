package dr.inference.model;

/**
 * Marker interface for matrix parameters that are constrained to be orthogonal
 * and are backed by a lower-dimensional chart parameterization.
 */
public interface OrthogonalMatrixProvider {

    Parameter getOrthogonalParameter();

    void fillOrthogonalMatrix(double[] outRowMajor);

    void fillOrthogonalTranspose(double[] outRowMajor);

    void applyOrthogonal(double[] in, double[] out);

    void applyOrthogonalTranspose(double[] in, double[] out);

    double[] pullBackGradient(double[][] gradientWrtMatrix);

    default double[] pullBackGradientFlat(final double[] gradientWrtMatrixRowMajor,
                                          final int dimension) {
        if (gradientWrtMatrixRowMajor == null || gradientWrtMatrixRowMajor.length != dimension * dimension) {
            throw new IllegalArgumentException(
                    "gradientWrtMatrixRowMajor must have length " + (dimension * dimension));
        }
        final double[][] gradientWrtMatrix = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(gradientWrtMatrixRowMajor, i * dimension, gradientWrtMatrix[i], 0, dimension);
        }
        return pullBackGradient(gradientWrtMatrix);
    }
}
