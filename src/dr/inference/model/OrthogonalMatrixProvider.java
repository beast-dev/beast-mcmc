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

    default void fillPullBackGradientFlat(final double[] gradientWrtMatrixRowMajor,
                                          final int dimension,
                                          final double[] out) {
        fillPullBackGradientFlat(gradientWrtMatrixRowMajor, dimension, out, 0);
    }

    default void fillPullBackGradientFlat(final double[] gradientWrtMatrixRowMajor,
                                          final int dimension,
                                          final double[] out,
                                          final int offset) {
        if (gradientWrtMatrixRowMajor == null || gradientWrtMatrixRowMajor.length != dimension * dimension) {
            throw new IllegalArgumentException(
                    "gradientWrtMatrixRowMajor must have length " + (dimension * dimension));
        }
        final int angleCount = getOrthogonalParameter().getDimension();
        if (out == null || out.length < offset + angleCount) {
            throw new IllegalArgumentException("out must have room for " + angleCount + " entries");
        }
        final double[][] gradientWrtMatrix = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(gradientWrtMatrixRowMajor, i * dimension, gradientWrtMatrix[i], 0, dimension);
        }
        final double[] gradient = pullBackGradient(gradientWrtMatrix);
        System.arraycopy(gradient, 0, out, offset, angleCount);
    }
}
