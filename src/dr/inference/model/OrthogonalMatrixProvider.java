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
}
