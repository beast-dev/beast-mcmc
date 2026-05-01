package dr.evomodel.continuous.ou;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;

/**
 * Parametrization-specific access to the OU diffusion matrix Q and pullback of
 * analytical gradients from dense matrix space to native parameters.
 */
public interface DiffusionMatrixParameterization {

    MatrixParameterInterface getMatrixParameter();

    int getDimension();

    void fillDiffusionMatrix(double[][] out);

    default void fillDiffusionMatrixFlat(final double[] out) {
        final int dimension = getDimension();
        final double[][] dense = new double[dimension][dimension];
        fillDiffusionMatrix(dense);
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(dense[i], 0, out, i * dimension, dimension);
        }
    }

    boolean supportsParameter(Parameter parameter);

    double[] pullBackGradient(Parameter parameter, double[] denseGradient);
}
