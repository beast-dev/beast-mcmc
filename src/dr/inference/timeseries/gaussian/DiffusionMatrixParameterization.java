package dr.inference.timeseries.gaussian;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;

/**
 * Parametrization-specific access to the diffusion matrix Q and pullback of
 * analytical gradients from the dense matrix space to native parameters.
 */
public interface DiffusionMatrixParameterization {

    MatrixParameterInterface getMatrixParameter();

    int getDimension();

    void fillDiffusionMatrix(double[][] out);

    boolean supportsParameter(Parameter parameter);

    double[] pullBackGradient(Parameter parameter, double[] denseGradient);
}
