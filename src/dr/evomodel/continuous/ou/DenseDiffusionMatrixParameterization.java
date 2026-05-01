package dr.evomodel.continuous.ou;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;

/**
 * Default dense OU diffusion-matrix parametrization.
 */
public final class DenseDiffusionMatrixParameterization implements DiffusionMatrixParameterization {

    private final MatrixParameterInterface matrixParameter;
    private final int dimension;

    public DenseDiffusionMatrixParameterization(final MatrixParameterInterface matrixParameter) {
        if (matrixParameter == null) {
            throw new IllegalArgumentException("matrixParameter must not be null");
        }
        if (matrixParameter.getRowDimension() != matrixParameter.getColumnDimension()) {
            throw new IllegalArgumentException("diffusion matrix must be square");
        }
        this.matrixParameter = matrixParameter;
        this.dimension = matrixParameter.getRowDimension();
    }

    @Override
    public MatrixParameterInterface getMatrixParameter() {
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
        return parameter == matrixParameter;
    }

    @Override
    public double[] pullBackGradient(final Parameter parameter, final double[] denseGradient) {
        if (parameter != matrixParameter) {
            throw new IllegalArgumentException("Unsupported dense diffusion parameter");
        }
        return denseGradient;
    }
}
