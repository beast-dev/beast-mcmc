package dr.inference.timeseries.gaussian;

import dr.inference.model.MatrixParameterInterface;

/**
 * Default dense selection-matrix parametrization.
 */
public class DenseSelectionMatrixParameterization implements SelectionMatrixParameterization {

    private final MatrixParameterInterface matrixParameter;
    private final int dimension;

    public DenseSelectionMatrixParameterization(final MatrixParameterInterface matrixParameter) {
        if (matrixParameter == null) {
            throw new IllegalArgumentException("matrixParameter must not be null");
        }
        if (matrixParameter.getRowDimension() != matrixParameter.getColumnDimension()) {
            throw new IllegalArgumentException("selection matrix must be square");
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
    public void fillSelectionMatrix(final double[][] out) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i][j] = matrixParameter.getParameterValue(i, j);
            }
        }
    }

    @Override
    public void fillTransitionMatrix(final double dt, final double[][] out) {
        final double[][] minusAdt = new double[dimension][dimension];
        fillScaledNegativeSelectionMatrix(dt, minusAdt);
        MatrixExponentialUtils.expm(minusAdt, out);
    }

    @Override
    public void accumulateGradientFromTransition(final double dt,
                                                 final double[] stationaryMean,
                                                 final double[][] dLogL_dF,
                                                 final double[] dLogL_df,
                                                 final double[] gradientAccumulator) {
        final double[][] minusAdt = new double[dimension][dimension];
        final double[][] totalUpstreamOnF = new double[dimension][dimension];
        final double[][] outer = new double[dimension][dimension];
        final double[][] dLogL_dMinusAdt = new double[dimension][dimension];

        fillScaledNegativeSelectionMatrix(dt, minusAdt);
        MatrixExponentialUtils.outerProduct(dLogL_df, stationaryMean, outer);
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                totalUpstreamOnF[i][j] = dLogL_dF[i][j] - outer[i][j];
            }
        }

        MatrixExponentialUtils.adjointExp(minusAdt, totalUpstreamOnF, dLogL_dMinusAdt);
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                gradientAccumulator[i * dimension + j] += -dt * dLogL_dMinusAdt[i][j];
            }
        }
    }

    protected final void fillScaledNegativeSelectionMatrix(final double dt, final double[][] out) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i][j] = -dt * matrixParameter.getParameterValue(i, j);
            }
        }
    }
}
