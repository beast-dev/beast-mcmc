package dr.evomodel.continuous.ou;

import dr.inference.model.MatrixParameterInterface;

/**
 * Default dense selection-matrix parametrization.
 */
public class DenseSelectionMatrixParameterization implements SelectionMatrixParameterization {

    private final MatrixParameterInterface matrixParameter;
    private final int dimension;
    private final ThreadLocal<Workspace> workspaceLocal;

    public DenseSelectionMatrixParameterization(final MatrixParameterInterface matrixParameter) {
        if (matrixParameter == null) {
            throw new IllegalArgumentException("matrixParameter must not be null");
        }
        if (matrixParameter.getRowDimension() != matrixParameter.getColumnDimension()) {
            throw new IllegalArgumentException("selection matrix must be square");
        }
        this.matrixParameter = matrixParameter;
        this.dimension = matrixParameter.getRowDimension();
        this.workspaceLocal = ThreadLocal.withInitial(() -> new Workspace(dimension));
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
    public void fillSelectionMatrixFlat(final double[] out) {
        checkFlatSquare(out, "selection matrix");
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                out[rowOffset + j] = matrixParameter.getParameterValue(i, j);
            }
        }
    }

    @Override
    public void fillTransitionMatrixFlat(final double dt, final double[] out) {
        checkFlatSquare(out, "transition matrix");
        final Workspace workspace = workspace();
        fillScaledNegativeSelectionMatrix(dt, workspace.minusAdt);
        MatrixExponentialUtils.expmFlat(workspace.minusAdt, out, dimension);
    }

    @Override
    public void accumulateGradientFromTransitionFlat(final double dt,
                                                     final double[] stationaryMean,
                                                     final double[] dLogL_dF,
                                                     final double[] dLogL_df,
                                                     final double[] gradientAccumulator) {
        checkFlatSquare(dLogL_dF, "transition adjoint");
        final Workspace workspace = workspace();

        fillScaledNegativeSelectionMatrix(dt, workspace.minusAdt);
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            final double offsetAdjoint = dLogL_df[i];
            for (int j = 0; j < dimension; ++j) {
                workspace.totalUpstreamOnF[rowOffset + j] =
                        dLogL_dF[rowOffset + j] - offsetAdjoint * stationaryMean[j];
            }
        }

        MatrixExponentialUtils.adjointExpFlat(
                workspace.minusAdt,
                workspace.totalUpstreamOnF,
                workspace.dLogL_dMinusAdt,
                dimension);
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                gradientAccumulator[rowOffset + j] += -dt * workspace.dLogL_dMinusAdt[rowOffset + j];
            }
        }
    }

    protected final void fillScaledNegativeSelectionMatrix(final double dt, final double[] out) {
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                out[rowOffset + j] = -dt * matrixParameter.getParameterValue(i, j);
            }
        }
    }

    private Workspace workspace() {
        return workspaceLocal.get();
    }

    private void checkFlatSquare(final double[] matrix, final String label) {
        if (matrix == null || matrix.length != dimension * dimension) {
            throw new IllegalArgumentException(
                    label + " must have length " + (dimension * dimension)
                            + " for a " + dimension + "x" + dimension + " matrix");
        }
    }

    private static final class Workspace {
        private final double[] minusAdt;
        private final double[] totalUpstreamOnF;
        private final double[] dLogL_dMinusAdt;

        private Workspace(final int dimension) {
            final int length = dimension * dimension;
            this.minusAdt = new double[length];
            this.totalUpstreamOnF = new double[length];
            this.dLogL_dMinusAdt = new double[length];
        }
    }
}
