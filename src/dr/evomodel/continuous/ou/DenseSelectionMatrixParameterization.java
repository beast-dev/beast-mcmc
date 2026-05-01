package dr.evomodel.continuous.ou;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
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
    public void fillSelectionMatrix(final double[][] out) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i][j] = matrixParameter.getParameterValue(i, j);
            }
        }
    }

    @Override
    public void fillTransitionMatrix(final double dt, final double[][] out) {
        final Workspace workspace = workspace();
        fillScaledNegativeSelectionMatrix(dt, workspace.minusAdt);
        MatrixExponentialUtils.expm(workspace.minusAdt, out);
    }

    @Override
    public void fillTransitionMatrixFlat(final double dt, final double[] out) {
        checkFlatSquare(out, "transition matrix");
        final Workspace workspace = workspace();
        fillTransitionMatrix(dt, workspace.transitionMatrix);
        MatrixOps.toFlat(workspace.transitionMatrix, out, dimension);
    }

    @Override
    public void accumulateGradientFromTransition(final double dt,
                                                 final double[] stationaryMean,
                                                 final double[][] dLogL_dF,
                                                 final double[] dLogL_df,
                                                 final double[] gradientAccumulator) {
        final Workspace workspace = workspace();

        fillScaledNegativeSelectionMatrix(dt, workspace.minusAdt);
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                workspace.totalUpstreamOnF[i][j] = dLogL_dF[i][j] - dLogL_df[i] * stationaryMean[j];
            }
        }

        MatrixExponentialUtils.adjointExp(
                workspace.minusAdt,
                workspace.totalUpstreamOnF,
                workspace.dLogL_dMinusAdt);
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                gradientAccumulator[i * dimension + j] += -dt * workspace.dLogL_dMinusAdt[i][j];
            }
        }
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
                workspace.totalUpstreamOnF[i][j] =
                        dLogL_dF[rowOffset + j] - offsetAdjoint * stationaryMean[j];
            }
        }

        MatrixExponentialUtils.adjointExp(
                workspace.minusAdt,
                workspace.totalUpstreamOnF,
                workspace.dLogL_dMinusAdt);
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                gradientAccumulator[rowOffset + j] += -dt * workspace.dLogL_dMinusAdt[i][j];
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
        private final double[][] minusAdt;
        private final double[][] transitionMatrix;
        private final double[][] totalUpstreamOnF;
        private final double[][] dLogL_dMinusAdt;

        private Workspace(final int dimension) {
            this.minusAdt = new double[dimension][dimension];
            this.transitionMatrix = new double[dimension][dimension];
            this.totalUpstreamOnF = new double[dimension][dimension];
            this.dLogL_dMinusAdt = new double[dimension][dimension];
        }
    }
}
