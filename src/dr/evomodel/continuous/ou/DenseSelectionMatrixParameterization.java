package dr.evomodel.continuous.ou;

import dr.inference.model.MatrixParameterInterface;

/**
 * Default dense selection-matrix parametrization.
 */
public class DenseSelectionMatrixParameterization extends AbstractSelectionMatrixParameterization {

    private final ThreadLocal<Workspace> workspaceLocal;

    public DenseSelectionMatrixParameterization(final MatrixParameterInterface matrixParameter) {
        super(matrixParameter);
        this.workspaceLocal = ThreadLocal.withInitial(() -> new Workspace(getDimension()));
    }

    @Override
    public void fillSelectionMatrixFlat(final double[] out) {
        checkFlatSquare(out, "selection matrix");
        final MatrixParameterInterface matrixParameter = getMatrixParameter();
        final int dimension = getDimension();
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
        MatrixExponentialUtils.expmFlat(workspace.minusAdt, out, getDimension());
    }

    @Override
    public void accumulateGradientFromTransitionFlat(final double dt,
                                                     final double[] stationaryMean,
                                                     final double[] dLogL_dF,
                                                     final double[] dLogL_df,
                                                     final double[] gradientAccumulator) {
        checkFlatSquare(dLogL_dF, "transition adjoint");
        checkVector(stationaryMean, "stationary mean");
        checkVector(dLogL_df, "transition offset adjoint");
        checkFlatSquare(gradientAccumulator, "selection gradient accumulator");
        final Workspace workspace = workspace();
        final int dimension = getDimension();

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
        final MatrixParameterInterface matrixParameter = getMatrixParameter();
        final int dimension = getDimension();
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
