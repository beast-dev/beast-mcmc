package dr.evomodel.continuous.ou.diagonal;

import dr.evomodel.continuous.ou.AbstractSelectionMatrixParameterization;
import dr.evomodel.continuous.ou.kernel.OUBranchKernel;
import dr.evomodel.continuous.ou.kernel.OUBranchKernelCapability;
import dr.inference.model.DiagonalMatrix;

import java.util.Arrays;

/**
 * Diagonal OU selection parametrization.
 */
public final class DiagonalSelectionMatrixParameterization
        extends AbstractSelectionMatrixParameterization
        implements OUBranchKernelCapability {

    private final DiagonalMatrix diagonalMatrix;

    public DiagonalSelectionMatrixParameterization(final DiagonalMatrix diagonalMatrix) {
        super(diagonalMatrix);
        this.diagonalMatrix = diagonalMatrix;
    }

    @Override
    public void fillSelectionMatrixFlat(final double[] out) {
        checkFlatSquare(out, "selection matrix");
        Arrays.fill(out, 0.0);
        for (int i = 0; i < getDimension(); ++i) {
            out[i * getDimension() + i] = diagonalMatrix.getParameterValue(i, i);
        }
    }

    @Override
    public void fillTransitionMatrixFlat(final double dt, final double[] out) {
        checkFlatSquare(out, "transition matrix");
        Arrays.fill(out, 0.0);
        for (int i = 0; i < getDimension(); ++i) {
            out[i * getDimension() + i] =
                    Math.exp(-dt * diagonalMatrix.getParameterValue(i, i));
        }
    }

    @Override
    public void accumulateGradientFromTransitionFlat(final double dt,
                                                     final double[] stationaryMean,
                                                     final double[] dLogL_dF,
                                                     final double[] dLogL_df,
                                                     final double[] gradientAccumulator) {
        checkVector(stationaryMean, "stationary mean");
        checkFlatSquare(dLogL_dF, "transition adjoint");
        checkVector(dLogL_df, "transition offset adjoint");
        checkFlatSquare(gradientAccumulator, "selection gradient accumulator");
        final int dimension = getDimension();
        for (int i = 0; i < dimension; ++i) {
            final int diagonalOffset = i * dimension + i;
            final double actualization =
                    Math.exp(-dt * diagonalMatrix.getParameterValue(i, i));
            final double upstream =
                    dLogL_dF[diagonalOffset] - dLogL_df[i] * stationaryMean[i];
            gradientAccumulator[diagonalOffset] += -dt * actualization * upstream;
        }
    }

    @Override
    public OUBranchKernel createBranchKernel() {
        return new DiagonalOUBranchKernel(diagonalMatrix);
    }
}
