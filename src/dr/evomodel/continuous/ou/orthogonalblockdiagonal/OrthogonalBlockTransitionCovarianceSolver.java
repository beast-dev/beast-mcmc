package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovSolver;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.model.MatrixParameterInterface;
import org.ejml.data.DenseMatrix64F;

final class OrthogonalBlockTransitionCovarianceSolver {

    private OrthogonalBlockTransitionCovarianceSolver() { }

    static void fillTransitionCovariance(final MatrixParameterInterface diffusionMatrix,
                                         final DenseMatrix64F rMatrix,
                                         final DenseMatrix64F rtMatrix,
                                         final double[] expD,
                                         final double[] blockDParams,
                                         final BlockDiagonalLyapunovSolver lyapunovSolver,
                                         final DenseMatrix64F qMatrix,
                                         final DenseMatrix64F qDBasis,
                                         final DenseMatrix64F stationaryCovDBasis,
                                         final DenseMatrix64F transitionCovDBasis,
                                         final DenseMatrix64F temp,
                                         final DenseMatrix64F out,
                                         final boolean symmetrizeTransitionCovDBasis,
                                         final int[] blockStarts,
                                         final int[] blockSizes) {
        fillDenseMatrix(diffusionMatrix, qMatrix);
        MatrixOps.symmetricSandwichTransposeLeft(rMatrix.data, qMatrix.data, qDBasis.data, temp.data, qMatrix.numRows);
        lyapunovSolver.solve(blockDParams, qDBasis, stationaryCovDBasis);
        OrthogonalBlockMatrixOps.multiplyBlockDiagonalLeft(
                expD, stationaryCovDBasis.data, temp.data, temp.numRows, blockStarts, blockSizes);
        OrthogonalBlockMatrixOps.multiplyRightBlockDiagonalTranspose(
                temp.data, expD, transitionCovDBasis.data, transitionCovDBasis.numRows, blockStarts, blockSizes);
        OrthogonalBlockDenseMatrixOps.subtract(stationaryCovDBasis, transitionCovDBasis, transitionCovDBasis);
        if (symmetrizeTransitionCovDBasis) {
            symmetrize(transitionCovDBasis);
        }
        if (symmetrizeTransitionCovDBasis) {
            MatrixOps.symmetricSandwichTransposeRight(
                    rMatrix.data,
                    transitionCovDBasis.data,
                    out.data,
                    temp.data,
                    out.numRows);
        } else {
            OrthogonalBlockDenseMatrixOps.mult(rMatrix, transitionCovDBasis, temp);
            OrthogonalBlockDenseMatrixOps.mult(temp, rtMatrix, out);
            symmetrize(out);
        }
    }

    static void copyPreparedCovariance(final OrthogonalBlockPreparedBranchBasis prepared,
                                       final OrthogonalBlockBranchGradientWorkspace workspace) {
        workspace.qMatrix.set(prepared.qMatrix);
        workspace.stationaryCovDBasis.set(prepared.stationaryCovDBasis);
        workspace.transitionCovDBasis.set(prepared.transitionCovDBasis);
        workspace.transitionCovariance.set(prepared.transitionCovariance);
    }

    static void storePreparedCovariance(final OrthogonalBlockPreparedBranchBasis prepared,
                                        final OrthogonalBlockBranchGradientWorkspace workspace) {
        prepared.qMatrix.set(workspace.qMatrix);
        prepared.stationaryCovDBasis.set(workspace.stationaryCovDBasis);
        prepared.transitionCovDBasis.set(workspace.transitionCovDBasis);
        prepared.transitionCovariance.set(workspace.transitionCovariance);
        prepared.covariancePrepared = true;
    }

    private static void fillDenseMatrix(final MatrixParameterInterface parameter, final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                data[i * dimension + j] = parameter.getParameterValue(i, j);
            }
        }
    }

    private static void symmetrize(final DenseMatrix64F matrix) {
        final int dimension = matrix.numRows;
        final double[] data = matrix.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final double value = 0.5 * (data[i * dimension + j] + data[j * dimension + i]);
                data[i * dimension + j] = value;
                data[j * dimension + i] = value;
            }
        }
    }

}
