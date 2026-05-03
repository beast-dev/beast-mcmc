package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovSolver;
import dr.inference.model.MatrixParameterInterface;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

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
        CommonOps.mult(rtMatrix, qMatrix, temp);
        CommonOps.mult(temp, rMatrix, qDBasis);
        lyapunovSolver.solve(blockDParams, qDBasis, stationaryCovDBasis);
        OrthogonalBlockMatrixOps.multiplyBlockDiagonalLeft(
                expD, stationaryCovDBasis.data, temp.data, temp.numRows, blockStarts, blockSizes);
        OrthogonalBlockMatrixOps.multiplyRightBlockDiagonalTranspose(
                temp.data, expD, transitionCovDBasis.data, transitionCovDBasis.numRows, blockStarts, blockSizes);
        CommonOps.subtract(stationaryCovDBasis, transitionCovDBasis, transitionCovDBasis);
        if (symmetrizeTransitionCovDBasis) {
            symmetrize(transitionCovDBasis);
        }
        CommonOps.mult(rMatrix, transitionCovDBasis, temp);
        if (symmetrizeTransitionCovDBasis) {
            multiplyRightTransposeSymmetric(temp, rMatrix, out);
        } else {
            CommonOps.mult(temp, rtMatrix, out);
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

    private static void multiplyRightTransposeSymmetric(final DenseMatrix64F left,
                                                        final DenseMatrix64F right,
                                                        final DenseMatrix64F out) {
        final int dimension = left.numRows;
        final double[] leftData = left.data;
        final double[] rightData = right.data;
        final double[] outData = out.data;
        for (int i = 0; i < dimension; ++i) {
            final int leftRowOffset = i * dimension;
            for (int j = i; j < dimension; ++j) {
                final int rightRowOffset = j * dimension;
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += leftData[leftRowOffset + k] * rightData[rightRowOffset + k];
                }
                outData[i * dimension + j] = sum;
                outData[j * dimension + i] = sum;
            }
        }
    }

}
