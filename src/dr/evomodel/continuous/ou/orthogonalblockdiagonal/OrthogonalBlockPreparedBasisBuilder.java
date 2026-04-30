package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;

final class OrthogonalBlockPreparedBasisBuilder {

    private OrthogonalBlockPreparedBasisBuilder() { }

    static void prepare(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                        final OrthogonalMatrixProvider orthogonalRotation,
                        final BlockDiagonalExpSolver expSolver,
                        final double dt,
                        final double[] stationaryMean,
                        final OrthogonalBlockPreparedBranchBasis prepared) {
        final int dimension = prepared.dimension;
        if (stationaryMean.length != dimension) {
            throw new IllegalArgumentException(
                    "stationaryMean must have length " + dimension + " but has " + stationaryMean.length);
        }

        prepared.dt = dt;
        System.arraycopy(stationaryMean, 0, prepared.stationaryMean, 0, dimension);
        orthogonalRotation.fillOrthogonalMatrix(prepared.rMatrix.data);
        orthogonalRotation.fillOrthogonalTranspose(prepared.rtMatrix.data);
        blockParameter.fillBlockDiagonalElements(prepared.blockDParams);
        expSolver.compute(prepared.blockDParams, dt, prepared.expD);
        multiplyRowMajor(
                prepared.rMatrix.data,
                prepared.expD.data,
                dimension,
                prepared.workMatrix);
        multiply(
                prepared.workMatrix,
                prepared.rtMatrix.data,
                dimension,
                prepared.transitionMatrix.data);
        prepared.covariancePrepared = false;
    }

    static void fillTransitionMatrix(final double[] rData,
                                     final double[] expDData,
                                     final double[] rtData,
                                     final int dimension,
                                     final double[][] workMatrix,
                                     final double[] transitionData) {
        multiplyRowMajor(rData, expDData, dimension, workMatrix);
        multiply(workMatrix, rtData, dimension, transitionData);
    }

    private static void multiplyRowMajor(final double[] leftRowMajor,
                                         final double[] rightRowMajor,
                                         final int dimension,
                                         final double[][] out) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += leftRowMajor[i * dimension + k] * rightRowMajor[k * dimension + j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void multiply(final double[][] left,
                                 final double[] rightRowMajor,
                                 final int dimension,
                                 final double[] outData) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += left[i][k] * rightRowMajor[k * dimension + j];
                }
                outData[i * dimension + j] = sum;
            }
        }
    }
}
