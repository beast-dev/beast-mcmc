package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
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
        expSolver.computeCompressed(prepared.blockDParams, dt, prepared.expD);
        OrthogonalBlockMatrixOps.multiplyRightBlockDiagonal(
                prepared.rMatrix.data,
                prepared.expD,
                dimension,
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes(),
                prepared.workMatrix);
        MatrixOps.matMul(
                prepared.workMatrix,
                prepared.rtMatrix.data,
                prepared.transitionMatrix.data,
                dimension);
        prepared.covariancePrepared = false;
    }

    static void fillTransitionMatrix(final double[] rData,
                                     final double[] expDData,
                                     final double[] rtData,
                                     final int dimension,
                                     final int[] blockStarts,
                                     final int[] blockSizes,
                                     final double[] workMatrix,
                                     final double[] transitionData) {
        OrthogonalBlockMatrixOps.multiplyRightBlockDiagonal(
                rData, expDData, dimension, blockStarts, blockSizes, workMatrix);
        MatrixOps.matMul(workMatrix, rtData, transitionData, dimension);
    }
}
