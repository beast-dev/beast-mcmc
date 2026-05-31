package dr.evomodel.continuous.ou.blockdiagonal;

import dr.inference.model.BlockDiagonalDecomposition;
import org.ejml.data.DenseMatrix64F;

/**
 * Reusable branch-local basis prepared for one block-diagonal OU transition.
 */
public class BlockDiagonalPreparedBranchBasis {
    final int dimension;
    final BlockDiagonalDecomposition decomposition;
    double dt;
    final double[] stationaryMean;
    final double[] blockDParams;
    final double[] expD;
    final DenseMatrix64F rMatrix;
    final DenseMatrix64F rinvMatrix;
    /**
     * Historical alias used by orthogonal callers. For a general decomposition this is R^{-1}.
     */
    final DenseMatrix64F rtMatrix;
    final DenseMatrix64F transitionMatrix;
    final DenseMatrix64F qMatrix;
    final DenseMatrix64F qDBasis;
    final DenseMatrix64F stationaryCovDBasis;
    final DenseMatrix64F transitionCovDBasis;
    final DenseMatrix64F transitionCovariance;
    final double[] workMatrix;
    boolean covariancePrepared;

    BlockDiagonalPreparedBranchBasis(final int dimension,
                                     final int compressedBlockDimension,
                                     final BlockDiagonalDecomposition decomposition) {
        this.dimension = dimension;
        if (decomposition.getDimension() != dimension) {
            throw new IllegalArgumentException(
                    "decomposition dimension must be " + dimension + " but is " + decomposition.getDimension());
        }
        this.decomposition = decomposition;
        this.dt = Double.NaN;
        this.stationaryMean = new double[dimension];
        this.blockDParams = decomposition.getBlockDiagonal();
        this.expD = new double[compressedBlockDimension];
        this.rMatrix = DenseMatrix64F.wrap(dimension, dimension, decomposition.getR());
        this.rinvMatrix = DenseMatrix64F.wrap(dimension, dimension, decomposition.getRInverse());
        this.rtMatrix = this.rinvMatrix;
        this.transitionMatrix = new DenseMatrix64F(dimension, dimension);
        this.qMatrix = new DenseMatrix64F(dimension, dimension);
        this.qDBasis = new DenseMatrix64F(dimension, dimension);
        this.stationaryCovDBasis = new DenseMatrix64F(dimension, dimension);
        this.transitionCovDBasis = new DenseMatrix64F(dimension, dimension);
        this.transitionCovariance = new DenseMatrix64F(dimension, dimension);
        this.workMatrix = new double[dimension * dimension];
        this.covariancePrepared = false;
    }

    public void invalidateCovariance() {
        covariancePrepared = false;
    }
}
