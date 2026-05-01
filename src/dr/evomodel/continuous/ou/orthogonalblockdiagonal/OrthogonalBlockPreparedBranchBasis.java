package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import org.ejml.data.DenseMatrix64F;

/**
 * Reusable branch-local basis prepared for one orthogonal block-diagonal OU transition.
 */
public final class OrthogonalBlockPreparedBranchBasis {
    final int dimension;
    double dt;
    final double[] stationaryMean;
    final double[] blockDParams;
    final DenseMatrix64F expD;
    final DenseMatrix64F rMatrix;
    final DenseMatrix64F rtMatrix;
    final DenseMatrix64F transitionMatrix;
    final DenseMatrix64F qMatrix;
    final DenseMatrix64F stationaryCovDBasis;
    final DenseMatrix64F transitionCovDBasis;
    final DenseMatrix64F transitionCovariance;
    final double[] workMatrix;
    boolean covariancePrepared;

    OrthogonalBlockPreparedBranchBasis(final int dimension,
                                       final int blockDParamDimension) {
        this.dimension = dimension;
        this.dt = Double.NaN;
        this.stationaryMean = new double[dimension];
        this.blockDParams = new double[blockDParamDimension];
        this.expD = new DenseMatrix64F(dimension, dimension);
        this.rMatrix = new DenseMatrix64F(dimension, dimension);
        this.rtMatrix = new DenseMatrix64F(dimension, dimension);
        this.transitionMatrix = new DenseMatrix64F(dimension, dimension);
        this.qMatrix = new DenseMatrix64F(dimension, dimension);
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
