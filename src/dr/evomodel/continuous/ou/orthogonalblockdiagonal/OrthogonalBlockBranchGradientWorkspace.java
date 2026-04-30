package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.continuous.ou.CanonicalBranchWorkspace;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalFrechetHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovAdjointHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovSolver;
import org.ejml.data.DenseMatrix64F;

/**
 * Caller-owned buffers for orthogonal block-diagonal branch gradient pullbacks.
 */
public final class OrthogonalBlockBranchGradientWorkspace implements CanonicalBranchWorkspace {
    final BlockDiagonalFrechetHelper frechetHelper;
    final BlockDiagonalLyapunovSolver lyapunovSolver;
    final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper;
    final DenseMatrix64F qMatrix;
    final DenseMatrix64F qDBasis;
    final DenseMatrix64F stationaryCovDBasis;
    final DenseMatrix64F transitionCovDBasis;
    final DenseMatrix64F transitionCovariance;
    final double[] transitionCovarianceArrayScratch;
    final double[] transitionOffsetScratch;
    final double[][] scaledNegativeBlockDScratch;
    final double[][] denseAdjointScratch;
    final DenseMatrix64F upstreamF;
    final DenseMatrix64F upstreamFD;
    final DenseMatrix64F gradD;
    final DenseMatrix64F gradR;
    final DenseMatrix64F gV;
    final DenseMatrix64F hDBasis;
    final DenseMatrix64F gS;
    final DenseMatrix64F yAdjoint;
    final DenseMatrix64F gECov;
    final DenseMatrix64F temp1;
    final DenseMatrix64F temp2;
    final DenseMatrix64F temp3;
    final double[] choleskyScratch;
    final double[] lowerInverseScratch;
    final double[] tempVector1;
    final double[] tempVector2;

    public OrthogonalBlockBranchGradientWorkspace(final int dimension,
                                                  final int[] blockStarts,
                                                  final int[] blockSizes) {
        final BlockDiagonalExpSolver.BlockStructure structure =
                new BlockDiagonalExpSolver.BlockStructure(dimension, blockStarts, blockSizes);
        this.frechetHelper = new BlockDiagonalFrechetHelper(structure);
        this.lyapunovSolver = new BlockDiagonalLyapunovSolver(dimension, blockStarts, blockSizes);
        this.lyapunovAdjointHelper = new BlockDiagonalLyapunovAdjointHelper(dimension, lyapunovSolver);
        this.qMatrix = new DenseMatrix64F(dimension, dimension);
        this.qDBasis = new DenseMatrix64F(dimension, dimension);
        this.stationaryCovDBasis = new DenseMatrix64F(dimension, dimension);
        this.transitionCovDBasis = new DenseMatrix64F(dimension, dimension);
        this.transitionCovariance = new DenseMatrix64F(dimension, dimension);
        this.transitionCovarianceArrayScratch = new double[dimension * dimension];
        this.transitionOffsetScratch = new double[dimension];
        this.scaledNegativeBlockDScratch = new double[dimension][dimension];
        this.denseAdjointScratch = new double[dimension][dimension];
        this.upstreamF = new DenseMatrix64F(dimension, dimension);
        this.upstreamFD = new DenseMatrix64F(dimension, dimension);
        this.gradD = new DenseMatrix64F(dimension, dimension);
        this.gradR = new DenseMatrix64F(dimension, dimension);
        this.gV = new DenseMatrix64F(dimension, dimension);
        this.hDBasis = new DenseMatrix64F(dimension, dimension);
        this.gS = new DenseMatrix64F(dimension, dimension);
        this.yAdjoint = new DenseMatrix64F(dimension, dimension);
        this.gECov = new DenseMatrix64F(dimension, dimension);
        this.temp1 = new DenseMatrix64F(dimension, dimension);
        this.temp2 = new DenseMatrix64F(dimension, dimension);
        this.temp3 = new DenseMatrix64F(dimension, dimension);
        this.choleskyScratch = new double[dimension * dimension];
        this.lowerInverseScratch = new double[dimension * dimension];
        this.tempVector1 = new double[dimension];
        this.tempVector2 = new double[dimension];
    }
}
