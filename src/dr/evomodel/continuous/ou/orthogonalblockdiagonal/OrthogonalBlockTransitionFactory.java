package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovSolver;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalTransitionAdjointUtils;
import dr.inference.model.MatrixParameterInterface;
import org.ejml.data.DenseMatrix64F;

/**
 * Forward-path transition and canonical-message builder for orthogonal block OU.
 */
final class OrthogonalBlockTransitionFactory {

    private final BlockDiagonalLyapunovSolver lyapunovSolver;
    private final DenseMatrix64F qMatrix;
    private final DenseMatrix64F qDBasis;
    private final DenseMatrix64F stationaryCovDBasis;
    private final DenseMatrix64F transitionCovDBasis;
    private final DenseMatrix64F transitionCovariance;
    private final DenseMatrix64F temp;
    private final double[] transitionMatrixScratch;
    private final double[] transitionCovarianceScratch;
    private final double[] precisionScratch;
    private final double[] transitionOffsetScratch;
    private final double[] choleskyScratch;
    private final double[] lowerInverseScratch;
    private final CanonicalTransitionAdjointUtils.Workspace canonicalAdjointWorkspace;

    OrthogonalBlockTransitionFactory(final int dimension,
                                     final BlockDiagonalLyapunovSolver lyapunovSolver) {
        this.lyapunovSolver = lyapunovSolver;
        this.qMatrix = new DenseMatrix64F(dimension, dimension);
        this.qDBasis = new DenseMatrix64F(dimension, dimension);
        this.stationaryCovDBasis = new DenseMatrix64F(dimension, dimension);
        this.transitionCovDBasis = new DenseMatrix64F(dimension, dimension);
        this.transitionCovariance = new DenseMatrix64F(dimension, dimension);
        this.temp = new DenseMatrix64F(dimension, dimension);
        this.transitionMatrixScratch = new double[dimension * dimension];
        this.transitionCovarianceScratch = new double[dimension * dimension];
        this.precisionScratch = new double[dimension * dimension];
        this.transitionOffsetScratch = new double[dimension];
        this.choleskyScratch = new double[dimension * dimension];
        this.lowerInverseScratch = new double[dimension * dimension];
        this.canonicalAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(dimension);
    }

    void fillTransitionCovariance(final MatrixParameterInterface diffusionMatrix,
                                  final OrthogonalBlockBasisCache basis,
                                  final double[][] out) {
        fillTransitionCovarianceMatrix(diffusionMatrix, basis);
        copyDenseMatrixToArray(transitionCovariance, out);
    }

    void fillCanonicalTransition(final MatrixParameterInterface diffusionMatrix,
                                 final double[] stationaryMean,
                                 final OrthogonalBlockBasisCache basis,
                                 final CanonicalGaussianTransition out) {
        fillTransitionCovarianceMatrix(diffusionMatrix, basis);
        OrthogonalBlockCanonicalTransitionAssembler.fillCanonicalTransition(
                basis.transitionMatrix,
                transitionCovariance,
                stationaryMean,
                transitionOffsetScratch,
                transitionCovarianceScratch,
                choleskyScratch,
                lowerInverseScratch,
                out);
    }

    void fillCanonicalLocalAdjoints(final MatrixParameterInterface diffusionMatrix,
                                    final double[] stationaryMean,
                                    final OrthogonalBlockBasisCache basis,
                                    final CanonicalBranchMessageContribution contribution,
                                    final CanonicalLocalTransitionAdjoints out) {
        fillTransitionCovarianceMatrix(diffusionMatrix, basis);
        copyDenseMatrixToFlat(basis.transitionMatrix, transitionMatrixScratch);
        OrthogonalBlockCanonicalTransitionAssembler.fillTransitionOffset(
                basis.transitionMatrix, stationaryMean, transitionOffsetScratch);
        OrthogonalBlockPositiveDefiniteInverter.copyAndInvertFlat(
                transitionCovariance,
                transitionCovarianceScratch,
                precisionScratch,
                choleskyScratch,
                lowerInverseScratch);
        CanonicalTransitionAdjointUtils.fillFromMoments(
                precisionScratch,
                transitionCovarianceScratch,
                transitionMatrixScratch,
                transitionOffsetScratch,
                contribution,
                canonicalAdjointWorkspace,
                out);
    }

    private void fillTransitionCovarianceMatrix(final MatrixParameterInterface diffusionMatrix,
                                                final OrthogonalBlockBasisCache basis) {
        OrthogonalBlockTransitionCovarianceSolver.fillTransitionCovariance(
                diffusionMatrix,
                basis.rMatrix,
                basis.rtMatrix,
                basis.expD,
                basis.blockDParams,
                lyapunovSolver,
                qMatrix,
                qDBasis,
                stationaryCovDBasis,
                transitionCovDBasis,
                temp,
                transitionCovariance,
                false);
    }

    private static void copyDenseMatrixToArray(final DenseMatrix64F source, final double[][] out) {
        final int dimension = source.numRows;
        final double[] data = source.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i][j] = data[i * dimension + j];
            }
        }
    }

    private static void copyDenseMatrixToFlat(final DenseMatrix64F source, final double[] out) {
        System.arraycopy(source.data, 0, out, 0, source.numRows * source.numCols);
    }
}
