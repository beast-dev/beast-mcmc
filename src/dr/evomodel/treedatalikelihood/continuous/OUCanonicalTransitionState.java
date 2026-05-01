package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumerics;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumericsOptions;
import dr.evomodel.treedatalikelihood.preorder.MatrixSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import org.ejml.data.DenseMatrix64F;

/**
 * Owns the branch-local OU transition buffers and SPD canonicalization.
 */
final class OUCanonicalTransitionState {

    private final OUGaussianBranchTransitionProvider branchTransitionProvider;
    private final int dimension;
    private final CanonicalGaussianTransition transition;
    private final DenseMatrix64F actualizationMatrix;
    private final DenseMatrix64F displacementVector;
    private final DenseMatrix64F branchPrecisionMatrix;
    private final DenseMatrix64F branchCovarianceMatrix;
    private final OUCanonicalParentAboveResolver parentAboveResolver;
    private final DenseMatrix64F scratchMatrix;
    private final DenseMatrix64F secondaryScratchMatrix;
    private final double[] spdMatrixScratch;
    private final double[] spdInverseScratch;
    private final double[] spdCholeskyScratch;
    private final double[] spdLowerInverseScratch;
    private final double[] transitionMatrixScratch;
    private final double[] transitionCovarianceScratch;
    private final double[] transitionOffsetScratch;
    private final double[] transitionPrecisionScratch;
    private final double[] transitionSymmetricScratch;
    private final double[] transitionCholeskyScratch;
    private final double[] transitionLowerInverseScratch;
    private final double[] transitionTransposeScratch;
    private final double[] matrixScratch;

    OUCanonicalTransitionState(final OUGaussianBranchTransitionProvider branchTransitionProvider,
                               final CanonicalGaussianTransition transition,
                               final CanonicalNumericsOptions numericsOptions) {
        this.branchTransitionProvider = branchTransitionProvider;
        this.dimension = branchTransitionProvider.getStateDimension();
        this.transition = transition;
        this.actualizationMatrix = new DenseMatrix64F(dimension, dimension);
        this.displacementVector = new DenseMatrix64F(dimension, 1);
        this.branchPrecisionMatrix = new DenseMatrix64F(dimension, dimension);
        this.branchCovarianceMatrix = new DenseMatrix64F(dimension, dimension);
        this.parentAboveResolver = new OUCanonicalParentAboveResolver(
                dimension,
                actualizationMatrix,
                displacementVector,
                branchPrecisionMatrix,
                branchCovarianceMatrix,
                numericsOptions);
        this.scratchMatrix = new DenseMatrix64F(dimension, dimension);
        this.secondaryScratchMatrix = new DenseMatrix64F(dimension, dimension);
        final int maxReducedDimension = 2 * dimension;
        this.spdMatrixScratch = new double[maxReducedDimension * maxReducedDimension];
        this.spdInverseScratch = new double[maxReducedDimension * maxReducedDimension];
        this.spdCholeskyScratch = new double[maxReducedDimension * maxReducedDimension];
        this.spdLowerInverseScratch = new double[maxReducedDimension * maxReducedDimension];
        this.transitionMatrixScratch = new double[dimension * dimension];
        this.transitionCovarianceScratch = new double[dimension * dimension];
        this.transitionOffsetScratch = new double[dimension];
        this.transitionPrecisionScratch = new double[dimension * dimension];
        this.transitionSymmetricScratch = new double[dimension * dimension];
        this.transitionCholeskyScratch = new double[dimension * dimension];
        this.transitionLowerInverseScratch = new double[dimension * dimension];
        this.transitionTransposeScratch = new double[dimension * dimension];
        this.matrixScratch = new double[dimension * dimension];
    }

    void fillFromStatistics(final MatrixSufficientStatistics branchStatistics) {
        final DenseMatrix64F actualization = branchStatistics.getRawActualization();
        final DenseMatrix64F displacement = branchStatistics.getRawDisplacement();
        final DenseMatrix64F precision = branchStatistics.getRawPrecision();

        if (!CanonicalNumerics.isFinite(actualization)
                || !CanonicalNumerics.isFinite(displacement)
                || !CanonicalNumerics.isFinite(precision)) {
            throw new IllegalStateException(
                    "Non-finite branch transition statistics in fillCanonicalTransition"
                            + "; actualizationFinite=" + CanonicalNumerics.isFinite(actualization)
                            + "; displacementFinite=" + CanonicalNumerics.isFinite(displacement)
                            + "; precisionFinite=" + CanonicalNumerics.isFinite(precision));
        }

        actualizationMatrix.set(actualization);
        displacementVector.set(displacement);
        branchPrecisionMatrix.set(precision);
        canonicalizeBranchPrecisionCovariancePair("fillCanonicalTransition");
        OUCanonicalTransitionBuilder.fillFromPrecisionMoments(
                actualizationMatrix, displacementVector, branchPrecisionMatrix, transition);
    }

    void fillFromKernel(final double branchLength, final double[] optimum) {
        final double[] transitionMatrix = transitionMatrixScratch;
        final double[] transitionCovariance = transitionCovarianceScratch;
        final double[] transitionOffset = transitionOffsetScratch;
        fillTransitionMomentsFromKernel(branchLength, optimum, transitionMatrix, transitionOffset, transitionCovariance);

        for (int i = 0; i < dimension; ++i) {
            displacementVector.unsafe_set(i, 0, transitionOffset[i]);
            final int iOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                actualizationMatrix.unsafe_set(i, j, transitionMatrix[iOffset + j]);
                branchCovarianceMatrix.unsafe_set(i, j, transitionCovariance[iOffset + j]);
            }
        }
        canonicalizeBranchCovariancePrecisionPair(
                "fillCanonicalTransitionFromKernel(branchLength=" + branchLength + ")");

        if (branchLength <= 0.0) {
            OUCanonicalTransitionBuilder.fillFromMoments(
                    transitionMatrix,
                    transitionOffset,
                    transitionCovariance,
                    transition,
                    transitionPrecisionScratch,
                    transitionSymmetricScratch,
                    transitionCholeskyScratch,
                    transitionLowerInverseScratch,
                    transitionTransposeScratch,
                    matrixScratch);
        } else {
            OUCanonicalTransitionBuilder.fillFromPrecisionMoments(
                    actualizationMatrix, displacementVector, branchPrecisionMatrix, transition);
        }
    }

    void fillTransitionMomentsFromKernel(final double branchLength,
                                         final double[] optimum,
                                         final double[] transitionMatrixOut,
                                         final double[] transitionOffsetOut,
                                         final double[] transitionCovarianceOut) {
        branchTransitionProvider.fillBranchTransitionMatrixFlat(branchLength, transitionMatrixOut);
        branchTransitionProvider.fillBranchTransitionCovarianceFlat(branchLength, transitionCovarianceOut);
        fillTransitionOffset(transitionMatrixOut, optimum, transitionOffsetOut);
    }

    DenseMatrix64F recoverOrUseParentAbovePrecision(final NormalSufficientStatistics aboveChild,
                                                    final NormalSufficientStatistics aboveParent) {
        return parentAboveResolver.recoverOrUseParentAbovePrecision(aboveChild, aboveParent);
    }

    DenseMatrix64F requireParentAbovePrecision(final NormalSufficientStatistics aboveParent) {
        return parentAboveResolver.requireParentAbovePrecision(aboveParent);
    }

    DenseMatrix64F recoverParentAbovePrecision(final NormalSufficientStatistics above) {
        return parentAboveResolver.recoverParentAbovePrecision(above);
    }

    DenseMatrix64F getAboveParentMeanVector() {
        return parentAboveResolver.getAboveParentMeanVector();
    }

    private void canonicalizeBranchPrecisionCovariancePair(final String context) {
        CanonicalNumerics.symmetrizeInPlace(branchPrecisionMatrix);
        safeInvertSymmetricPositiveDefinite(
                branchPrecisionMatrix,
                branchCovarianceMatrix,
                context + ":branchPrecision->branchCovariance");
        safeInvertSymmetricPositiveDefinite(
                branchCovarianceMatrix,
                branchPrecisionMatrix,
                context + ":branchCovariance->branchPrecision");
        CanonicalNumerics.symmetrizeInPlace(branchPrecisionMatrix);
        CanonicalNumerics.symmetrizeInPlace(branchCovarianceMatrix);
    }

    private void canonicalizeBranchCovariancePrecisionPair(final String context) {
        CanonicalNumerics.symmetrizeInPlace(branchCovarianceMatrix);
        safeInvertSymmetricPositiveDefinite(
                branchCovarianceMatrix,
                branchPrecisionMatrix,
                context + ":branchCovariance->branchPrecision");
        safeInvertSymmetricPositiveDefinite(
                branchPrecisionMatrix,
                branchCovarianceMatrix,
                context + ":branchPrecision->branchCovariance");
        CanonicalNumerics.symmetrizeInPlace(branchPrecisionMatrix);
        CanonicalNumerics.symmetrizeInPlace(branchCovarianceMatrix);
    }

    private void safeInvertSymmetricPositiveDefinite(final DenseMatrix64F source,
                                                     final DenseMatrix64F inverseOut,
                                                     final String context) {
        CanonicalNumerics.copyAndInvertPositiveDefiniteFlat(
                source,
                spdMatrixScratch,
                spdInverseScratch,
                spdCholeskyScratch,
                spdLowerInverseScratch,
                CanonicalNumericsOptions.OU_TREE);
        System.arraycopy(spdInverseScratch, 0, inverseOut.data, 0, dimension * dimension);
    }

    private static void fillTransitionOffset(final double[] transitionMatrix,
                                             final double[] optimum,
                                             final double[] out) {
        final int d = out.length;
        for (int i = 0; i < d; ++i) {
            double sum = optimum[i];
            final int iOffset = i * d;
            for (int j = 0; j < d; ++j) {
                sum -= transitionMatrix[iOffset + j] * optimum[j];
            }
            out[i] = sum;
        }
    }
}
