package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumerics;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.SufficientStatisticsTipObservationPattern;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import org.ejml.data.DenseMatrix64F;

final class OUCanonicalBranchContributionAssembler {

    private final int dimension;
    private final CanonicalGradientFallbackPolicy fallbackPolicy;
    private final CanonicalGaussianTransition transition;
    private final OUCanonicalTransitionState transitionState;
    private final OUCanonicalParentAboveMessages parentAboveMessages;
    private final OUCanonicalFrozenLocalFactorEvaluator frozenLocalFactorEvaluator;
    private final CanonicalBranchMessageContribution contribution;
    private final CanonicalGaussianState pairPosterior;
    private final SufficientStatisticsTipObservationPattern observationPattern;

    private final double[] aboveInformation;
    private final double[] belowInformation;
    private final double[][] reducedPrecisionScratch;
    private final double[][] reducedCovarianceScratch;
    private final double[] reducedInformationScratch;
    private final double[] reducedMeanScratch;
    private final double[][] reducedCholeskyScratch;
    private final double[][] reducedLowerInverseScratch;

    OUCanonicalBranchContributionAssembler(final int dimension,
                                           final CanonicalGradientFallbackPolicy fallbackPolicy,
                                           final CanonicalGaussianTransition transition,
                                           final OUCanonicalTransitionState transitionState,
                                           final OUCanonicalParentAboveMessages parentAboveMessages) {
        this.dimension = dimension;
        this.fallbackPolicy = fallbackPolicy;
        this.transition = transition;
        this.transitionState = transitionState;
        this.parentAboveMessages = parentAboveMessages;
        this.contribution = new CanonicalBranchMessageContribution(dimension);
        this.pairPosterior = new CanonicalGaussianState(2 * dimension);
        this.observationPattern = new SufficientStatisticsTipObservationPattern(dimension);
        this.frozenLocalFactorEvaluator = new OUCanonicalFrozenLocalFactorEvaluator(
                dimension,
                transition,
                observationPattern);
        this.aboveInformation = new double[dimension];
        this.belowInformation = new double[dimension];
        final int maxReducedDimension = 2 * dimension;
        this.reducedPrecisionScratch = new double[maxReducedDimension][maxReducedDimension];
        this.reducedCovarianceScratch = new double[maxReducedDimension][maxReducedDimension];
        this.reducedInformationScratch = new double[maxReducedDimension];
        this.reducedMeanScratch = new double[maxReducedDimension];
        this.reducedCholeskyScratch = new double[maxReducedDimension][maxReducedDimension];
        this.reducedLowerInverseScratch = new double[maxReducedDimension][maxReducedDimension];
    }

    CanonicalBranchMessageContribution prepareFromCurrentTransition(final BranchSufficientStatistics statistics) {
        final int observedCount = fallbackPolicy.isExactTipShortcutDisabled()
                ? 0
                : classifyExactObservationPattern(statistics.getBelow());
        if (observedCount > 0) {
            fillContributionForPartiallyObservedTip(
                    statistics.getAboveParent(),
                    statistics.getBelow(),
                    observedCount);
        } else {
            fillPairPosterior(
                    statistics.getAboveParent(),
                    statistics.getBelow());
            fillContributionFromPairPosterior();
        }
        return contribution;
    }

    double evaluateFrozenLocalLogFactor(final double branchLength,
                                        final double[] optimum,
                                        final NormalSufficientStatistics aboveParent,
                                        final NormalSufficientStatistics below) {
        transitionState.fillFromKernel(branchLength, optimum);
        final DenseMatrix64F parentAbovePrecision = parentAboveMessages.require(aboveParent);
        final DenseMatrix64F parentAboveMean = transitionState.getAboveParentMeanVector();
        return frozenLocalFactorEvaluator.evaluate(parentAbovePrecision, parentAboveMean, below);
    }

    private void fillPairPosterior(final NormalSufficientStatistics aboveParent,
                                   final NormalSufficientStatistics below) {
        final DenseMatrix64F abovePrecision = parentAboveMessages.require(aboveParent);
        final DenseMatrix64F aboveMean = transitionState.getAboveParentMeanVector();
        final DenseMatrix64F belowPrecision = below.getRawPrecision();
        final DenseMatrix64F belowMean = below.getRawMean();

        fillInformation(abovePrecision, aboveMean, aboveInformation);
        fillInformation(belowPrecision, belowMean, belowInformation);

        final int pairDimension = pairPosterior.getDimension();
        for (int i = 0; i < dimension; ++i) {
            pairPosterior.information[i] = transition.informationX[i] + aboveInformation[i];
            pairPosterior.information[dimension + i] = transition.informationY[i] + belowInformation[i];
            final int pairRowOffset = i * pairDimension;
            final int transitionRowOffset = i * dimension;
            final int lowerRowOffset = (dimension + i) * pairDimension;
            for (int j = 0; j < dimension; ++j) {
                pairPosterior.precision[pairRowOffset + j] =
                        transition.precisionXX[transitionRowOffset + j] + abovePrecision.unsafe_get(i, j);
                pairPosterior.precision[pairRowOffset + dimension + j] =
                        transition.precisionXY[transitionRowOffset + j];
                pairPosterior.precision[lowerRowOffset + j] =
                        transition.precisionYX[transitionRowOffset + j];
                pairPosterior.precision[lowerRowOffset + dimension + j] =
                        transition.precisionYY[transitionRowOffset + j] + belowPrecision.unsafe_get(i, j);
            }
        }
        pairPosterior.logNormalizer = 0.0;

        if (!CanonicalNumerics.isFinite(pairPosterior.precision)
                || !CanonicalNumerics.isFinite(pairPosterior.information)) {
            throw new IllegalStateException(
                    "Non-finite pair posterior before canonical->moments conversion"
                            + "; abovePrecisionFinite=" + CanonicalNumerics.isFinite(abovePrecision)
                            + "; belowPrecisionFinite=" + CanonicalNumerics.isFinite(belowPrecision)
                            + "; belowMeanFinite=" + CanonicalNumerics.isFinite(belowMean)
                            + "; transitionXXFinite=" + CanonicalNumerics.isFinite(transition.precisionXX)
                            + "; transitionXYFinite=" + CanonicalNumerics.isFinite(transition.precisionXY)
                            + "; transitionYYFinite=" + CanonicalNumerics.isFinite(transition.precisionYY)
                            + "; transitionInfoXFinite=" + CanonicalNumerics.isFinite(transition.informationX)
                            + "; transitionInfoYFinite=" + CanonicalNumerics.isFinite(transition.informationY));
        }
    }

    private void fillContributionForPartiallyObservedTip(final NormalSufficientStatistics aboveParent,
                                                         final NormalSufficientStatistics below,
                                                         final int observedCount) {
        final DenseMatrix64F abovePrecision = parentAboveMessages.require(aboveParent);
        final DenseMatrix64F aboveMean = transitionState.getAboveParentMeanVector();
        final DenseMatrix64F observedChild = below.getRawMean();

        fillInformation(abovePrecision, aboveMean, aboveInformation);
        final int missingCount = collectObservationPartition(observedCount);
        final int reducedDimension = dimension + missingCount;
        final double[][] reducedPrecision = reducedPrecisionScratch;
        final double[][] reducedCovariance = reducedCovarianceScratch;
        final double[] reducedInformation = reducedInformationScratch;
        final double[] reducedMean = reducedMeanScratch;

        fillReducedCanonicalSystem(abovePrecision, observedChild, missingCount, reducedPrecision, reducedInformation);
        invertSymmetricPositiveDefinite(reducedPrecision, reducedDimension, reducedCovariance);
        for (int i = 0; i < reducedDimension; ++i) {
            double sum = 0.0;
            for (int j = 0; j < reducedDimension; ++j) {
                sum += reducedCovariance[i][j] * reducedInformation[j];
            }
            reducedMean[i] = sum;
        }

        OUCanonicalBranchContributionBuilder.fillPartiallyObservedTipContribution(
                dimension,
                observationPattern,
                observedChild,
                reducedMean,
                reducedCovariance,
                contribution);
    }

    private void fillContributionFromPairPosterior() {
        fillMomentsFromCanonicalLocally(pairPosterior, reducedMeanScratch, reducedCovarianceScratch);

        OUCanonicalBranchContributionBuilder.fillPairPosteriorContribution(
                dimension,
                reducedMeanScratch,
                reducedCovarianceScratch,
                contribution);
    }

    private int classifyExactObservationPattern(final NormalSufficientStatistics below) {
        return observationPattern.classify(below);
    }

    private int collectObservationPartition(final int observedCount) {
        return observationPattern.collectPartition(observedCount);
    }

    private void fillReducedCanonicalSystem(final DenseMatrix64F abovePrecision,
                                            final DenseMatrix64F observedChild,
                                            final int missingCount,
                                            final double[][] reducedPrecision,
                                            final double[] reducedInformation) {
        for (int i = 0; i < dimension; ++i) {
            double information = transition.informationX[i] + aboveInformation[i];
            final int iOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                reducedPrecision[i][j] = transition.precisionXX[iOffset + j] + abovePrecision.unsafe_get(i, j);
            }
            for (int observed = 0; observed < dimension - missingCount; ++observed) {
                final int observedIndex = observationPattern.observedIndex(observed);
                information -= transition.precisionXY[iOffset + observedIndex] * observedChild.unsafe_get(observedIndex, 0);
            }
            reducedInformation[i] = information;
            for (int missing = 0; missing < missingCount; ++missing) {
                reducedPrecision[i][dimension + missing] =
                        transition.precisionXY[iOffset + observationPattern.missingIndex(missing)];
            }
        }

        for (int missing = 0; missing < missingCount; ++missing) {
            final int childIndex = observationPattern.missingIndex(missing);
            final int row = dimension + missing;
            double information = transition.informationY[childIndex];
            final int childOffset = childIndex * dimension;
            for (int observed = 0; observed < dimension - missingCount; ++observed) {
                final int observedIndex = observationPattern.observedIndex(observed);
                information -= transition.precisionYY[childOffset + observedIndex] * observedChild.unsafe_get(observedIndex, 0);
            }
            reducedInformation[row] = information;
            for (int j = 0; j < dimension; ++j) {
                reducedPrecision[row][j] = transition.precisionYX[childOffset + j];
            }
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                reducedPrecision[row][dimension + otherMissing] =
                        transition.precisionYY[childOffset + observationPattern.missingIndex(otherMissing)];
            }
        }
    }

    private void fillMomentsFromCanonicalLocally(final CanonicalGaussianState state,
                                                 final double[] meanOut,
                                                 final double[][] covarianceOut) {
        final int dimensionUsed = state.getDimension();
        invertFlatSymmetricPositiveDefinite(state.precision, dimensionUsed, covarianceOut);
        for (int i = 0; i < dimensionUsed; ++i) {
            double sum = 0.0;
            for (int j = 0; j < dimensionUsed; ++j) {
                sum += covarianceOut[i][j] * state.information[j];
            }
            meanOut[i] = sum;
        }
    }

    private double invertSymmetricPositiveDefinite(final double[][] matrix,
                                                   final int dimensionUsed,
                                                   final double[][] inverseOut) {
        return CanonicalNumerics.invertSymmetricPositiveDefinite(
                matrix,
                dimensionUsed,
                inverseOut,
                reducedCholeskyScratch,
                reducedPrecisionScratch,
                reducedLowerInverseScratch);
    }

    private double invertFlatSymmetricPositiveDefinite(final double[] matrix,
                                                       final int dimensionUsed,
                                                       final double[][] inverseOut) {
        final double[][] square = reducedPrecisionScratch;
        for (int i = 0; i < dimensionUsed; ++i) {
            System.arraycopy(matrix, i * dimensionUsed, square[i], 0, dimensionUsed);
        }
        return invertSymmetricPositiveDefinite(square, dimensionUsed, inverseOut);
    }

    private static void fillInformation(final DenseMatrix64F precision,
                                        final DenseMatrix64F mean,
                                        final double[] out) {
        final int d = out.length;
        for (int i = 0; i < d; ++i) {
            double sum = 0.0;
            for (int j = 0; j < d; ++j) {
                sum += precision.unsafe_get(i, j) * mean.unsafe_get(j, 0);
            }
            out[i] = sum;
        }
    }
}
