package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumerics;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.SufficientStatisticsTipObservationPattern;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import org.ejml.data.DenseMatrix64F;

/**
 * Evaluates frozen branch-local canonical log factors for branch-length scores.
 */
final class OUCanonicalFrozenLocalFactorEvaluator {

    private final int dimension;
    private final CanonicalGaussianTransition transition;
    private final SufficientStatisticsTipObservationPattern observationPattern;
    private final double[] aboveInformation;
    private final double[] belowInformation;
    private final double[][] reducedPrecisionScratch;
    private final double[][] reducedCovarianceScratch;
    private final double[] reducedInformationScratch;
    private final double[] reducedMeanScratch;
    private final double[][] reducedCholeskyScratch;
    private final double[][] reducedLowerInverseScratch;

    OUCanonicalFrozenLocalFactorEvaluator(final int dimension,
                                          final CanonicalGaussianTransition transition,
                                          final SufficientStatisticsTipObservationPattern observationPattern) {
        this.dimension = dimension;
        this.transition = transition;
        this.observationPattern = observationPattern;
        final int maxReducedDimension = 2 * dimension;
        this.aboveInformation = new double[dimension];
        this.belowInformation = new double[dimension];
        this.reducedPrecisionScratch = new double[maxReducedDimension][maxReducedDimension];
        this.reducedCovarianceScratch = new double[maxReducedDimension][maxReducedDimension];
        this.reducedInformationScratch = new double[maxReducedDimension];
        this.reducedMeanScratch = new double[maxReducedDimension];
        this.reducedCholeskyScratch = new double[maxReducedDimension][maxReducedDimension];
        this.reducedLowerInverseScratch = new double[maxReducedDimension][maxReducedDimension];
    }

    double evaluate(final DenseMatrix64F parentAbovePrecision,
                    final DenseMatrix64F parentAboveMean,
                    final NormalSufficientStatistics below) {
        fillInformation(parentAbovePrecision, parentAboveMean, aboveInformation);

        final int observedCount = observationPattern.classify(below);
        if (observedCount == dimension) {
            return evaluateFullyObserved(parentAbovePrecision, below.getRawMean());
        }

        if (observedCount > 0) {
            return evaluatePartiallyObserved(parentAbovePrecision, below.getRawMean(), observedCount);
        }

        fillInformation(below.getRawPrecision(), below.getRawMean(), belowInformation);
        final int doubled = 2 * dimension;
        final double[][] precision = reducedPrecisionScratch;
        final double[] information = reducedInformationScratch;
        for (int i = 0; i < dimension; ++i) {
            information[i] = transition.informationX[i] + aboveInformation[i];
            information[dimension + i] = transition.informationY[i] + belowInformation[i];
            final int iOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                precision[i][j] = transition.precisionXX[iOffset + j] + parentAbovePrecision.unsafe_get(i, j);
                precision[i][dimension + j] = transition.precisionXY[iOffset + j];
                precision[dimension + i][j] = transition.precisionYX[iOffset + j];
                precision[dimension + i][dimension + j] =
                        transition.precisionYY[iOffset + j] + below.getRawPrecision().unsafe_get(i, j);
            }
        }
        return normalizedLogNormalizer(precision, information, doubled) - transition.logNormalizer;
    }

    private double evaluateFullyObserved(final DenseMatrix64F parentAbovePrecision,
                                         final DenseMatrix64F observed) {
        final double[][] precision = reducedPrecisionScratch;
        final double[] information = reducedInformationScratch;
        for (int i = 0; i < dimension; ++i) {
            double info = transition.informationX[i] + aboveInformation[i];
            final int iOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                precision[i][j] = transition.precisionXX[iOffset + j] + parentAbovePrecision.unsafe_get(i, j);
                info -= transition.precisionXY[iOffset + j] * observed.unsafe_get(j, 0);
            }
            information[i] = info;
        }
        return normalizedLogNormalizer(precision, information, dimension)
                - transition.logNormalizer
                - 0.5 * quadraticForm(transition.precisionYY, observed)
                + dot(transition.informationY, observed);
    }

    private double evaluatePartiallyObserved(final DenseMatrix64F parentAbovePrecision,
                                             final DenseMatrix64F observedChild,
                                             final int observedCount) {
        final int missingCount = observationPattern.collectPartition(observedCount);
        final int reducedDimension = dimension + missingCount;
        final double[][] reducedPrecision = reducedPrecisionScratch;
        final double[] reducedInformation = reducedInformationScratch;
        fillReducedCanonicalSystem(parentAbovePrecision, observedChild, missingCount, reducedPrecision, reducedInformation);

        double observedQuadratic = 0.0;
        double observedLinear = 0.0;
        for (int oi = 0; oi < observedCount; ++oi) {
            final int i = observationPattern.observedIndex(oi);
            final double yi = observedChild.unsafe_get(i, 0);
            observedLinear += transition.informationY[i] * yi;
            final int iOffset = i * dimension;
            for (int oj = 0; oj < observedCount; ++oj) {
                final int j = observationPattern.observedIndex(oj);
                observedQuadratic += yi * transition.precisionYY[iOffset + j] * observedChild.unsafe_get(j, 0);
            }
        }

        return normalizedLogNormalizer(reducedPrecision, reducedInformation, reducedDimension)
                - transition.logNormalizer
                - 0.5 * observedQuadratic
                + observedLinear;
    }

    void fillReducedCanonicalSystem(final DenseMatrix64F abovePrecision,
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

    private double normalizedLogNormalizer(final double[][] precision,
                                           final double[] information,
                                           final int dimensionUsed) {
        final double[][] inverse = reducedCovarianceScratch;
        final double[] mean = reducedMeanScratch;
        final double logDet = invertSymmetricPositiveDefinite(precision, dimensionUsed, inverse);
        for (int i = 0; i < dimensionUsed; ++i) {
            double sum = 0.0;
            for (int j = 0; j < dimensionUsed; ++j) {
                sum += inverse[i][j] * information[j];
            }
            mean[i] = sum;
        }
        return 0.5 * (dimensionUsed * Math.log(2.0 * Math.PI) - logDet + dot(information, mean, dimensionUsed));
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

    private static double quadraticForm(final double[] matrix, final DenseMatrix64F vector) {
        double sum = 0.0;
        final int dimensionUsed = vector.numRows;
        for (int i = 0; i < dimensionUsed; ++i) {
            final int iOffset = i * dimensionUsed;
            for (int j = 0; j < dimensionUsed; ++j) {
                sum += vector.unsafe_get(i, 0) * matrix[iOffset + j] * vector.unsafe_get(j, 0);
            }
        }
        return sum;
    }

    private static double dot(final double[] left, final DenseMatrix64F right) {
        double sum = 0.0;
        for (int i = 0; i < left.length; ++i) {
            sum += left[i] * right.unsafe_get(i, 0);
        }
        return sum;
    }

    private static double dot(final double[] left, final double[] right, final int dimensionUsed) {
        double sum = 0.0;
        for (int i = 0; i < dimensionUsed; ++i) {
            sum += left[i] * right[i];
        }
        return sum;
    }
}
