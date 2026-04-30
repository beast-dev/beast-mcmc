package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContribution;
import org.ejml.data.DenseMatrix64F;

/**
 * Converts branch-local posterior moments into canonical transition
 * contribution adjoints.
 */
final class OUCanonicalBranchContributionBuilder {

    private OUCanonicalBranchContributionBuilder() { }

    static void clear(final CanonicalBranchMessageContribution contribution) {
        zero(contribution.dLogL_dInformationX);
        zero(contribution.dLogL_dInformationY);
        zero(contribution.dLogL_dPrecisionXX);
        zero(contribution.dLogL_dPrecisionXY);
        zero(contribution.dLogL_dPrecisionYX);
        zero(contribution.dLogL_dPrecisionYY);
        contribution.dLogL_dLogNormalizer = 0.0;
    }

    static void fillObservedTipContribution(final int dimension,
                                            final double[] parentMean,
                                            final double[][] parentCovariance,
                                            final DenseMatrix64F observedChild,
                                            final CanonicalBranchMessageContribution contribution) {
        for (int i = 0; i < dimension; ++i) {
            final double xi = parentMean[i];
            final double yi = observedChild.unsafe_get(i, 0);
            contribution.dLogL_dInformationX[i] = xi;
            contribution.dLogL_dInformationY[i] = yi;
            for (int j = 0; j < dimension; ++j) {
                final double xj = parentMean[j];
                final double yj = observedChild.unsafe_get(j, 0);
                final double exx = parentCovariance[i][j] + xi * xj;
                final int ij = i * dimension + j;
                contribution.dLogL_dPrecisionXX[ij] = -0.5 * exx;
                contribution.dLogL_dPrecisionXY[ij] = -0.5 * (xi * yj);
                contribution.dLogL_dPrecisionYX[ij] = -0.5 * (yi * xj);
                contribution.dLogL_dPrecisionYY[ij] = -0.5 * (yi * yj);
            }
        }
        contribution.dLogL_dLogNormalizer = -1.0;
    }

    static void fillPartiallyObservedTipContribution(final int dimension,
                                                     final OUCanonicalObservationPattern observationPattern,
                                                     final DenseMatrix64F observedChild,
                                                     final double[] reducedMean,
                                                     final double[][] reducedCovariance,
                                                     final CanonicalBranchMessageContribution contribution) {
        clear(contribution);
        for (int i = 0; i < dimension; ++i) {
            contribution.dLogL_dInformationX[i] = reducedMean[i];
            contribution.dLogL_dInformationY[i] = observationPattern.isObserved(i)
                    ? observedChild.unsafe_get(i, 0)
                    : reducedMean[observationPattern.reducedIndexByTrait(i)];
        }

        for (int i = 0; i < dimension; ++i) {
            final double xi = reducedMean[i];
            final int reducedI = observationPattern.reducedIndexByTrait(i);
            for (int j = 0; j < dimension; ++j) {
                final double xj = reducedMean[j];
                final int reducedJ = observationPattern.reducedIndexByTrait(j);
                final double yi = contribution.dLogL_dInformationY[i];
                final double yj = contribution.dLogL_dInformationY[j];

                final double exx = reducedCovariance[i][j] + xi * xj;
                final double exy = observationPattern.isObserved(j)
                        ? xi * yj
                        : reducedCovariance[i][reducedJ] + xi * yj;
                final double eyx = observationPattern.isObserved(i)
                        ? yi * xj
                        : reducedCovariance[reducedI][j] + yi * xj;
                final double eyy = (observationPattern.isObserved(i) || observationPattern.isObserved(j))
                        ? yi * yj
                        : reducedCovariance[reducedI][reducedJ] + yi * yj;

                final int ij = i * dimension + j;
                contribution.dLogL_dPrecisionXX[ij] = -0.5 * exx;
                contribution.dLogL_dPrecisionXY[ij] = -0.5 * exy;
                contribution.dLogL_dPrecisionYX[ij] = -0.5 * eyx;
                contribution.dLogL_dPrecisionYY[ij] = -0.5 * eyy;
            }
        }
        contribution.dLogL_dLogNormalizer = -1.0;
    }

    static void fillPairPosteriorContribution(final int dimension,
                                              final double[] pairMean,
                                              final double[][] pairCovariance,
                                              final CanonicalBranchMessageContribution contribution) {
        for (int i = 0; i < dimension; ++i) {
            final double currentMeanI = pairMean[i];
            final double nextMeanI = pairMean[dimension + i];
            contribution.dLogL_dInformationX[i] = currentMeanI;
            contribution.dLogL_dInformationY[i] = nextMeanI;
            for (int j = 0; j < dimension; ++j) {
                final double currentMeanJ = pairMean[j];
                final double nextMeanJ = pairMean[dimension + j];
                final double currentSecondMoment =
                        pairCovariance[i][j] + currentMeanI * currentMeanJ;
                final double crossSecondMoment =
                        pairCovariance[dimension + i][j] + nextMeanI * currentMeanJ;
                final double nextSecondMoment =
                        pairCovariance[dimension + i][dimension + j] + nextMeanI * nextMeanJ;
                final int ij = i * dimension + j;
                contribution.dLogL_dPrecisionXX[ij] = -0.5 * currentSecondMoment;
                contribution.dLogL_dPrecisionXY[ij] =
                        -0.5 * (pairCovariance[i][dimension + j] + currentMeanI * nextMeanJ);
                contribution.dLogL_dPrecisionYX[ij] = -0.5 * crossSecondMoment;
                contribution.dLogL_dPrecisionYY[ij] = -0.5 * nextSecondMoment;
            }
        }
        contribution.dLogL_dLogNormalizer = -1.0;
    }

    private static void zero(final double[] vector) {
        for (int i = 0; i < vector.length; ++i) {
            vector[i] = 0.0;
        }
    }
}
