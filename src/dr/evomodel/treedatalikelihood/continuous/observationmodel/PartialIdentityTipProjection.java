package dr.evomodel.treedatalikelihood.continuous.observationmodel;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTransitionMomentProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;

public final class PartialIdentityTipProjection {

    private static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);

    private final int dimension;
    private final TipObservationPartition partition;
    private final double[] shiftedObservation;
    private final double[] precisionTimesShifted;
    private final double[] varianceFlat;
    private final double[] precisionFlat;
    private final double[] projectedPrecision;
    private final double[] cholesky;
    private final double[] transitionMatrixFlat;
    private final double[] covarianceFlat;
    private final double[] transitionOffset;

    public PartialIdentityTipProjection(final int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        this.dimension = dimension;
        this.partition = new TipObservationPartition(dimension);
        this.shiftedObservation = new double[dimension];
        this.precisionTimesShifted = new double[dimension];
        this.varianceFlat = new double[dimension * dimension];
        this.precisionFlat = new double[dimension * dimension];
        this.projectedPrecision = new double[dimension * dimension];
        this.cholesky = new double[dimension * dimension];
        this.transitionMatrixFlat = new double[dimension * dimension];
        this.covarianceFlat = new double[dimension * dimension];
        this.transitionOffset = new double[dimension];
    }

    public void projectObservedChildToParent(final IdentityCanonicalTipObservationModel tipObservation,
                                             final CanonicalTransitionMomentProvider transitionMomentProvider,
                                             final double branchLength,
                                             final CanonicalGaussianState out) {
        final int observedCount = partition.update(tipObservation);
        if (observedCount == 0) {
            CanonicalGaussianMessageOps.clearState(out);
            return;
        }

        transitionMomentProvider.fillTransitionMomentsFlat(
                branchLength,
                transitionMatrixFlat,
                transitionOffset,
                covarianceFlat);

        for (int observed = 0; observed < observedCount; ++observed) {
            final int observedTrait = partition.observedIndex(observed);
            shiftedObservation[observed] =
                    tipObservation.valueAt(observedTrait) - transitionOffset[observedTrait];
            final int rowOffset = observed * observedCount;
            for (int otherObserved = 0; otherObserved < observedCount; ++otherObserved) {
                varianceFlat[rowOffset + otherObserved] =
                        covarianceFlat[observedTrait * dimension + partition.observedIndex(otherObserved)];
            }
        }

        final double logDetVariance = MatrixOps.invertSPDCompact(
                varianceFlat,
                precisionFlat,
                observedCount,
                transitionOffset,
                cholesky);

        for (int row = 0; row < observedCount; ++row) {
            double sum = 0.0;
            final int rowOffset = row * observedCount;
            for (int k = 0; k < observedCount; ++k) {
                sum += precisionFlat[rowOffset + k] * shiftedObservation[k];
            }
            precisionTimesShifted[row] = sum;
        }

        for (int row = 0; row < observedCount; ++row) {
            final int precisionRowOffset = row * observedCount;
            final int projectedRowOffset = row * dimension;
            for (int col = 0; col < dimension; ++col) {
                double sum = 0.0;
                for (int k = 0; k < observedCount; ++k) {
                    sum += precisionFlat[precisionRowOffset + k]
                            * transitionMatrixFlat[partition.observedIndex(k) * dimension + col];
                }
                projectedPrecision[projectedRowOffset + col] = sum;
            }
        }

        for (int i = 0; i < dimension; ++i) {
            double information = 0.0;
            for (int observed = 0; observed < observedCount; ++observed) {
                information += transitionMatrixFlat[partition.observedIndex(observed) * dimension + i]
                        * precisionTimesShifted[observed];
            }
            out.information[i] = information;

            for (int j = 0; j < dimension; ++j) {
                double precision = 0.0;
                for (int observed = 0; observed < observedCount; ++observed) {
                    precision += transitionMatrixFlat[partition.observedIndex(observed) * dimension + i]
                            * projectedPrecision[observed * dimension + j];
                }
                out.precision[i * dimension + j] = precision;
            }
        }
        MatrixOps.symmetrize(out.precision, dimension);

        double quadratic = 0.0;
        for (int observed = 0; observed < observedCount; ++observed) {
            quadratic += shiftedObservation[observed] * precisionTimesShifted[observed];
        }
        out.logNormalizer = 0.5 * (observedCount * LOG_TWO_PI + logDetVariance + quadratic);
    }
}
