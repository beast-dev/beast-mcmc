package dr.evomodel.treedatalikelihood.continuous.observationmodel;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;

import java.util.Arrays;

public final class PartialIdentityTipContribution {

    private final int dimension;

    public PartialIdentityTipContribution(final int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        this.dimension = dimension;
    }

    public void fillContributionForPartiallyObservedTip(
            final CanonicalGaussianState aboveState,
            final CanonicalGaussianTransition transition,
            final CanonicalTipObservation tipObservation,
            final PartialIdentityTipObservationWorkspace workspace,
            final CanonicalBranchMessageContribution contribution) {
        final TipObservationPartition partition = workspace.partition;
        final int observedCount = partition.update(tipObservation);
        final int missingCount = partition.getMissingCount();
        final int reducedDimension = dimension + missingCount;

        for (int i = 0; i < dimension; ++i) {
            double info = transition.informationX[i] + aboveState.information[i];
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                workspace.reducedPrecision[i * reducedDimension + j] =
                        transition.precisionXX[iOff + j] + aboveState.precision[iOff + j];
            }
            for (int observed = 0; observed < observedCount; ++observed) {
                final int observedIndex = partition.observedIndex(observed);
                info -= transition.precisionXY[iOff + observedIndex] * tipObservation.values[observedIndex];
            }
            workspace.reducedInformation[i] = info;
            for (int missing = 0; missing < missingCount; ++missing) {
                workspace.reducedPrecision[i * reducedDimension + dimension + missing] =
                        transition.precisionXY[iOff + partition.missingIndex(missing)];
            }
        }

        for (int missing = 0; missing < missingCount; ++missing) {
            final int childIndex = partition.missingIndex(missing);
            final int row = dimension + missing;
            final int childOff = childIndex * dimension;
            double info = transition.informationY[childIndex];
            for (int observed = 0; observed < observedCount; ++observed) {
                final int observedIndex = partition.observedIndex(observed);
                info -= transition.precisionYY[childOff + observedIndex] * tipObservation.values[observedIndex];
            }
            workspace.reducedInformation[row] = info;
            for (int j = 0; j < dimension; ++j) {
                workspace.reducedPrecision[row * reducedDimension + j] =
                        transition.precisionYX[childOff + j];
            }
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                workspace.reducedPrecision[row * reducedDimension + dimension + otherMissing] =
                        transition.precisionYY[childOff + partition.missingIndex(otherMissing)];
            }
        }

        MatrixOps.safeInvertPrecision(
                workspace.reducedPrecision,
                workspace.reducedCovariance,
                reducedDimension);

        for (int i = 0; i < reducedDimension; ++i) {
            double sum = 0.0;
            for (int j = 0; j < reducedDimension; ++j) {
                sum += workspace.reducedCovariance[i * reducedDimension + j]
                        * workspace.reducedInformation[j];
            }
            workspace.reducedMean[i] = Double.isNaN(sum) ? 0.0 : sum;
        }

        clearContribution(contribution);
        for (int i = 0; i < dimension; ++i) {
            contribution.dLogL_dInformationX[i] = workspace.reducedMean[i];
            contribution.dLogL_dInformationY[i] = tipObservation.observed[i]
                    ? tipObservation.values[i]
                    : workspace.reducedMean[partition.reducedIndexByTrait(i)];
        }

        for (int i = 0; i < dimension; ++i) {
            final double xi = workspace.reducedMean[i];
            final int reducedI = partition.reducedIndexByTrait(i);
            for (int j = 0; j < dimension; ++j) {
                final double xj = workspace.reducedMean[j];
                final int reducedJ = partition.reducedIndexByTrait(j);
                final double yi = contribution.dLogL_dInformationY[i];
                final double yj = contribution.dLogL_dInformationY[j];

                final double exx = workspace.reducedCovariance[i * reducedDimension + j] + xi * xj;
                final double exy = tipObservation.observed[j]
                        ? xi * yj
                        : workspace.reducedCovariance[i * reducedDimension + reducedJ] + xi * yj;
                final double eyx = tipObservation.observed[i]
                        ? yi * xj
                        : workspace.reducedCovariance[reducedI * reducedDimension + j] + yi * xj;
                final double eyy = (tipObservation.observed[i] || tipObservation.observed[j])
                        ? yi * yj
                        : workspace.reducedCovariance[reducedI * reducedDimension + reducedJ] + yi * yj;

                final int ij = i * dimension + j;
                contribution.dLogL_dPrecisionXX[ij] = -0.5 * exx;
                contribution.dLogL_dPrecisionXY[ij] = -0.5 * exy;
                contribution.dLogL_dPrecisionYX[ij] = -0.5 * eyx;
                contribution.dLogL_dPrecisionYY[ij] = -0.5 * eyy;
            }
        }
        contribution.dLogL_dLogNormalizer = -1.0;
    }

    public void fillContributionForFixedParentPartiallyObservedTip(
            final CanonicalGaussianTransition transition,
            final CanonicalTipObservation tipObservation,
            final double[] fixedRootValue,
            final PartialIdentityTipObservationWorkspace workspace,
            final CanonicalGaussianState conditionedChildState,
            final CanonicalBranchMessageContribution contribution) {
        final TipObservationPartition partition = workspace.partition;
        final int observedCount = partition.update(tipObservation);
        final int missingCount = partition.getMissingCount();

        CanonicalGaussianMessageOps.conditionOnObservedFirstBlock(
                transition, fixedRootValue, conditionedChildState);

        for (int missing = 0; missing < missingCount; ++missing) {
            final int missingTrait = partition.missingIndex(missing);
            double info = conditionedChildState.information[missingTrait];
            final int missingTraitOff = missingTrait * dimension;
            for (int observed = 0; observed < observedCount; ++observed) {
                final int observedTrait = partition.observedIndex(observed);
                info -= conditionedChildState.precision[missingTraitOff + observedTrait]
                        * tipObservation.values[observedTrait];
            }
            workspace.reducedInformation[missing] = info;
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                workspace.reducedPrecision[missing * missingCount + otherMissing] =
                        conditionedChildState.precision[missingTraitOff + partition.missingIndex(otherMissing)];
            }
        }

        MatrixOps.safeInvertPrecision(
                workspace.reducedPrecision,
                workspace.reducedCovariance,
                missingCount);

        for (int missing = 0; missing < missingCount; ++missing) {
            double sum = 0.0;
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                sum += workspace.reducedCovariance[missing * missingCount + otherMissing]
                        * workspace.reducedInformation[otherMissing];
            }
            workspace.missingMean[missing] = sum;
        }

        clearContribution(contribution);
        for (int i = 0; i < dimension; ++i) {
            contribution.dLogL_dInformationX[i] = fixedRootValue[i];
            contribution.dLogL_dInformationY[i] = tipObservation.observed[i]
                    ? tipObservation.values[i]
                    : workspace.missingMean[partition.reducedIndexByTrait(i) - dimension];
        }

        for (int i = 0; i < dimension; ++i) {
            final double xi = fixedRootValue[i];
            final double yi = contribution.dLogL_dInformationY[i];
            final int missingI = tipObservation.observed[i] ? -1 : partition.reducedIndexByTrait(i) - dimension;
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                final double xj = fixedRootValue[j];
                final double yj = contribution.dLogL_dInformationY[j];
                final int missingJ = tipObservation.observed[j] ? -1 : partition.reducedIndexByTrait(j) - dimension;
                final double eyy = (missingI < 0 || missingJ < 0)
                        ? yi * yj
                        : workspace.reducedCovariance[missingI * missingCount + missingJ] + yi * yj;

                contribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * (xi * xj);
                contribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                contribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                contribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * eyy;
            }
        }
        contribution.dLogL_dLogNormalizer = -1.0;
    }

    private static void clearContribution(final CanonicalBranchMessageContribution contribution) {
        Arrays.fill(contribution.dLogL_dInformationX, 0.0);
        Arrays.fill(contribution.dLogL_dInformationY, 0.0);
        Arrays.fill(contribution.dLogL_dPrecisionXX, 0.0);
        Arrays.fill(contribution.dLogL_dPrecisionXY, 0.0);
        Arrays.fill(contribution.dLogL_dPrecisionYX, 0.0);
        Arrays.fill(contribution.dLogL_dPrecisionYY, 0.0);
        contribution.dLogL_dLogNormalizer = -1.0;
    }
}
