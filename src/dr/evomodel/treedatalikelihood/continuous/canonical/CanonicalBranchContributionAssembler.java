package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;

import java.util.Arrays;

/**
 * Builds branch-local canonical transition adjoints from current tree messages.
 */
final class CanonicalBranchContributionAssembler {

    private final Tree tree;
    private final int dim;
    final CanonicalTreeStateStore stateStore;
    private final BranchContributionStrategyFactory strategyFactory;

    CanonicalBranchContributionAssembler(final Tree tree,
                                         final int dim,
                                         final CanonicalTreeStateStore stateStore) {
        this.tree = tree;
        this.dim = dim;
        this.stateStore = stateStore;
        this.strategyFactory = new BranchContributionStrategyFactory(tree, dim, stateStore);
    }

    boolean fillLocalAdjointsForBranch(final int childIndex,
                                       final CanonicalBranchTransitionProvider transitionProvider,
                                       final BranchGradientWorkspace workspace) {
        final BranchContributionStrategy strategy = strategyFactory.select(childIndex);
        if (!strategy.contributes()) {
            return false;
        }
        final CanonicalGaussianTransition transition =
                transitionFor(childIndex, transitionProvider, workspace);
        if (!strategy.fillContribution(childIndex, transition, this, workspace)) {
            return false;
        }

        fillAdjointsFromContribution(transition, workspace);
        return true;
    }

    private void fillAdjointsFromContribution(final CanonicalGaussianTransition transition,
                                              final BranchGradientWorkspace workspace) {
        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                transition,
                workspace.contribution,
                workspace.transitionAdjointWorkspace,
                workspace.adjoints);
    }

    private CanonicalGaussianTransition transitionFor(final int childIndex,
                                                      final CanonicalBranchTransitionProvider transitionProvider,
                                                      final BranchGradientWorkspace workspace) {
        final CanonicalGaussianTransition transitionView =
                transitionProvider.getCanonicalTransitionView(childIndex);
        if (transitionView != null) {
            return transitionView;
        }
        transitionProvider.fillCanonicalTransition(childIndex, workspace.transition);
        return workspace.transition;
    }

    int collectObservationPartition(final CanonicalTipObservation tipObservation,
                                    final BranchGradientWorkspace workspace) {
        int observedCount = 0;
        int missingCount = 0;

        for (int i = 0; i < dim; i++) {
            if (tipObservation.observed[i]) {
                workspace.observedIndexScratch[observedCount++] = i;
                workspace.reducedIndexByTraitScratch[i] = -1;
            } else {
                workspace.missingIndexScratch[missingCount++] = i;
                workspace.reducedIndexByTraitScratch[i] = dim + missingCount - 1;
            }
        }

        if (observedCount != tipObservation.observedCount || observedCount + missingCount != dim) {
            throw new UnsupportedOperationException(
                    "Canonical tip observation partition is inconsistent with observedCount.");
        }
        return observedCount;
    }

    void fillContributionForObservedTip(final CanonicalGaussianState aboveState,
                                        final CanonicalGaussianTransition transition,
                                        final CanonicalTipObservation tipObservation,
                                        final BranchGradientWorkspace workspace) {
        for (int i = 0; i < dim; ++i) {
            double info = transition.informationX[i] + aboveState.information[i];
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                info -= transition.precisionXY[iOff + j] * tipObservation.values[j];
                workspace.parentPosterior.precision[iOff + j] =
                        transition.precisionXX[iOff + j] + aboveState.precision[iOff + j];
            }
            workspace.parentPosterior.information[i] = info;
        }
        workspace.parentPosterior.logNormalizer = 0.0;

        workspace.adjoint.fillMomentsFromCanonical(
                workspace.parentPosterior, workspace.mean, workspace.covariance, dim);

        for (int i = 0; i < dim; ++i) {
            final double xi = workspace.mean[i];
            final double yi = tipObservation.values[i];
            workspace.contribution.dLogL_dInformationX[i] = xi;
            workspace.contribution.dLogL_dInformationY[i] = yi;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = workspace.mean[j];
                final double yj = tipObservation.values[j];
                final double exx = workspace.covariance[i][j] + xi * xj;
                workspace.contribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * exx;
                workspace.contribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                workspace.contribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                workspace.contribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * (yi * yj);
            }
        }
        workspace.contribution.dLogL_dLogNormalizer = -1.0;
    }

    void fillContributionForFixedParentObservedTip(final CanonicalGaussianTransition transition,
                                                   final CanonicalTipObservation tipObservation,
                                                   final BranchGradientWorkspace workspace) {
        clearContribution(workspace);
        for (int i = 0; i < dim; ++i) {
            final double xi = stateStore.fixedRootValue[i];
            final double yi = tipObservation.values[i];
            workspace.contribution.dLogL_dInformationX[i] = xi;
            workspace.contribution.dLogL_dInformationY[i] = yi;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = stateStore.fixedRootValue[j];
                final double yj = tipObservation.values[j];
                workspace.contribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * (xi * xj);
                workspace.contribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                workspace.contribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                workspace.contribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * (yi * yj);
            }
        }
        workspace.contribution.dLogL_dLogNormalizer = -1.0;
    }

    void fillContributionForFixedParentInternalNode(final CanonicalGaussianTransition transition,
                                                    final CanonicalGaussianState childMessage,
                                                    final BranchGradientWorkspace workspace) {
        CanonicalGaussianMessageOps.conditionOnObservedFirstBlock(
                transition, stateStore.fixedRootValue, workspace.state);
        CanonicalGaussianMessageOps.combineStates(workspace.state, childMessage, workspace.parentPosterior);
        workspace.adjoint.fillMomentsFromCanonical(
                workspace.parentPosterior, workspace.mean, workspace.covariance, dim);
        fillContributionFromFixedParentChildMoments(workspace.mean, workspace.covariance, workspace);
    }

    void fillContributionForFixedParentPartiallyObservedTip(final CanonicalGaussianTransition transition,
                                                            final CanonicalTipObservation tipObservation,
                                                            final int observedCount,
                                                            final BranchGradientWorkspace workspace) {
        final int missingCount = dim - observedCount;
        CanonicalGaussianMessageOps.conditionOnObservedFirstBlock(
                transition, stateStore.fixedRootValue, workspace.state);

        for (int missing = 0; missing < missingCount; ++missing) {
            final int missingTrait = workspace.missingIndexScratch[missing];
            double info = workspace.state.information[missingTrait];
            final int missingTraitOff = missingTrait * dim;
            for (int observed = 0; observed < observedCount; ++observed) {
                final int observedTrait = workspace.observedIndexScratch[observed];
                info -= workspace.state.precision[missingTraitOff + observedTrait] * tipObservation.values[observedTrait];
            }
            workspace.reducedInformationScratch[missing] = info;
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                workspace.reducedPrecisionFlatScratch[missing * missingCount + otherMissing] =
                        workspace.state.precision[missingTraitOff + workspace.missingIndexScratch[otherMissing]];
            }
        }

        MatrixOps.safeInvertPrecision(
                workspace.reducedPrecisionFlatScratch,
                workspace.reducedCovarianceFlatScratch,
                missingCount);

        for (int missing = 0; missing < missingCount; ++missing) {
            double sum = 0.0;
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                sum += workspace.reducedCovarianceFlatScratch[missing * missingCount + otherMissing]
                        * workspace.reducedInformationScratch[otherMissing];
            }
            workspace.mean2[missing] = sum;
        }

        clearContribution(workspace);
        for (int i = 0; i < dim; ++i) {
            workspace.contribution.dLogL_dInformationX[i] = stateStore.fixedRootValue[i];
            workspace.contribution.dLogL_dInformationY[i] = tipObservation.observed[i]
                    ? tipObservation.values[i]
                    : workspace.mean2[workspace.reducedIndexByTraitScratch[i] - dim];
        }

        for (int i = 0; i < dim; ++i) {
            final double xi = stateStore.fixedRootValue[i];
            final double yi = workspace.contribution.dLogL_dInformationY[i];
            final int missingI = tipObservation.observed[i] ? -1 : workspace.reducedIndexByTraitScratch[i] - dim;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = stateStore.fixedRootValue[j];
                final double yj = workspace.contribution.dLogL_dInformationY[j];
                final int missingJ = tipObservation.observed[j] ? -1 : workspace.reducedIndexByTraitScratch[j] - dim;
                final double eyy = (missingI < 0 || missingJ < 0)
                        ? yi * yj
                        : workspace.reducedCovarianceFlatScratch[missingI * missingCount + missingJ] + yi * yj;

                workspace.contribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * (xi * xj);
                workspace.contribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                workspace.contribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                workspace.contribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * eyy;
            }
        }
        workspace.contribution.dLogL_dLogNormalizer = -1.0;
    }

    private void fillContributionFromFixedParentChildMoments(final double[] childMean,
                                                             final double[][] childCovariance,
                                                             final BranchGradientWorkspace workspace) {
        clearContribution(workspace);
        for (int i = 0; i < dim; ++i) {
            final double xi = stateStore.fixedRootValue[i];
            final double yi = childMean[i];
            workspace.contribution.dLogL_dInformationX[i] = xi;
            workspace.contribution.dLogL_dInformationY[i] = yi;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = stateStore.fixedRootValue[j];
                final double yj = childMean[j];
                final double eyy = childCovariance[i][j] + yi * yj;
                workspace.contribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * (xi * xj);
                workspace.contribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                workspace.contribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                workspace.contribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * eyy;
            }
        }
        workspace.contribution.dLogL_dLogNormalizer = -1.0;
    }

    void fillContributionForPartiallyObservedTip(final CanonicalGaussianState aboveState,
                                                 final CanonicalGaussianTransition transition,
                                                 final CanonicalTipObservation tipObservation,
                                                 final int observedCount,
                                                 final BranchGradientWorkspace workspace) {
        final int missingCount = dim - observedCount;
        final int reducedDimension = dim + missingCount;

        for (int i = 0; i < dim; ++i) {
            double info = transition.informationX[i] + aboveState.information[i];
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                workspace.reducedPrecisionFlatScratch[i * reducedDimension + j] =
                        transition.precisionXX[iOff + j] + aboveState.precision[iOff + j];
            }
            for (int observed = 0; observed < observedCount; ++observed) {
                final int observedIndex = workspace.observedIndexScratch[observed];
                info -= transition.precisionXY[iOff + observedIndex] * tipObservation.values[observedIndex];
            }
            workspace.reducedInformationScratch[i] = info;
            for (int missing = 0; missing < missingCount; ++missing) {
                workspace.reducedPrecisionFlatScratch[i * reducedDimension + dim + missing] =
                        transition.precisionXY[iOff + workspace.missingIndexScratch[missing]];
            }
        }

        for (int missing = 0; missing < missingCount; ++missing) {
            final int childIndex = workspace.missingIndexScratch[missing];
            final int row = dim + missing;
            final int childOff = childIndex * dim;
            double info = transition.informationY[childIndex];
            for (int observed = 0; observed < observedCount; ++observed) {
                final int observedIndex = workspace.observedIndexScratch[observed];
                info -= transition.precisionYY[childOff + observedIndex] * tipObservation.values[observedIndex];
            }
            workspace.reducedInformationScratch[row] = info;
            for (int j = 0; j < dim; ++j) {
                workspace.reducedPrecisionFlatScratch[row * reducedDimension + j] =
                        transition.precisionYX[childOff + j];
            }
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                workspace.reducedPrecisionFlatScratch[row * reducedDimension + dim + otherMissing] =
                        transition.precisionYY[childOff + workspace.missingIndexScratch[otherMissing]];
            }
        }

        MatrixOps.safeInvertPrecision(
                workspace.reducedPrecisionFlatScratch,
                workspace.reducedCovarianceFlatScratch,
                reducedDimension);

        for (int i = 0; i < reducedDimension; ++i) {
            double sum = 0.0;
            for (int j = 0; j < reducedDimension; ++j) {
                sum += workspace.reducedCovarianceFlatScratch[i * reducedDimension + j]
                        * workspace.reducedInformationScratch[j];
            }
            workspace.reducedMeanScratch[i] = Double.isNaN(sum) ? 0.0 : sum;
        }

        clearContribution(workspace);
        for (int i = 0; i < dim; ++i) {
            workspace.contribution.dLogL_dInformationX[i] = workspace.reducedMeanScratch[i];
            workspace.contribution.dLogL_dInformationY[i] = tipObservation.observed[i]
                    ? tipObservation.values[i]
                    : workspace.reducedMeanScratch[workspace.reducedIndexByTraitScratch[i]];
        }

        for (int i = 0; i < dim; ++i) {
            final double xi = workspace.reducedMeanScratch[i];
            final int reducedI = workspace.reducedIndexByTraitScratch[i];
            for (int j = 0; j < dim; ++j) {
                final double xj = workspace.reducedMeanScratch[j];
                final int reducedJ = workspace.reducedIndexByTraitScratch[j];
                final double yi = workspace.contribution.dLogL_dInformationY[i];
                final double yj = workspace.contribution.dLogL_dInformationY[j];

                final double exx = workspace.reducedCovarianceFlatScratch[i * reducedDimension + j] + xi * xj;
                final double exy = tipObservation.observed[j]
                        ? xi * yj
                        : workspace.reducedCovarianceFlatScratch[i * reducedDimension + reducedJ] + xi * yj;
                final double eyx = tipObservation.observed[i]
                        ? yi * xj
                        : workspace.reducedCovarianceFlatScratch[reducedI * reducedDimension + j] + yi * xj;
                final double eyy = (tipObservation.observed[i] || tipObservation.observed[j])
                        ? yi * yj
                        : workspace.reducedCovarianceFlatScratch[reducedI * reducedDimension + reducedJ] + yi * yj;

                final int ij = i * dim + j;
                workspace.contribution.dLogL_dPrecisionXX[ij] = -0.5 * exx;
                workspace.contribution.dLogL_dPrecisionXY[ij] = -0.5 * exy;
                workspace.contribution.dLogL_dPrecisionYX[ij] = -0.5 * eyx;
                workspace.contribution.dLogL_dPrecisionYY[ij] = -0.5 * eyy;
            }
        }
        workspace.contribution.dLogL_dLogNormalizer = -1.0;
    }

    private void clearContribution(final BranchGradientWorkspace workspace) {
        Arrays.fill(workspace.contribution.dLogL_dInformationX, 0.0);
        Arrays.fill(workspace.contribution.dLogL_dInformationY, 0.0);
        Arrays.fill(workspace.contribution.dLogL_dPrecisionXX, 0.0);
        Arrays.fill(workspace.contribution.dLogL_dPrecisionXY, 0.0);
        Arrays.fill(workspace.contribution.dLogL_dPrecisionYX, 0.0);
        Arrays.fill(workspace.contribution.dLogL_dPrecisionYY, 0.0);
        workspace.contribution.dLogL_dLogNormalizer = -1.0;
    }
}
