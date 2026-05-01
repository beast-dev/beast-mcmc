package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.PartialIdentityTipContribution;

import java.util.Arrays;

/**
 * Builds branch-local canonical transition adjoints from current tree messages.
 */
final class CanonicalBranchContributionAssembler {

    private final Tree tree;
    private final int dim;
    final CanonicalTreeStateStore stateStore;
    private final BranchContributionStrategyFactory strategyFactory;
    private final PartialIdentityTipContribution partialIdentityContribution;

    CanonicalBranchContributionAssembler(final Tree tree,
                                         final int dim,
                                         final CanonicalTreeStateStore stateStore) {
        this.tree = tree;
        this.dim = dim;
        this.stateStore = stateStore;
        this.strategyFactory = new BranchContributionStrategyFactory(tree, dim, stateStore);
        this.partialIdentityContribution = new PartialIdentityTipContribution(dim);
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
                                                            final BranchGradientWorkspace workspace) {
        partialIdentityContribution.fillContributionForFixedParentPartiallyObservedTip(
                transition,
                tipObservation,
                stateStore.fixedRootValue,
                workspace.partialObservationWorkspace,
                workspace.state,
                workspace.contribution);
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
                                                 final BranchGradientWorkspace workspace) {
        partialIdentityContribution.fillContributionForPartiallyObservedTip(
                aboveState,
                transition,
                tipObservation,
                workspace.partialObservationWorkspace,
                workspace.contribution);
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
