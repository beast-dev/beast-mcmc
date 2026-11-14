package dr.evomodel.continuous.canonical.contribution;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedTransitionCapability;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalTransitionMomentProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalOUTransitionProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalPreparedBranchSnapshot;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalPreparedBranchSnapshotProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientMessageBundle;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientMessageProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.CanonicalTipObservationModel;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.IdentityCanonicalTipObservationModel;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.LinearExactTipContribution;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.LinearExactTipObservationModel;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.PartialIdentityTipContribution;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.TipObservationMode;

import java.util.Arrays;

/**
 * Builds branch-local canonical transition adjoints from current tree messages.
 */
public final class CanonicalBranchContributionAssembler {

    private final int dim;
    private final BranchContributionStrategyFactory strategyFactory;
    private final PartialIdentityTipContribution partialIdentityContribution;
    private final LinearExactTipContribution linearExactContribution;

    public CanonicalBranchContributionAssembler(final int dim) {
        this.dim = dim;
        this.strategyFactory = new BranchContributionStrategyFactory();
        this.partialIdentityContribution = new PartialIdentityTipContribution(dim);
        this.linearExactContribution = new LinearExactTipContribution(dim);
    }

    public boolean fillLocalAdjointsForBranch(final int childIndex,
                                              final CanonicalBranchTransitionProvider transitionProvider,
                                              final BranchGradientMessageProvider messageProvider,
                                              final BranchGradientWorkspace workspace) {
        messageProvider.fillBranchMessages(childIndex, workspace.branchMessages);
        return fillLocalAdjointsForBranch(workspace.branchMessages, transitionProvider, workspace);
    }

    public boolean fillLocalAdjointsForBranch(final BranchGradientMessageBundle messages,
                                             final CanonicalBranchTransitionProvider transitionProvider,
                                             final BranchGradientWorkspace workspace) {
        final int childIndex = messages.getChildNodeIndex();
        workspace.adjoints.clearObservedTipValueAdjoint();
        workspace.adjoints.clearObservationBasisAdjoint();
        if (fillLinearExactAdjointsIfNeeded(messages, transitionProvider, workspace)) {
            return true;
        }
        final BranchContributionStrategy strategy = strategyFactory.select(messages);
        if (!strategy.contributes()) {
            return false;
        }
        final CanonicalGaussianTransition transition =
                transitionFor(childIndex, transitionProvider, workspace);
        if (!strategy.fillContribution(messages, transition, this, workspace)) {
            return false;
        }

        fillAdjointsFromContribution(childIndex, transitionProvider, transition, workspace);
        return true;
    }

    private boolean fillLinearExactAdjointsIfNeeded(
            final BranchGradientMessageBundle messages,
            final CanonicalBranchTransitionProvider transitionProvider,
            final BranchGradientWorkspace workspace) {
        if (!messages.isExternalChild()
                || messages.getTipObservationModel().getMode() != TipObservationMode.PARTIAL_LINEAR_EXACT) {
            return false;
        }
        final CanonicalTipObservationModel observationModel = messages.getTipObservationModel();
        if (!(observationModel instanceof LinearExactTipObservationModel)) {
            throw new IllegalStateException("Partial linear exact mode requires "
                    + LinearExactTipObservationModel.class.getSimpleName()
                    + " but found " + observationModel.getClass().getName());
        }
        if (!(transitionProvider instanceof CanonicalTransitionMomentProvider)) {
            throw new UnsupportedOperationException(
                    "Partial linear exact observations require a CanonicalTransitionMomentProvider.");
        }
        final LinearExactTipObservationModel linearTip =
                (LinearExactTipObservationModel) observationModel;
        final CanonicalTransitionMomentProvider momentProvider =
                (CanonicalTransitionMomentProvider) transitionProvider;
        final int childIndex = messages.getChildNodeIndex();
        final double branchLength = transitionProvider.getEffectiveBranchLength(childIndex);
        if (messages.isFixedParent()) {
            linearExactContribution.fillLocalAdjointsForFixedParentLinearTip(
                    linearTip,
                    childIndex,
                    momentProvider,
                    branchLength,
                    messages.getFixedRootValue(),
                    workspace);
        } else {
            linearExactContribution.fillLocalAdjointsForLinearTip(
                    messages.getAbovePreOrderAtParent(),
                    linearTip,
                    childIndex,
                    momentProvider,
                    branchLength,
                    workspace);
        }
        return true;
    }

    private void fillAdjointsFromContribution(final int childIndex,
                                              final CanonicalBranchTransitionProvider transitionProvider,
                                              final CanonicalGaussianTransition transition,
                                              final BranchGradientWorkspace workspace) {
        if (fillAdjointsFromPreparedMoments(childIndex, transitionProvider, transition, workspace)) {
            return;
        }
        if (fillAdjointsFromProviderMoments(childIndex, transitionProvider, transition, workspace)) {
            return;
        }
        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                transition,
                workspace.contribution,
                workspace.transitionAdjointWorkspace,
                workspace.adjoints);
    }

    private boolean fillAdjointsFromPreparedMoments(final int childIndex,
                                                    final CanonicalBranchTransitionProvider transitionProvider,
                                                    final CanonicalGaussianTransition transition,
                                                    final BranchGradientWorkspace workspace) {
        if (!(transitionProvider instanceof CanonicalOUTransitionProvider)
                || !(transitionProvider instanceof CanonicalPreparedBranchSnapshotProvider)) {
            return false;
        }

        final CanonicalPreparedBranchSnapshot snapshot =
                ((CanonicalPreparedBranchSnapshotProvider) transitionProvider)
                        .getPreparedBranchSnapshot(childIndex);
        if (snapshot == null || snapshot.getPreparedBranchHandle() == null) {
            return false;
        }

        final OUProcessModel processModel =
                ((CanonicalOUTransitionProvider) transitionProvider)
                        .getProcessModel(snapshot.getProcessIndex());
        if (!(processModel.getSelectionMatrixParameterization()
                instanceof CanonicalPreparedTransitionCapability)) {
            return false;
        }

        final CanonicalPreparedTransitionCapability preparedTransition =
                (CanonicalPreparedTransitionCapability) processModel.getSelectionMatrixParameterization();
        final CanonicalBranchWorkspace specializedWorkspace =
                workspace.ensureSpecializedBranchWorkspace(preparedTransition);
        if (!preparedTransition.fillTransitionMomentsPreparedFlat(
                snapshot.getPreparedBranchHandle(),
                processModel.getDiffusionMatrix(),
                specializedWorkspace,
                workspace.transitionMatrixFlat,
                workspace.mean2,
                workspace.covariance2)) {
            return false;
        }

        CanonicalTransitionAdjointUtils.fillFromMoments(
                transition.precisionYY,
                workspace.covariance2,
                workspace.transitionMatrixFlat,
                workspace.mean2,
                workspace.contribution,
                workspace.transitionAdjointWorkspace,
                workspace.adjoints);
        return true;
    }

    private boolean fillAdjointsFromProviderMoments(final int childIndex,
                                                    final CanonicalBranchTransitionProvider transitionProvider,
                                                    final CanonicalGaussianTransition transition,
                                                    final BranchGradientWorkspace workspace) {
        if (!(transitionProvider instanceof CanonicalTransitionMomentProvider)
                || workspace.transitionMatrixFlat == null
                || workspace.mean2 == null
                || workspace.covariance2 == null) {
            return false;
        }

        final CanonicalTransitionMomentProvider momentProvider =
                (CanonicalTransitionMomentProvider) transitionProvider;
        if (!momentProvider.useTransitionMomentsForBranchAdjoints(childIndex)) {
            return false;
        }

        momentProvider.fillTransitionMomentsFlat(
                childIndex,
                transitionProvider.getEffectiveBranchLength(childIndex),
                workspace.transitionMatrixFlat,
                workspace.mean2,
                workspace.covariance2);
        CanonicalTransitionAdjointUtils.fillFromMoments(
                transition.precisionYY,
                workspace.covariance2,
                workspace.transitionMatrixFlat,
                workspace.mean2,
                workspace.contribution,
                workspace.transitionAdjointWorkspace,
                workspace.adjoints);
        return true;
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
                                        final IdentityCanonicalTipObservationModel tipObservation,
                                        final BranchGradientWorkspace workspace) {
        for (int i = 0; i < dim; ++i) {
            double info = transition.informationX[i] + aboveState.information[i];
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                info -= transition.precisionXY[iOff + j] * tipObservation.valueAt(j);
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
            final double yi = tipObservation.valueAt(i);
            workspace.contribution.dLogL_dInformationX[i] = xi;
            workspace.contribution.dLogL_dInformationY[i] = yi;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = workspace.mean[j];
                final double yj = tipObservation.valueAt(j);
                final double exx = workspace.covariance[iOff + j] + xi * xj;
                workspace.contribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * exx;
                workspace.contribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                workspace.contribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                workspace.contribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * (yi * yj);
            }
        }
        workspace.contribution.dLogL_dLogNormalizer = -1.0;
        fillObservedTipValueAdjoint(transition, tipObservation, workspace.mean, workspace.adjoints);
    }

    void fillContributionForFixedParentObservedTip(final CanonicalGaussianTransition transition,
                                                   final IdentityCanonicalTipObservationModel tipObservation,
                                                   final double[] fixedRootValue,
                                                   final BranchGradientWorkspace workspace) {
        clearContribution(workspace);
        for (int i = 0; i < dim; ++i) {
            final double xi = fixedRootValue[i];
            final double yi = tipObservation.valueAt(i);
            workspace.contribution.dLogL_dInformationX[i] = xi;
            workspace.contribution.dLogL_dInformationY[i] = yi;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = fixedRootValue[j];
                final double yj = tipObservation.valueAt(j);
                workspace.contribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * (xi * xj);
                workspace.contribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                workspace.contribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                workspace.contribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * (yi * yj);
            }
        }
        workspace.contribution.dLogL_dLogNormalizer = -1.0;
        fillObservedTipValueAdjoint(transition, tipObservation, fixedRootValue, workspace.adjoints);
    }

    private void fillObservedTipValueAdjoint(final CanonicalGaussianTransition transition,
                                             final IdentityCanonicalTipObservationModel tipObservation,
                                             final double[] parentMean,
                                             final dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalLocalTransitionAdjoints adjoints) {
        for (int k = 0; k < dim; ++k) {
            final double yk = tipObservation.valueAt(k);
            adjoints.observedTipValue[k] = yk;
            double precisionYYy = 0.0;
            double precisionXYTransposeMean = 0.0;
            for (int i = 0; i < dim; ++i) {
                precisionYYy += transition.precisionYY[k * dim + i] * tipObservation.valueAt(i);
                precisionXYTransposeMean += transition.precisionXY[i * dim + k] * parentMean[i];
            }
            adjoints.dLogL_dObservedTipValue[k] =
                    -precisionXYTransposeMean - precisionYYy + transition.informationY[k];
        }
        adjoints.hasObservedTipValueAdjoint = true;
    }

    void fillContributionForFixedParentInternalNode(final CanonicalGaussianTransition transition,
                                                    final CanonicalGaussianState childMessage,
                                                    final double[] fixedRootValue,
                                                    final BranchGradientWorkspace workspace) {
        CanonicalGaussianMessageOps.conditionOnObservedFirstBlock(
                transition, fixedRootValue, workspace.state);
        CanonicalGaussianMessageOps.combineStates(workspace.state, childMessage, workspace.parentPosterior);
        workspace.adjoint.fillMomentsFromCanonical(
                workspace.parentPosterior, workspace.mean, workspace.covariance, dim);
        fillContributionFromFixedParentChildMoments(
                workspace.mean, workspace.covariance, fixedRootValue, workspace);
    }

    void fillContributionForFixedParentPartiallyObservedTip(final CanonicalGaussianTransition transition,
                                                            final IdentityCanonicalTipObservationModel tipObservation,
                                                            final double[] fixedRootValue,
                                                            final BranchGradientWorkspace workspace) {
        partialIdentityContribution.fillContributionForFixedParentPartiallyObservedTip(
                transition,
                tipObservation,
                fixedRootValue,
                workspace.partialObservationWorkspace,
                workspace.state,
                workspace.contribution);
    }

    private void fillContributionFromFixedParentChildMoments(final double[] childMean,
                                                             final double[] childCovariance,
                                                             final double[] fixedRootValue,
                                                             final BranchGradientWorkspace workspace) {
        clearContribution(workspace);
        for (int i = 0; i < dim; ++i) {
            final double xi = fixedRootValue[i];
            final double yi = childMean[i];
            workspace.contribution.dLogL_dInformationX[i] = xi;
            workspace.contribution.dLogL_dInformationY[i] = yi;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = fixedRootValue[j];
                final double yj = childMean[j];
                final double eyy = childCovariance[iOff + j] + yi * yj;
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
                                                 final IdentityCanonicalTipObservationModel tipObservation,
                                                 final BranchGradientWorkspace workspace) {
        partialIdentityContribution.fillContributionForPartiallyObservedTip(
                aboveState,
                transition,
                tipObservation,
                workspace.partialObservationWorkspace,
                workspace.contribution);
    }

    IdentityCanonicalTipObservationModel identityTipObservationModel(final BranchGradientMessageBundle messages) {
        final CanonicalTipObservationModel observationModel = messages.getTipObservationModel();
        if (observationModel instanceof IdentityCanonicalTipObservationModel) {
            return (IdentityCanonicalTipObservationModel) observationModel;
        }
        throw new IllegalStateException("Expected identity tip observation model for node "
                + messages.getChildNodeIndex()
                + " but found " + observationModel.getClass().getName());
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
