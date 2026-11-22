package dr.evomodel.treedatalikelihood.continuousNew.canonical;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuousNew.CanonicalGradientFallbackPolicy;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientInputs;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.CanonicalBranchLengthGradientEngine;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.traversal.CanonicalTreeStateStore;
import dr.evomodel.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.CanonicalTipObservationModel;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.TipObservationMode;

import java.util.Arrays;

final class CanonicalBranchLengthGradientRunner {

    private static final double BRANCH_LENGTH_FD_RELATIVE_STEP = 1.0e-6;
    private static final double BRANCH_LENGTH_FD_ABSOLUTE_STEP = 1.0e-8;

    private final Tree tree;
    private final int nodeCount;
    private final CanonicalGradientFallbackPolicy fallbackPolicy;
    private final CanonicalTreeStateStore stateStore;
    private final BranchGradientWorkspace workspace;
    private final CanonicalTraversalRunner traversalRunner;
    private final CanonicalGradientPreparationRunner gradientPreparationRunner;
    private final BranchGradientInputs preparedBranchGradientInputs;
    private final CanonicalBranchLengthGradientEngine branchLengthGradientEngine;

    CanonicalBranchLengthGradientRunner(final Tree tree,
                                        final CanonicalGradientFallbackPolicy fallbackPolicy,
                                        final CanonicalTreeStateStore stateStore,
                                        final BranchGradientWorkspace workspace,
                                        final CanonicalTraversalRunner traversalRunner,
                                        final CanonicalGradientPreparationRunner gradientPreparationRunner,
                                        final BranchGradientInputs preparedBranchGradientInputs,
                                        final CanonicalBranchLengthGradientEngine branchLengthGradientEngine) {
        this.tree = tree;
        this.nodeCount = tree.getNodeCount();
        this.fallbackPolicy = fallbackPolicy;
        this.stateStore = stateStore;
        this.workspace = workspace;
        this.traversalRunner = traversalRunner;
        this.gradientPreparationRunner = gradientPreparationRunner;
        this.preparedBranchGradientInputs = preparedBranchGradientInputs;
        this.branchLengthGradientEngine = branchLengthGradientEngine;
    }

    void compute(final CanonicalBranchTransitionProvider transitionProvider,
                 final double[] gradT) {
        traversalRunner.requireGradientState();

        final CanonicalOUTransitionProvider ouProvider =
                CanonicalOUProviderSupport.requireOUProvider(transitionProvider);
        if (!fallbackPolicy.useBranchLengthFiniteDifference()) {
            gradientPreparationRunner.prepare(
                    transitionProvider,
                    CanonicalTransitionCachePhases.BRANCH_LENGTH_GRADIENT,
                    preparedBranchGradientInputs);
            computeAnalytic(ouProvider, preparedBranchGradientInputs, gradT);
            return;
        }

        try (CanonicalCachePhaseScope ignored = CanonicalCachePhaseScope.push(
                transitionProvider, CanonicalTransitionCachePhases.BRANCH_LENGTH_GRADIENT)) {
            computeFiniteDifference(transitionProvider, ouProvider, gradT);
        }
    }

    void compute(final CanonicalBranchTransitionProvider transitionProvider,
                 final BranchGradientInputs inputs,
                 final double[] gradT) {
        final CanonicalOUTransitionProvider ouProvider =
                CanonicalOUProviderSupport.requireOUProvider(transitionProvider);
        if (!fallbackPolicy.useBranchLengthFiniteDifference()) {
            computeAnalytic(ouProvider, inputs, gradT);
            return;
        }

        traversalRunner.requireGradientState();
        try (CanonicalCachePhaseScope ignored = CanonicalCachePhaseScope.push(
                transitionProvider, CanonicalTransitionCachePhases.BRANCH_LENGTH_GRADIENT)) {
            computeFiniteDifference(transitionProvider, ouProvider, gradT);
        }
    }

    private void computeAnalytic(final CanonicalOUTransitionProvider ouProvider,
                                 final BranchGradientInputs inputs,
                                 final double[] gradT) {
        branchLengthGradientEngine.compute(
                ouProvider,
                inputs,
                gradT);
    }

    private void computeFiniteDifference(final CanonicalBranchTransitionProvider transitionProvider,
                                         final CanonicalOUTransitionProvider ouProvider,
                                         final double[] gradT) {
        Arrays.fill(gradT, 0.0);
        final int rootIndex = tree.getRoot().getNumber();
        for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
            if (childIndex == rootIndex) {
                continue;
            }
            final CanonicalTipObservationModel observationModel =
                    tree.isExternal(tree.getNode(childIndex)) ? stateStore.tipObservationModels[childIndex] : null;
            if (observationModel != null && observationModel.isEmpty()) {
                gradT[childIndex] = 0.0;
                continue;
            }
            final double branchLength = transitionProvider.getEffectiveBranchLength(childIndex);
            gradT[childIndex] = finiteDifferenceBranchLengthGradient(childIndex, ouProvider, branchLength);
        }
    }

    private double finiteDifferenceBranchLengthGradient(final int childIndex,
                                                        final CanonicalOUTransitionProvider provider,
                                                        final double branchLength) {
        final double step = Math.max(BRANCH_LENGTH_FD_ABSOLUTE_STEP,
                BRANCH_LENGTH_FD_RELATIVE_STEP * Math.max(1.0, Math.abs(branchLength)));

        if (branchLength > step) {
            final double plus = evaluateFrozenLocalLogFactor(childIndex, provider, branchLength + step);
            final double minus = evaluateFrozenLocalLogFactor(childIndex, provider, branchLength - step);
            return (plus - minus) / (2.0 * step);
        }

        final double plus = evaluateFrozenLocalLogFactor(childIndex, provider, branchLength + step);
        final double base = evaluateFrozenLocalLogFactor(childIndex, provider, branchLength);
        return (plus - base) / step;
    }

    private double evaluateFrozenLocalLogFactor(final int childIndex,
                                                final CanonicalOUTransitionProvider provider,
                                                final double branchLength) {
        provider.fillCanonicalTransitionForBranchLength(childIndex, branchLength, workspace.transition);
        if (tree.isExternal(tree.getNode(childIndex))) {
            final CanonicalTipObservationModel observationModel = stateStore.tipObservationModels[childIndex];
            final CanonicalTransitionMomentProvider momentProvider =
                    observationModel.getMode() == TipObservationMode.PARTIAL_EXACT_IDENTITY
                            ? provider
                            : null;
            observationModel.fillParentMessage(
                    workspace.transition,
                    childIndex,
                    momentProvider,
                    branchLength,
                    workspace.tipParentMessageWorkspace,
                    workspace.gaussianWorkspace,
                    workspace.state);
        } else {
            CanonicalGaussianMessageOps.pushBackward(
                    stateStore.postOrder[childIndex],
                    workspace.transition,
                    workspace.gaussianWorkspace,
                    workspace.state);
        }
        CanonicalGaussianMessageOps.combineStates(
                stateStore.branchAboveParent[childIndex], workspace.state, workspace.combinedState);
        return CanonicalGaussianMessageOps.normalizationShift(workspace.combinedState, workspace.gaussianWorkspace);
    }
}
