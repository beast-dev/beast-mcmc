package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.traversal.CanonicalTreeStateStore;
import dr.evomodel.treedatalikelihood.continuous.canonical.traversal.CanonicalTreeTraversal;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;

final class CanonicalTraversalRunner {

    private final CanonicalTreeTraversal treeTraversal;
    private final CanonicalTreeStateStore stateStore;
    private final BranchGradientWorkspace workspace;

    CanonicalTraversalRunner(final CanonicalTreeTraversal treeTraversal,
                             final CanonicalTreeStateStore stateStore,
                             final BranchGradientWorkspace workspace) {
        this.treeTraversal = treeTraversal;
        this.stateStore = stateStore;
        this.workspace = workspace;
    }

    double computePostOrderLogLikelihood(final CanonicalBranchTransitionProvider transitionProvider,
                                         final CanonicalRootPrior rootPrior) {
        try (CanonicalCachePhaseScope ignored = CanonicalCachePhaseScope.push(
                transitionProvider, CanonicalTransitionCachePhases.POSTORDER)) {
            return treeTraversal.computePostOrderLogLikelihood(
                    transitionProvider,
                    rootPrior,
                    stateStore,
                    workspace);
        }
    }

    void computePreOrder(final CanonicalBranchTransitionProvider transitionProvider,
                         final CanonicalRootPrior rootPrior) {
        try (CanonicalCachePhaseScope ignored = CanonicalCachePhaseScope.push(
                transitionProvider, CanonicalTransitionCachePhases.PREORDER)) {
            treeTraversal.computePreOrder(
                    transitionProvider,
                    rootPrior,
                    stateStore,
                    workspace);
        }
    }

    void requireGradientState() {
        if (!stateStore.hasPostOrderState || !stateStore.hasPreOrderState) {
            throw new IllegalStateException(
                    "Canonical gradients require both computePostOrderLogLikelihood and computePreOrder to have been called.");
        }
    }
}
