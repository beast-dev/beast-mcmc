package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.traversal.CanonicalTreeStateStore;
import dr.evomodel.treedatalikelihood.continuous.canonical.traversal.CanonicalTreeTraversal;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.util.CanonicalTraversalTimer;

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
            final long timingStart = CanonicalTraversalTimer.start();
            final long requestsBefore = CanonicalCachePhaseScope.transitionCacheRequests(
                    transitionProvider, CanonicalTransitionCachePhases.POSTORDER);
            final long missesBefore = CanonicalCachePhaseScope.transitionCacheMisses(
                    transitionProvider, CanonicalTransitionCachePhases.POSTORDER);
            try {
                return treeTraversal.computePostOrderLogLikelihood(
                        transitionProvider,
                        rootPrior,
                        stateStore,
                        workspace);
            } finally {
                final long requestsAfter = CanonicalCachePhaseScope.transitionCacheRequests(
                        transitionProvider, CanonicalTransitionCachePhases.POSTORDER);
                final long missesAfter = CanonicalCachePhaseScope.transitionCacheMisses(
                        transitionProvider, CanonicalTransitionCachePhases.POSTORDER);
                CanonicalTraversalTimer.finishPostorder(
                        timingStart,
                        requestsAfter - requestsBefore,
                        missesAfter - missesBefore);
            }
        }
    }

    void computePreOrder(final CanonicalBranchTransitionProvider transitionProvider,
                         final CanonicalRootPrior rootPrior) {
        try (CanonicalCachePhaseScope ignored = CanonicalCachePhaseScope.push(
                transitionProvider, CanonicalTransitionCachePhases.PREORDER)) {
            final long timingStart = CanonicalTraversalTimer.start();
            final long requestsBefore = CanonicalCachePhaseScope.transitionCacheRequests(
                    transitionProvider, CanonicalTransitionCachePhases.PREORDER);
            final long missesBefore = CanonicalCachePhaseScope.transitionCacheMisses(
                    transitionProvider, CanonicalTransitionCachePhases.PREORDER);
            try {
                treeTraversal.computePreOrder(
                        transitionProvider,
                        rootPrior,
                        stateStore,
                        workspace);
            } finally {
                final long requestsAfter = CanonicalCachePhaseScope.transitionCacheRequests(
                        transitionProvider, CanonicalTransitionCachePhases.PREORDER);
                final long missesAfter = CanonicalCachePhaseScope.transitionCacheMisses(
                        transitionProvider, CanonicalTransitionCachePhases.PREORDER);
                CanonicalTraversalTimer.finishPreorder(
                        timingStart,
                        requestsAfter - requestsBefore,
                        missesAfter - missesBefore);
            }
        }
    }

    void requireGradientState() {
        if (!stateStore.hasPostOrderState || !stateStore.hasPreOrderState) {
            throw new IllegalStateException(
                    "Canonical gradients require both computePostOrderLogLikelihood and computePreOrder to have been called.");
        }
    }
}
