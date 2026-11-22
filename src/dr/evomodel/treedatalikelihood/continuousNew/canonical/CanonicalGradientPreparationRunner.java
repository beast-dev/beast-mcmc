package dr.evomodel.treedatalikelihood.continuousNew.canonical;

import dr.evomodel.treedatalikelihood.continuousNew.CanonicalDebugOptions;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientInputs;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.CanonicalBranchAdjointPreparer;

final class CanonicalGradientPreparationRunner {

    private final CanonicalDebugOptions debugOptions;
    private final CanonicalTraversalRunner traversalRunner;
    private final CanonicalBranchAdjointPreparer branchAdjointPreparer;

    CanonicalGradientPreparationRunner(final CanonicalDebugOptions debugOptions,
                                       final CanonicalTraversalRunner traversalRunner,
                                       final CanonicalBranchAdjointPreparer branchAdjointPreparer) {
        this.debugOptions = debugOptions;
        this.traversalRunner = traversalRunner;
        this.branchAdjointPreparer = branchAdjointPreparer;
    }

    void prepare(final CanonicalBranchTransitionProvider transitionProvider,
                 final String phase,
                 final BranchGradientInputs out) {
        try (CanonicalCachePhaseScope ignored = CanonicalCachePhaseScope.push(transitionProvider, phase)) {
            traversalRunner.requireGradientState();
            final long missesBefore = CanonicalCachePhaseScope.transitionCacheMisses(transitionProvider, phase);
            branchAdjointPreparer.prepare(transitionProvider, out);
            assertNoGradientTransitionMisses(transitionProvider, phase, missesBefore);
        }
    }

    private void assertNoGradientTransitionMisses(final CanonicalBranchTransitionProvider transitionProvider,
                                                  final String phase,
                                                  final long missesBefore) {
        if (!debugOptions.isAssertNoGradientCacheMissesEnabled()) {
            return;
        }
        final long missesAfter = CanonicalCachePhaseScope.transitionCacheMisses(transitionProvider, phase);
        if (missesAfter != missesBefore) {
            throw new IllegalStateException(
                    "Canonical transition cache rebuilt " + (missesAfter - missesBefore)
                            + " branch transitions during " + phase + ".");
        }
    }
}
