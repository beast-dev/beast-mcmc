package dr.evomodel.treedatalikelihood.continuous.canonical;

final class CanonicalCachePhaseScope implements AutoCloseable {

    private final CanonicalTransitionCacheDiagnostics diagnostics;
    private final String previousPhase;

    private CanonicalCachePhaseScope(final CanonicalBranchTransitionProvider transitionProvider,
                                     final String phase) {
        if (transitionProvider instanceof CanonicalTransitionCacheDiagnostics) {
            this.diagnostics = (CanonicalTransitionCacheDiagnostics) transitionProvider;
            this.previousPhase = diagnostics.pushDiagnosticPhase(phase);
        } else {
            this.diagnostics = null;
            this.previousPhase = null;
        }
    }

    static CanonicalCachePhaseScope push(final CanonicalBranchTransitionProvider transitionProvider,
                                         final String phase) {
        return new CanonicalCachePhaseScope(transitionProvider, phase);
    }

    static long transitionCacheMisses(final CanonicalBranchTransitionProvider transitionProvider,
                                      final String phase) {
        if (transitionProvider instanceof CanonicalTransitionCacheDiagnostics) {
            return ((CanonicalTransitionCacheDiagnostics) transitionProvider)
                    .getTransitionCacheMissCount(phase);
        }
        return 0L;
    }

    static long transitionCacheRequests(final CanonicalBranchTransitionProvider transitionProvider,
                                        final String phase) {
        if (transitionProvider instanceof CanonicalTransitionCacheDiagnostics) {
            return ((CanonicalTransitionCacheDiagnostics) transitionProvider)
                    .getTransitionCacheRequestCount(phase);
        }
        return 0L;
    }

    @Override
    public void close() {
        if (diagnostics != null) {
            diagnostics.popDiagnosticPhase(previousPhase);
        }
    }
}
