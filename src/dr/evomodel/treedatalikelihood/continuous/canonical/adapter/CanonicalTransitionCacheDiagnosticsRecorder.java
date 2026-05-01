package dr.evomodel.treedatalikelihood.continuous.canonical.adapter;

import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTransitionCachePhases;

import java.util.concurrent.atomic.AtomicLong;

final class CanonicalTransitionCacheDiagnosticsRecorder {

    private static final String PHASE_OTHER = "other";

    private final boolean enabled;
    private final ThreadLocal<String> phase = new ThreadLocal<>();
    private final AtomicLong requests = new AtomicLong();
    private final AtomicLong hits = new AtomicLong();
    private final AtomicLong misses = new AtomicLong();
    private final AtomicLong clears = new AtomicLong();
    private final AtomicLong diffusionChangedClears = new AtomicLong();
    private final AtomicLong stationaryMeanChangedClears = new AtomicLong();
    private final AtomicLong selectionChangedClears = new AtomicLong();
    private final AtomicLong branchLengthChangedClears = new AtomicLong();
    private final AtomicLong modelChangedClears = new AtomicLong();
    private final AtomicLong restoreClears = new AtomicLong();
    private final AtomicLong explicitClears = new AtomicLong();
    private final AtomicLong transitionRebuilds = new AtomicLong();
    private final AtomicLong preparedBasisRebuilds = new AtomicLong();
    private final AtomicLong preparedCovarianceInvalidations = new AtomicLong();
    private final AtomicLong snapshotRefreshes = new AtomicLong();
    private final AtomicLong stores = new AtomicLong();
    private final AtomicLong restores = new AtomicLong();
    private final AtomicLong accepts = new AtomicLong();
    private final AtomicLong postOrderRequests = new AtomicLong();
    private final AtomicLong postOrderMisses = new AtomicLong();
    private final AtomicLong preOrderRequests = new AtomicLong();
    private final AtomicLong preOrderMisses = new AtomicLong();
    private final AtomicLong gradientPrepRequests = new AtomicLong();
    private final AtomicLong gradientPrepMisses = new AtomicLong();
    private final AtomicLong branchLengthRequests = new AtomicLong();
    private final AtomicLong branchLengthMisses = new AtomicLong();
    private final AtomicLong otherRequests = new AtomicLong();
    private final AtomicLong otherMisses = new AtomicLong();

    CanonicalTransitionCacheDiagnosticsRecorder(final boolean enabled) {
        this.enabled = enabled;
    }

    String pushPhase(final String nextPhase) {
        if (!enabled) {
            return null;
        }
        final String previous = phase.get();
        phase.set(nextPhase);
        return previous;
    }

    void popPhase(final String previous) {
        if (!enabled) {
            return;
        }
        if (previous == null) {
            phase.remove();
        } else {
            phase.set(previous);
        }
    }

    String currentPhase() {
        if (!enabled) {
            return PHASE_OTHER;
        }
        final String current = phase.get();
        return current == null ? PHASE_OTHER : current;
    }

    void recordRequest(final String currentPhase) {
        if (!enabled) {
            return;
        }
        requests.incrementAndGet();
        phaseRequestCounter(currentPhase).incrementAndGet();
    }

    void recordHit() {
        if (enabled) {
            hits.incrementAndGet();
        }
    }

    void recordMiss(final String currentPhase) {
        if (!enabled) {
            return;
        }
        misses.incrementAndGet();
        phaseMissCounter(currentPhase).incrementAndGet();
    }

    void recordClear(final CanonicalTransitionCacheInvalidationReason reason) {
        if (!enabled) {
            return;
        }
        clears.incrementAndGet();
        if (reason == CanonicalTransitionCacheInvalidationReason.DIFFUSION_CHANGED) {
            diffusionChangedClears.incrementAndGet();
        } else if (reason == CanonicalTransitionCacheInvalidationReason.STATIONARY_MEAN_CHANGED) {
            stationaryMeanChangedClears.incrementAndGet();
        } else if (reason == CanonicalTransitionCacheInvalidationReason.SELECTION_CHANGED) {
            selectionChangedClears.incrementAndGet();
        } else if (reason == CanonicalTransitionCacheInvalidationReason.BRANCH_LENGTH_CHANGED) {
            branchLengthChangedClears.incrementAndGet();
        } else if (reason == CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED) {
            modelChangedClears.incrementAndGet();
        } else if (reason == CanonicalTransitionCacheInvalidationReason.RESTORE_STATE) {
            restoreClears.incrementAndGet();
        } else {
            explicitClears.incrementAndGet();
        }
    }

    void recordSnapshotRefresh() {
        if (enabled) {
            snapshotRefreshes.incrementAndGet();
        }
    }

    void recordTransitionRebuild() {
        if (enabled) {
            transitionRebuilds.incrementAndGet();
        }
    }

    void recordPreparedBasisRebuild() {
        if (enabled) {
            preparedBasisRebuilds.incrementAndGet();
        }
    }

    void recordPreparedCovarianceInvalidation() {
        if (enabled) {
            preparedCovarianceInvalidations.incrementAndGet();
        }
    }

    void recordStore() {
        if (enabled) {
            stores.incrementAndGet();
        }
    }

    void recordRestore() {
        if (enabled) {
            restores.incrementAndGet();
        }
    }

    void recordAccept() {
        if (enabled) {
            accepts.incrementAndGet();
        }
    }

    void report(final String label) {
        if (!enabled) {
            return;
        }
        System.err.println("[canonical-transition-cache] " + label
                + " requests=" + requests.get()
                + " hits=" + hits.get()
                + " misses=" + misses.get()
                + " clears=" + clears.get()
                + " clearReasons=diffusion/mean/selection/branch/model/restore/explicit="
                + diffusionChangedClears.get()
                + "/" + stationaryMeanChangedClears.get()
                + "/" + selectionChangedClears.get()
                + "/" + branchLengthChangedClears.get()
                + "/"
                + modelChangedClears.get()
                + "/" + restoreClears.get()
                + "/" + explicitClears.get()
                + " rebuilds=transition/basis/covarianceInvalidation="
                + transitionRebuilds.get()
                + "/" + preparedBasisRebuilds.get()
                + "/" + preparedCovarianceInvalidations.get()
                + " snapshots=" + snapshotRefreshes.get()
                + " store/restore/accept=" + stores.get()
                + "/" + restores.get()
                + "/" + accepts.get()
                + " postorder=" + postOrderRequests.get()
                + "/" + postOrderMisses.get()
                + " preorder=" + preOrderRequests.get()
                + "/" + preOrderMisses.get()
                + " gradientPrep=" + gradientPrepRequests.get()
                + "/" + gradientPrepMisses.get()
                + " branchLength=" + branchLengthRequests.get()
                + "/" + branchLengthMisses.get()
                + " other=" + otherRequests.get()
                + "/" + otherMisses.get());
    }

    long phaseMisses(final String currentPhase) {
        if (!enabled) {
            return 0L;
        }
        return phaseMissCounter(currentPhase).get();
    }

    long phaseRequests(final String currentPhase) {
        if (!enabled) {
            return 0L;
        }
        return phaseRequestCounter(currentPhase).get();
    }

    long transitionRebuilds() {
        return enabled ? transitionRebuilds.get() : 0L;
    }

    long preparedBasisRebuilds() {
        return enabled ? preparedBasisRebuilds.get() : 0L;
    }

    long clearCount(final CanonicalTransitionCacheInvalidationReason reason) {
        if (!enabled) {
            return 0L;
        }
        if (reason == CanonicalTransitionCacheInvalidationReason.DIFFUSION_CHANGED) {
            return diffusionChangedClears.get();
        }
        if (reason == CanonicalTransitionCacheInvalidationReason.STATIONARY_MEAN_CHANGED) {
            return stationaryMeanChangedClears.get();
        }
        if (reason == CanonicalTransitionCacheInvalidationReason.SELECTION_CHANGED) {
            return selectionChangedClears.get();
        }
        if (reason == CanonicalTransitionCacheInvalidationReason.BRANCH_LENGTH_CHANGED) {
            return branchLengthChangedClears.get();
        }
        if (reason == CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED) {
            return modelChangedClears.get();
        }
        if (reason == CanonicalTransitionCacheInvalidationReason.RESTORE_STATE) {
            return restoreClears.get();
        }
        return explicitClears.get();
    }

    private AtomicLong phaseRequestCounter(final String currentPhase) {
        if (CanonicalTransitionCachePhases.POSTORDER.equals(currentPhase)) {
            return postOrderRequests;
        }
        if (CanonicalTransitionCachePhases.PREORDER.equals(currentPhase)) {
            return preOrderRequests;
        }
        if (CanonicalTransitionCachePhases.GRADIENT_PREP.equals(currentPhase)) {
            return gradientPrepRequests;
        }
        if (CanonicalTransitionCachePhases.BRANCH_LENGTH_GRADIENT.equals(currentPhase)) {
            return branchLengthRequests;
        }
        return otherRequests;
    }

    private AtomicLong phaseMissCounter(final String currentPhase) {
        if (CanonicalTransitionCachePhases.POSTORDER.equals(currentPhase)) {
            return postOrderMisses;
        }
        if (CanonicalTransitionCachePhases.PREORDER.equals(currentPhase)) {
            return preOrderMisses;
        }
        if (CanonicalTransitionCachePhases.GRADIENT_PREP.equals(currentPhase)) {
            return gradientPrepMisses;
        }
        if (CanonicalTransitionCachePhases.BRANCH_LENGTH_GRADIENT.equals(currentPhase)) {
            return branchLengthMisses;
        }
        return otherMisses;
    }
}
