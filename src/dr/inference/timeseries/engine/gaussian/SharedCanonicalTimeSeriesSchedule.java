package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;

/**
 * Shared, data-independent canonical branch schedule for independent time series.
 *
 * <p>The pair-state precision for a branch depends on the OU parameters, grid,
 * observation precision, and missingness pattern, but not on the observed values.
 * Parallel time-series likelihoods can therefore reuse the pair covariance blocks
 * across independent series and only recompute the data-dependent information
 * products.</p>
 */
public final class SharedCanonicalTimeSeriesSchedule {

    private final Entry[] entries;
    private volatile boolean enabled;
    private boolean dirty;
    private long builds;
    private long hits;
    private long mismatches;

    public SharedCanonicalTimeSeriesSchedule(final int timeCount,
                                             final int stateDimension) {
        if (timeCount < 1) {
            throw new IllegalArgumentException("timeCount must be positive");
        }
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be positive");
        }
        final int branchCount = Math.max(0, timeCount - 1);
        this.entries = new Entry[branchCount];
        for (int i = 0; i < branchCount; ++i) {
            entries[i] = new Entry(stateDimension);
        }
        this.enabled = true;
        this.dirty = true;
    }

    public synchronized void makeDirty() {
        if (dirty) {
            enabled = true;
            return;
        }
        enabled = true;
        for (int i = 0; i < entries.length; ++i) {
            entries[i].makeDirty();
        }
        dirty = true;
    }

    boolean fillContribution(final int branchIndex,
                             final CanonicalGaussianState pairState,
                             final CanonicalBranchMessageContributionUtils.Workspace workspace,
                             final CanonicalBranchMessageContribution out) {
        if (!enabled || branchIndex < 0 || branchIndex >= entries.length) {
            return false;
        }
        final Entry entry = entries[branchIndex];
        synchronized (entry) {
            if (!entry.cache.isValid()) {
                CanonicalBranchMessageContributionUtils.fillPairMomentCache(
                        pairState,
                        workspace,
                        entry.cache);
                ++builds;
                markClean();
            } else if (!CanonicalBranchMessageContributionUtils.cachedPrecisionMatches(pairState, entry.cache)) {
                ++mismatches;
                enabled = false;
                return false;
            } else {
                ++hits;
            }
            CanonicalBranchMessageContributionUtils.fillFromCachedPairMoments(
                    pairState,
                    entry.cache,
                    out);
            return true;
        }
    }

    private synchronized void markClean() {
        dirty = false;
    }

    public boolean isEnabled() {
        return enabled;
    }

    public synchronized long getBuildCount() {
        return builds;
    }

    public synchronized long getHitCount() {
        return hits;
    }

    public synchronized long getMismatchCount() {
        return mismatches;
    }

    private static final class Entry {
        final CanonicalBranchMessageContributionUtils.PairMomentCache cache;

        Entry(final int stateDimension) {
            this.cache = new CanonicalBranchMessageContributionUtils.PairMomentCache(stateDimension);
        }

        void makeDirty() {
            synchronized (this) {
                cache.makeDirty();
            }
        }
    }
}
