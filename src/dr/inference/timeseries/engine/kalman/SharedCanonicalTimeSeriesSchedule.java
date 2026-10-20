package dr.inference.timeseries.engine.kalman;

import java.util.concurrent.atomic.AtomicLong;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;

/**
 * Shared, data-independent canonical branch schedule for independent time series.
 *
 * <p>The smoothed branch covariance blocks depend on the OU parameters, grid,
 * observation precision, and missingness pattern, but not on the observed values.
 * Parallel time-series likelihoods can therefore reuse the covariance part of the
 * RTS branch moments across independent series and only recompute the
 * data-dependent mean products.</p>
 *
 * <p><b>Locking.</b> {@link #makeDirty()} holds the schedule monitor and then takes each
 * {@link Entry} monitor to invalidate it.  Evaluation threads in
 * {@link #fillContributionFromRtsMoments} hold only the entry monitor and never take the
 * schedule monitor: the {@code dirty} flag they clear is volatile and the counters are
 * atomic, so the lock order is always schedule-then-entry and the two paths cannot
 * deadlock.  Cache validity never depends on {@code dirty}; every hit is re-validated
 * bitwise against the current covariances and a mismatch disables the schedule until the
 * next {@link #makeDirty()}.  The flag only lets {@link #makeDirty()} skip re-invalidating
 * entries that nobody has rebuilt since the previous invalidation.</p>
 */
public final class SharedCanonicalTimeSeriesSchedule {

    private final Entry[] entries;
    private volatile boolean enabled;
    /**
     * True when every entry has been invalidated and none rebuilt since.  Set under the
     * schedule monitor in {@link #makeDirty()}, cleared by evaluation threads without it.
     */
    private volatile boolean dirty;
    private final AtomicLong builds = new AtomicLong();
    private final AtomicLong hits = new AtomicLong();
    private final AtomicLong mismatches = new AtomicLong();

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
        enabled = true;
        if (dirty) {
            return;
        }
        // Raise the flag before touching the entries: a build that races with this sweep
        // clears it again, so the next makeDirty() repeats the full invalidation.
        dirty = true;
        for (int i = 0; i < entries.length; ++i) {
            entries[i].makeDirty();
        }
    }

    boolean fillContributionFromRtsMoments(final int branchIndex,
                                           final MomentBranchSmootherStats current,
                                           final MomentBranchSmootherStats next,
                                           final CanonicalBranchMessageContribution out) {
        if (!enabled || branchIndex < 0 || branchIndex >= entries.length) {
            return false;
        }
        final Entry entry = entries[branchIndex];
        synchronized (entry) {
            MatrixOps.matMul(
                    current.smootherGain,
                    next.smoothedCovariance,
                    entry.crossCovariance,
                    entry.stateDimension);
            if (!entry.momentCacheValid) {
                entry.cacheMomentCovariances(
                        current.smoothedCovariance,
                        entry.crossCovariance,
                        next.smoothedCovariance);
                builds.incrementAndGet();
                // Volatile write only; taking the schedule monitor here would invert the
                // schedule-then-entry lock order used by makeDirty().
                dirty = false;
            } else if (!entry.momentCovariancesMatch(
                    current.smoothedCovariance,
                    entry.crossCovariance,
                    next.smoothedCovariance)) {
                mismatches.incrementAndGet();
                enabled = false;
                return false;
            } else {
                hits.incrementAndGet();
            }
            entry.fillContribution(current.smoothedMean, next.smoothedMean, out);
            return true;
        }
    }

    public boolean isEnabled() {
        return enabled;
    }

    public long getBuildCount() {
        return builds.get();
    }

    public long getHitCount() {
        return hits.get();
    }

    public long getMismatchCount() {
        return mismatches.get();
    }

    private static final class Entry {
        final int stateDimension;
        final double[] covarianceXX;
        final double[] covarianceXY;
        final double[] covarianceYY;
        final double[] crossCovariance;
        boolean momentCacheValid;

        Entry(final int stateDimension) {
            this.stateDimension = stateDimension;
            final int matrixLength = stateDimension * stateDimension;
            this.covarianceXX = new double[matrixLength];
            this.covarianceXY = new double[matrixLength];
            this.covarianceYY = new double[matrixLength];
            this.crossCovariance = new double[matrixLength];
            this.momentCacheValid = false;
        }

        void makeDirty() {
            synchronized (this) {
                momentCacheValid = false;
            }
        }

        void cacheMomentCovariances(final double[] currentCovariance,
                                    final double[] crossCovariance,
                                    final double[] nextCovariance) {
            System.arraycopy(currentCovariance, 0, covarianceXX, 0, covarianceXX.length);
            System.arraycopy(crossCovariance, 0, covarianceXY, 0, covarianceXY.length);
            System.arraycopy(nextCovariance, 0, covarianceYY, 0, covarianceYY.length);
            momentCacheValid = true;
        }

        boolean momentCovariancesMatch(final double[] currentCovariance,
                                       final double[] crossCovariance,
                                       final double[] nextCovariance) {
            return arraysMatch(covarianceXX, currentCovariance)
                    && arraysMatch(covarianceXY, crossCovariance)
                    && arraysMatch(covarianceYY, nextCovariance);
        }

        void fillContribution(final double[] currentMean,
                              final double[] nextMean,
                              final CanonicalBranchMessageContribution out) {
            for (int i = 0; i < stateDimension; ++i) {
                out.dLogL_dInformationX[i] = currentMean[i];
                out.dLogL_dInformationY[i] = nextMean[i];
                for (int j = 0; j < stateDimension; ++j) {
                    final int ij = i * stateDimension + j;
                    out.dLogL_dPrecisionXX[ij] = -0.5
                            * (covarianceXX[ij] + currentMean[i] * currentMean[j]);
                    out.dLogL_dPrecisionXY[ij] = -0.5
                            * (covarianceXY[ij] + currentMean[i] * nextMean[j]);
                    out.dLogL_dPrecisionYX[ij] = -0.5
                            * (covarianceXY[j * stateDimension + i]
                            + nextMean[i] * currentMean[j]);
                    out.dLogL_dPrecisionYY[ij] = -0.5
                            * (covarianceYY[ij] + nextMean[i] * nextMean[j]);
                }
            }
            out.dLogL_dLogNormalizer = -1.0;
        }

        private static boolean arraysMatch(final double[] left,
                                           final double[] right) {
            for (int i = 0; i < left.length; ++i) {
                if (Double.doubleToLongBits(left[i]) != Double.doubleToLongBits(right[i])) {
                    return false;
                }
            }
            return true;
        }
    }
}
