package dr.inference.timeseries.engine.kalman;

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
                ++builds;
                markClean();
            } else if (!entry.momentCovariancesMatch(
                    current.smoothedCovariance,
                    entry.crossCovariance,
                    next.smoothedCovariance)) {
                ++mismatches;
                enabled = false;
                return false;
            } else {
                ++hits;
            }
            entry.fillContribution(current.smoothedMean, next.smoothedMean, out);
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
