package dr.inference.timeseries.representation;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedTransitionCapability;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianBranchTransitionKernel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianBranchTransitionKernel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps;
import dr.inference.timeseries.core.TimeGrid;

/**
 * Small cache keyed by exact {@code dt} values for independent time series sharing one OU process.
 */
final class RepeatedDeltaCanonicalTransitionCache {

    private final GaussianBranchTransitionKernel kernel;
    private final CanonicalGaussianBranchTransitionKernel canonicalKernel;
    private final OUProcessModel processModel;
    private final CanonicalPreparedTransitionCapability preparedTransition;
    private final ThreadLocal<CanonicalBranchWorkspace> branchWorkspace;
    private final ThreadLocal<double[]> stationaryMeanScratch;
    private Entry[] entries;
    private int entryCount;
    private long momentRequests;
    private long momentBuilds;
    private long canonicalRequests;
    private long canonicalBuilds;
    private boolean dirty;

    RepeatedDeltaCanonicalTransitionCache(final GaussianBranchTransitionKernel kernel,
                                          final OUProcessModel processModel) {
        if (kernel == null) {
            throw new IllegalArgumentException("kernel must not be null");
        }
        this.kernel = kernel;
        this.canonicalKernel = kernel instanceof CanonicalGaussianBranchTransitionKernel
                ? (CanonicalGaussianBranchTransitionKernel) kernel
                : null;
        this.processModel = processModel;
        this.preparedTransition = processModel != null
                && processModel.getSelectionMatrixParameterization() instanceof CanonicalPreparedTransitionCapability
                ? (CanonicalPreparedTransitionCapability) processModel.getSelectionMatrixParameterization()
                : null;
        this.branchWorkspace = preparedTransition == null
                ? null
                : ThreadLocal.withInitial(preparedTransition::createBranchWorkspace);
        this.stationaryMeanScratch = preparedTransition == null
                ? null
                : ThreadLocal.withInitial(() -> new double[kernel.getStateDimension()]);
        this.entries = new Entry[0];
        this.entryCount = 0;
        this.dirty = true;
    }

    synchronized void prepareTimeGrid(final TimeGrid timeGrid) {
        for (int t = 0; t < timeGrid.getTimeCount() - 1; ++t) {
            final double dt = validatedDelta(timeGrid, t, t + 1);
            ensureEntry(dt);
        }
    }

    void fillTransitionMatrix(final double dt, final double[][] out) {
        final Entry entry = ensureEntry(dt);
        synchronized (entry) {
            ++momentRequests;
            ensureMoments(entry, dt);
            GaussianMatrixOps.copyFlatToMatrix(entry.transitionMatrixFlat, out, entry.dimension);
        }
    }

    void fillTransitionOffset(final double dt, final double[] out) {
        final Entry entry = ensureEntry(dt);
        synchronized (entry) {
            ++momentRequests;
            ensureMoments(entry, dt);
            System.arraycopy(entry.transitionOffset, 0, out, 0, entry.dimension);
        }
    }

    void fillTransitionCovariance(final double dt, final double[][] out) {
        final Entry entry = ensureEntry(dt);
        synchronized (entry) {
            ++momentRequests;
            ensureMoments(entry, dt);
            GaussianMatrixOps.copyFlatToMatrix(entry.transitionCovarianceFlat, out, entry.dimension);
        }
    }

    void fillCanonicalTransition(final double dt, final CanonicalGaussianTransition out) {
        if (canonicalKernel == null) {
            throw new UnsupportedOperationException("Kernel does not support canonical transitions");
        }
        final Entry entry = ensureEntry(dt);
        synchronized (entry) {
            ++canonicalRequests;
            ensureCanonical(entry, dt);
            out.copyFrom(entry.canonicalTransition);
        }
    }

    CanonicalPreparedBranchHandle getThreadPreparedCanonicalBranch(final double dt,
                                                                   final double[] stationaryMean) {
        if (preparedTransition == null) {
            return null;
        }
        final Entry entry = ensureEntry(dt);
        final ThreadPreparedBranch prepared = entry.threadPrepared.get();
        if (prepared.version != entry.version) {
            preparedTransition.prepareBranch(dt, stationaryMean, prepared.handle);
            prepared.version = entry.version;
            markClean();
        }
        return prepared.handle;
    }

    synchronized RepeatedDeltaCacheStatistics getStatistics() {
        return new RepeatedDeltaCacheStatistics(
                entryCount,
                momentRequests,
                momentBuilds,
                canonicalRequests,
                canonicalBuilds);
    }

    synchronized void makeDirty() {
        if (dirty) {
            return;
        }
        for (int i = 0; i < entryCount; ++i) {
            entries[i].makeDirty();
        }
        dirty = true;
    }

    private Entry ensureEntry(final double dt) {
        final long bits = Double.doubleToLongBits(dt);
        Entry[] localEntries = entries;
        final int localCount = entryCount;
        for (int i = 0; i < localCount; ++i) {
            if (localEntries[i].dtBits == bits) {
                return localEntries[i];
            }
        }
        return createEntry(bits);
    }

    private synchronized Entry createEntry(final long bits) {
        for (int i = 0; i < entryCount; ++i) {
            if (entries[i].dtBits == bits) {
                return entries[i];
            }
        }
        if (entryCount == entries.length) {
            final int newLength = entries.length == 0 ? 4 : entries.length * 2;
            final Entry[] expanded = new Entry[newLength];
            System.arraycopy(entries, 0, expanded, 0, entries.length);
            entries = expanded;
        }
        final Entry entry = new Entry(bits, kernel.getStateDimension(), preparedTransition);
        entries[entryCount++] = entry;
        return entry;
    }

    private void ensureCanonical(final Entry entry, final double dt) {
        if (entry.canonicalValid) {
            return;
        }
        ++canonicalBuilds;
        if (preparedTransition == null) {
            canonicalKernel.fillCanonicalTransition(dt, entry.canonicalTransition);
        } else {
            ensurePrepared(entry, dt);
            preparedTransition.fillCanonicalTransitionPrepared(
                    entry.prepared,
                    processModel.getDiffusionMatrix(),
                    branchWorkspace.get(),
                    entry.canonicalTransition);
        }
        entry.canonicalValid = true;
        markClean();
    }

    private void ensureMoments(final Entry entry, final double dt) {
        if (entry.momentsValid) {
            return;
        }
        ++momentBuilds;
        boolean filledPrepared = false;
        if (preparedTransition != null) {
            ensurePrepared(entry, dt);
            filledPrepared = preparedTransition.fillTransitionMomentsPreparedFlat(
                    entry.prepared,
                    processModel.getDiffusionMatrix(),
                    branchWorkspace.get(),
                    entry.transitionMatrixFlat,
                    entry.transitionOffset,
                    entry.transitionCovarianceFlat);
        }
        if (!filledPrepared) {
            kernel.fillTransitionMatrix(dt, entry.transitionMatrix);
            kernel.fillTransitionOffset(dt, entry.transitionOffset);
            kernel.fillTransitionCovariance(dt, entry.transitionCovariance);
            GaussianMatrixOps.copyMatrixToFlat(entry.transitionMatrix, entry.transitionMatrixFlat, entry.dimension);
            GaussianMatrixOps.copyMatrixToFlat(entry.transitionCovariance, entry.transitionCovarianceFlat, entry.dimension);
        }
        entry.momentsValid = true;
        markClean();
    }

    private void ensurePrepared(final Entry entry, final double dt) {
        if (entry.preparedValid) {
            return;
        }
        final double[] stationaryMean = stationaryMeanScratch.get();
        kernel.getInitialMean(stationaryMean);
        preparedTransition.prepareBranch(dt, stationaryMean, entry.prepared);
        entry.preparedValid = true;
    }

    private synchronized void markClean() {
        dirty = false;
    }

    static double validatedDelta(final TimeGrid timeGrid, final int fromIndex, final int toIndex) {
        final double dt = timeGrid.getDelta(fromIndex, toIndex);
        if (!(dt > 0.0)) {
            throw new IllegalArgumentException("Time increments must be strictly positive");
        }
        return dt;
    }

    private static final class Entry {
        final long dtBits;
        final int dimension;
        final CanonicalGaussianTransition canonicalTransition;
        final double[] transitionMatrixFlat;
        final double[] transitionOffset;
        final double[] transitionCovarianceFlat;
        final double[][] transitionMatrix;
        final double[][] transitionCovariance;
        final CanonicalPreparedBranchHandle prepared;
        final ThreadLocal<ThreadPreparedBranch> threadPrepared;
        boolean canonicalValid;
        boolean momentsValid;
        boolean preparedValid;
        volatile int version;

        Entry(final long dtBits,
              final int dimension,
              final CanonicalPreparedTransitionCapability preparedTransition) {
            this.dtBits = dtBits;
            this.dimension = dimension;
            this.canonicalTransition = new CanonicalGaussianTransition(dimension);
            this.transitionMatrixFlat = new double[dimension * dimension];
            this.transitionOffset = new double[dimension];
            this.transitionCovarianceFlat = new double[dimension * dimension];
            this.transitionMatrix = preparedTransition == null ? new double[dimension][dimension] : null;
            this.transitionCovariance = preparedTransition == null ? new double[dimension][dimension] : null;
            this.prepared = preparedTransition == null ? null : preparedTransition.createPreparedBranchHandle();
            this.threadPrepared = preparedTransition == null
                    ? null
                    : new ThreadLocal<ThreadPreparedBranch>() {
                        @Override
                        protected ThreadPreparedBranch initialValue() {
                            return new ThreadPreparedBranch(preparedTransition.createPreparedBranchHandle());
                        }
                    };
            this.canonicalValid = false;
            this.momentsValid = false;
            this.preparedValid = false;
            this.version = 0;
        }

        void makeDirty() {
            canonicalValid = false;
            momentsValid = false;
            preparedValid = false;
            ++version;
            if (prepared != null) {
                prepared.invalidateCovariance();
            }
        }
    }

    private static final class ThreadPreparedBranch {
        final CanonicalPreparedBranchHandle handle;
        int version;

        ThreadPreparedBranch(final CanonicalPreparedBranchHandle handle) {
            this.handle = handle;
            this.version = -1;
        }
    }
}
