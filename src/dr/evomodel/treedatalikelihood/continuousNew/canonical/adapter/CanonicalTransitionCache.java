/*
 * CanonicalTransitionCache.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter;

import dr.evomodel.continuous.ou.canonical.CanonicalOUKernel;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedTransitionSupport;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalPreparedBranchSnapshot;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalTransitionCachePhases;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.util.TaskPool;

/**
 * Per-branch canonical transition cache for OU providers.
 *
 * <p>The common homogeneous case is detected once at construction and keeps the
 * original one-kernel, one-prepared-transition hot path. The branch-process path
 * supports branch-specific process selection for one-segment branches; multi-segment
 * composition is intentionally rejected until transition composition is added with
 * dedicated tests.
 */
final class CanonicalTransitionCache {

    private final int dimension;
    private final CanonicalOUBranchProcessProvider branchProcess;
    private final boolean homogeneousFastPath;
    private final CanonicalOUKernel homogeneousKernel;
    private final CanonicalPreparedTransitionSupport homogeneousPreparedTransition;
    private final double[] homogeneousStationaryMeanScratch;
    private final OUProcessModel[] processModels;
    private final CanonicalPreparedTransitionSupport[] preparedTransitions;
    private final double[][] stationaryMeanScratchByProcess;
    private final boolean[] stationaryMeanSnapshotValidByProcess;
    private final CanonicalTransitionCacheEntry[] entries;
    private final CanonicalTransitionCacheDiagnosticsRecorder diagnostics;
    private final CanonicalTransitionCacheEpoch epoch;
    private volatile boolean homogeneousStationaryMeanSnapshotValid;

    CanonicalTransitionCache(final CanonicalOUBranchProcessProvider branchProcess,
                             final CanonicalTransitionCacheOptions options) {
        this.branchProcess = branchProcess;
        this.dimension = branchProcess.getDimension();
        final int nodeCount = branchProcess.getNodeCount();
        this.homogeneousFastPath = isHomogeneousSingleSegment(branchProcess);
        this.entries = new CanonicalTransitionCacheEntry[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            this.entries[i] = new CanonicalTransitionCacheEntry();
        }
        if (homogeneousFastPath) {
            final OUProcessModel processModel = branchProcess.getProcessModel(0);
            this.homogeneousKernel = processModel.getCanonicalKernel();
            this.homogeneousPreparedTransition =
                    CanonicalPreparedTransitionSupport.create(processModel);
            this.homogeneousStationaryMeanScratch =
                    homogeneousPreparedTransition == null ? null : new double[dimension];
            this.processModels = null;
            this.preparedTransitions = null;
            this.stationaryMeanScratchByProcess = null;
            this.stationaryMeanSnapshotValidByProcess = null;
        } else {
            this.homogeneousKernel = null;
            this.homogeneousPreparedTransition = null;
            this.homogeneousStationaryMeanScratch = null;
            final int processCount = branchProcess.getProcessModelCount();
            this.processModels = new OUProcessModel[processCount];
            this.preparedTransitions = new CanonicalPreparedTransitionSupport[processCount];
            this.stationaryMeanScratchByProcess = new double[processCount][];
            this.stationaryMeanSnapshotValidByProcess = new boolean[processCount];
            for (int i = 0; i < processCount; ++i) {
                final OUProcessModel processModel = branchProcess.getProcessModel(i);
                processModels[i] = processModel;
                preparedTransitions[i] = CanonicalPreparedTransitionSupport.create(processModel);
                stationaryMeanScratchByProcess[i] =
                        preparedTransitions[i] == null ? null : new double[dimension];
            }
        }
        this.diagnostics = new CanonicalTransitionCacheDiagnosticsRecorder(options.isDiagnosticsEnabled());
        this.epoch = new CanonicalTransitionCacheEpoch();
        this.homogeneousStationaryMeanSnapshotValid = false;
    }

    void fillTransition(final int childNodeIndex, final CanonicalGaussianTransition out) {
        out.copyFrom(ensureTransition(childNodeIndex));
    }

    CanonicalGaussianTransition getTransitionView(final int childNodeIndex) {
        return ensureTransition(childNodeIndex);
    }

    void preloadTransitions(final int rootNodeIndex,
                            final TaskPool taskPool,
                            final int chunkSize) {
        if (taskPool == null || taskPool.getNumThreads() <= 1) {
            return;
        }
        if (homogeneousPreparedTransition != null) {
            ensureHomogeneousStationaryMeanSnapshot();
        }
        taskPool.forkDynamic(
                entries.length,
                chunkSize,
                (childNodeIndex, thread) -> {
                    final String previousPhase = pushDiagnosticPhase(CanonicalTransitionCachePhases.POSTORDER);
                    try {
                        if (childNodeIndex != rootNodeIndex) {
                            ensureTransition(childNodeIndex);
                        }
                    } finally {
                        popDiagnosticPhase(previousPhase);
                    }
                });
    }

    CanonicalPreparedBranchHandle getPreparedBranchHandle(final int childNodeIndex) {
        final CanonicalPreparedBranchSnapshot snapshot = getPreparedBranchSnapshot(childNodeIndex);
        return snapshot == null ? null : snapshot.getPreparedBranchHandle();
    }

    CanonicalPreparedBranchSnapshot getPreparedBranchSnapshot(final int childNodeIndex) {
        ensureTransition(childNodeIndex);
        return entries[childNodeIndex].snapshot;
    }

    void clear(final CanonicalTransitionCacheInvalidationReason reason) {
        diagnostics.recordClear(reason);
        epoch.advance(reason);
        if (reason != CanonicalTransitionCacheInvalidationReason.DIFFUSION_CHANGED) {
            clearStationaryMeanSnapshots();
        }
        for (int i = 0; i < entries.length; i++) {
            if (entries[i].invalidate(reason)) {
                diagnostics.recordPreparedCovarianceInvalidation();
            }
        }
    }

    String pushDiagnosticPhase(final String phase) {
        return diagnostics.pushPhase(phase);
    }

    void popDiagnosticPhase(final String previous) {
        diagnostics.popPhase(previous);
    }

    void recordSnapshotRefresh() {
        diagnostics.recordSnapshotRefresh();
    }

    void recordStore() {
        diagnostics.recordStore();
    }

    void recordRestore() {
        diagnostics.recordRestore();
    }

    void recordAccept() {
        diagnostics.recordAccept();
    }

    void report(final String label) {
        diagnostics.report(label);
    }

    long getMissCount(final String phase) {
        return diagnostics.phaseMisses(phase);
    }

    long getRequestCount(final String phase) {
        return diagnostics.phaseRequests(phase);
    }

    long getTransitionRebuildCount() {
        return diagnostics.transitionRebuilds();
    }

    long getPreparedBranchRebuildCount() {
        return diagnostics.preparedBasisRebuilds();
    }

    long getClearCount(final String reason) {
        return diagnostics.clearCount(CanonicalTransitionCacheInvalidationReason.valueOf(reason));
    }

    private CanonicalGaussianTransition ensureTransition(final int childNodeIndex) {
        recordRequest();
        final CanonicalTransitionCacheEntry entry = entries[childNodeIndex];
        final double effectiveBranchLength = branchProcess.getEffectiveBranchLength(childNodeIndex);
        if (!entry.isTransitionCurrent(effectiveBranchLength, epoch)) {
            recordMiss();
            diagnostics.recordTransitionRebuild();
            if (entry.hasBranchLengthChanged(effectiveBranchLength)) {
                diagnostics.recordClear(CanonicalTransitionCacheInvalidationReason.BRANCH_LENGTH_CHANGED);
            }
            if (entry.transition == null) {
                entry.transition = new CanonicalGaussianTransition(dimension);
            }
            if (Double.doubleToLongBits(entry.effectiveBranchLength)
                    != Double.doubleToLongBits(effectiveBranchLength)) {
                entry.preparedBasisValid = false;
            }
            final CanonicalPreparedBranchHandle prepared =
                    fillCachedTransition(childNodeIndex, entry, effectiveBranchLength);
            entry.updateSnapshot(childNodeIndex, effectiveBranchLength, entry.transition, prepared);
            entry.markTransitionCurrent(effectiveBranchLength, epoch);
        } else {
            recordHit();
        }
        return entry.transition;
    }

    private CanonicalPreparedBranchHandle fillCachedTransition(final int childNodeIndex,
                                                               final CanonicalTransitionCacheEntry entry,
                                                               final double effectiveBranchLength) {
        if (homogeneousFastPath) {
            return fillCachedHomogeneousTransition(entry, effectiveBranchLength);
        }
        return fillCachedBranchProcessTransition(childNodeIndex, entry);
    }

    private CanonicalPreparedBranchHandle fillCachedHomogeneousTransition(
            final CanonicalTransitionCacheEntry entry,
            final double effectiveBranchLength) {
        if (homogeneousPreparedTransition == null) {
            homogeneousKernel.fillCanonicalTransition(effectiveBranchLength, entry.transition);
            entry.processIndex = 0;
            return null;
        }

        CanonicalPreparedBranchHandle prepared =
                entry.preparedBranchHandle;
        if (prepared == null) {
            prepared = homogeneousPreparedTransition.createPreparedBranchHandle();
            entry.preparedBranchHandle = prepared;
            entry.preparedBasisValid = false;
        }
        ensureHomogeneousStationaryMeanSnapshot();
        if (!entry.isPreparedBasisCurrent(epoch)) {
            diagnostics.recordPreparedBasisRebuild();
            homogeneousPreparedTransition.prepareBranch(
                    effectiveBranchLength,
                    homogeneousStationaryMeanScratch,
                    prepared);
            entry.markPreparedBasisCurrent(epoch);
        }
        homogeneousPreparedTransition.fillCanonicalTransitionPrepared(prepared, entry.transition);
        entry.processIndex = 0;
        return prepared;
    }

    private CanonicalPreparedBranchHandle fillCachedBranchProcessTransition(
            final int childNodeIndex,
            final CanonicalTransitionCacheEntry entry) {
        final int segmentCount = branchProcess.getSegmentCount(childNodeIndex);
        if (segmentCount != 1) {
            throw new UnsupportedOperationException(
                    "Canonical OU transition cache currently supports one branch segment; node "
                            + childNodeIndex + " has " + segmentCount + " segments.");
        }
        final int processIndex = branchProcess.getSegmentProcessIndex(childNodeIndex, 0);
        final double effectiveTime = branchProcess.getSegmentEffectiveTime(childNodeIndex, 0);
        final CanonicalPreparedTransitionSupport preparedTransition =
                preparedTransitions[processIndex];
        if (entry.processIndex != processIndex) {
            entry.preparedBranchHandle = null;
            entry.preparedBasisValid = false;
            entry.processIndex = processIndex;
        }
        if (preparedTransition == null) {
            processModels[processIndex].getCanonicalKernel().fillCanonicalTransition(
                    effectiveTime,
                    entry.transition);
            return null;
        }
        CanonicalPreparedBranchHandle prepared = entry.preparedBranchHandle;
        if (prepared == null) {
            prepared = preparedTransition.createPreparedBranchHandle();
            entry.preparedBranchHandle = prepared;
            entry.preparedBasisValid = false;
        }
        ensureStationaryMeanSnapshot(processIndex);
        if (!entry.isPreparedBasisCurrent(epoch)) {
            diagnostics.recordPreparedBasisRebuild();
            preparedTransition.prepareBranch(
                    effectiveTime,
                    stationaryMeanScratchByProcess[processIndex],
                    prepared);
            entry.markPreparedBasisCurrent(epoch);
        }
        preparedTransition.fillCanonicalTransitionPrepared(prepared, entry.transition);
        return prepared;
    }

    private synchronized void ensureHomogeneousStationaryMeanSnapshot() {
        if (homogeneousStationaryMeanScratch == null
                || homogeneousStationaryMeanSnapshotValid) {
            return;
        }
        homogeneousKernel.getInitialMean(homogeneousStationaryMeanScratch);
        homogeneousStationaryMeanSnapshotValid = true;
    }

    private synchronized void ensureStationaryMeanSnapshot(final int processIndex) {
        if (stationaryMeanSnapshotValidByProcess[processIndex]) {
            return;
        }
        processModels[processIndex].getCanonicalKernel()
                .getInitialMean(stationaryMeanScratchByProcess[processIndex]);
        stationaryMeanSnapshotValidByProcess[processIndex] = true;
    }

    private void clearStationaryMeanSnapshots() {
        homogeneousStationaryMeanSnapshotValid = false;
        if (stationaryMeanSnapshotValidByProcess != null) {
            for (int i = 0; i < stationaryMeanSnapshotValidByProcess.length; ++i) {
                stationaryMeanSnapshotValidByProcess[i] = false;
            }
        }
    }

    private static boolean isHomogeneousSingleSegment(
            final CanonicalOUBranchProcessProvider branchProcess) {
        if (branchProcess.getProcessModelCount() != 1) {
            return false;
        }
        for (int childNodeIndex = 0; childNodeIndex < branchProcess.getNodeCount(); ++childNodeIndex) {
            final int segmentCount = branchProcess.getSegmentCount(childNodeIndex);
            if (segmentCount == 0) {
                continue;
            }
            if (segmentCount != 1
                    || branchProcess.getSegmentProcessIndex(childNodeIndex, 0) != 0) {
                return false;
            }
        }
        return true;
    }

    private void recordRequest() {
        diagnostics.recordRequest(currentDiagnosticPhase());
    }

    private void recordHit() {
        diagnostics.recordHit();
    }

    private void recordMiss() {
        diagnostics.recordMiss(currentDiagnosticPhase());
    }

    private String currentDiagnosticPhase() {
        return diagnostics.currentPhase();
    }

}
