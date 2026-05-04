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

package dr.evomodel.treedatalikelihood.continuous.canonical.adapter;

import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalOUKernel;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedTransitionCapability;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalPreparedBranchSnapshot;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTransitionCachePhases;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.util.TaskPool;

/**
 * Per-branch canonical transition cache for homogeneous OU providers.
 */
final class CanonicalTransitionCache {

    interface BranchLengthProvider {
        double getEffectiveBranchLength(int childNodeIndex);
    }

    private final int dimension;
    private final OUProcessModel processModel;
    private final CanonicalOUKernel kernel;
    private final BranchLengthProvider branchLengthProvider;
    private final CanonicalTransitionCacheEntry[] entries;
    private final CanonicalPreparedTransitionCapability preparedTransition;
    private final ThreadLocal<CanonicalBranchWorkspace> specializedWorkspace;
    private final double[] stationaryMeanScratch;
    private final CanonicalTransitionCacheDiagnosticsRecorder diagnostics;
    private final CanonicalTransitionCacheEpoch epoch;
    private volatile boolean stationaryMeanSnapshotValid;

    CanonicalTransitionCache(final int dimension,
                             final int nodeCount,
                             final OUProcessModel processModel,
                             final CanonicalPreparedTransitionCapability preparedTransition,
                             final BranchLengthProvider branchLengthProvider,
                             final CanonicalTransitionCacheOptions options) {
        this.dimension = dimension;
        this.processModel = processModel;
        this.kernel = processModel.getCanonicalKernel();
        this.branchLengthProvider = branchLengthProvider;
        this.entries = new CanonicalTransitionCacheEntry[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            this.entries[i] = new CanonicalTransitionCacheEntry();
        }
        this.preparedTransition = preparedTransition;
        this.specializedWorkspace =
                preparedTransition == null
                        ? null
                        : ThreadLocal.withInitial(preparedTransition::createBranchWorkspace);
        this.stationaryMeanScratch = preparedTransition == null ? null : new double[dimension];
        this.diagnostics = new CanonicalTransitionCacheDiagnosticsRecorder(options.isDiagnosticsEnabled());
        this.epoch = new CanonicalTransitionCacheEpoch();
        this.stationaryMeanSnapshotValid = false;
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
        if (preparedTransition != null) {
            ensureStationaryMeanSnapshot();
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
            stationaryMeanSnapshotValid = false;
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
        final double effectiveBranchLength = branchLengthProvider.getEffectiveBranchLength(childNodeIndex);
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
                    fillCachedTransition(entry, effectiveBranchLength);
            entry.updateSnapshot(childNodeIndex, effectiveBranchLength, entry.transition, prepared);
            entry.markTransitionCurrent(effectiveBranchLength, epoch);
        } else {
            recordHit();
        }
        return entry.transition;
    }

    private CanonicalPreparedBranchHandle fillCachedTransition(final CanonicalTransitionCacheEntry entry,
                                                               final double effectiveBranchLength) {
        if (preparedTransition == null) {
            kernel.fillCanonicalTransition(effectiveBranchLength, entry.transition);
            return null;
        }

        CanonicalPreparedBranchHandle prepared =
                entry.preparedBranchHandle;
        if (prepared == null) {
            prepared = preparedTransition.createPreparedBranchHandle();
            entry.preparedBranchHandle = prepared;
            entry.preparedBasisValid = false;
        }
        ensureStationaryMeanSnapshot();
        if (!entry.isPreparedBasisCurrent(epoch)) {
            diagnostics.recordPreparedBasisRebuild();
            preparedTransition.prepareBranch(effectiveBranchLength, stationaryMeanScratch, prepared);
            entry.markPreparedBasisCurrent(epoch);
        }
        preparedTransition.fillCanonicalTransitionPrepared(
                prepared,
                processModel.getDiffusionMatrix(),
                specializedWorkspace.get(),
                entry.transition);
        return prepared;
    }

    private synchronized void ensureStationaryMeanSnapshot() {
        if (stationaryMeanSnapshotValid) {
            return;
        }
        kernel.getInitialMean(stationaryMeanScratch);
        stationaryMeanSnapshotValid = true;
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
