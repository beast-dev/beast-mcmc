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

package dr.evomodel.treedatalikelihood.continuous.adapter;

import dr.evomodel.continuous.ou.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.CanonicalPreparedBranchHandle;
import dr.evomodel.continuous.ou.CanonicalPreparedTransitionCapability;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockPreparedBranchBasis;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockPreparedBranchHandle;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalPreparedBranchSnapshot;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;

/**
 * Per-branch canonical transition cache for homogeneous OU providers.
 */
final class CanonicalTransitionCache {

    interface BranchLengthProvider {
        double getEffectiveBranchLength(int childNodeIndex);
    }

    private final int dimension;
    private final OUProcessModel processModel;
    private final BranchLengthProvider branchLengthProvider;
    private final BranchCacheEntry[] entries;
    private final CanonicalPreparedTransitionCapability preparedTransition;
    private final ThreadLocal<CanonicalBranchWorkspace> specializedWorkspace;
    private final double[] stationaryMeanScratch;
    private final CanonicalTransitionCacheDiagnosticsRecorder diagnostics;
    private final CacheEpoch epoch;
    private boolean stationaryMeanSnapshotValid;

    CanonicalTransitionCache(final int dimension,
                             final int nodeCount,
                             final OUProcessModel processModel,
                             final CanonicalPreparedTransitionCapability preparedTransition,
                             final BranchLengthProvider branchLengthProvider,
                             final CanonicalTransitionCacheOptions options) {
        this.dimension = dimension;
        this.processModel = processModel;
        this.branchLengthProvider = branchLengthProvider;
        this.entries = new BranchCacheEntry[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            this.entries[i] = new BranchCacheEntry();
        }
        this.preparedTransition = preparedTransition;
        this.specializedWorkspace =
                preparedTransition == null
                        ? null
                        : ThreadLocal.withInitial(preparedTransition::createBranchWorkspace);
        this.stationaryMeanScratch = preparedTransition == null ? null : new double[dimension];
        this.diagnostics = new CanonicalTransitionCacheDiagnosticsRecorder(options.isDiagnosticsEnabled());
        this.epoch = new CacheEpoch();
        this.stationaryMeanSnapshotValid = false;
    }

    void fillTransition(final int childNodeIndex, final CanonicalGaussianTransition out) {
        out.copyFrom(ensureTransition(childNodeIndex));
    }

    CanonicalGaussianTransition getTransitionView(final int childNodeIndex) {
        return ensureTransition(childNodeIndex);
    }

    OrthogonalBlockPreparedBranchBasis
    getOrthogonalPreparedBranchBasis(final int childNodeIndex) {
        final CanonicalPreparedBranchHandle handle = getPreparedBranchHandle(childNodeIndex);
        if (handle instanceof OrthogonalBlockPreparedBranchHandle) {
            return ((OrthogonalBlockPreparedBranchHandle) handle).getBasis();
        }
        return null;
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

    private CanonicalGaussianTransition ensureTransition(final int childNodeIndex) {
        recordRequest();
        final BranchCacheEntry entry = entries[childNodeIndex];
        final double effectiveBranchLength = branchLengthProvider.getEffectiveBranchLength(childNodeIndex);
        if (!entry.isTransitionCurrent(effectiveBranchLength, epoch)) {
            recordMiss();
            diagnostics.recordTransitionRebuild();
            if (entry.transition == null) {
                entry.transition = new CanonicalGaussianTransition(dimension);
            }
            if (Double.doubleToLongBits(entry.effectiveBranchLength)
                    != Double.doubleToLongBits(effectiveBranchLength)) {
                entry.preparedBasisValid = false;
            }
            final CanonicalPreparedBranchHandle prepared =
                    fillCachedTransition(entry, effectiveBranchLength);
            entry.snapshot = new CanonicalPreparedBranchSnapshot(
                    childNodeIndex,
                    effectiveBranchLength,
                    entry.transition,
                    prepared);
            entry.markTransitionCurrent(effectiveBranchLength, epoch);
        } else {
            recordHit();
        }
        return entry.transition;
    }

    private CanonicalPreparedBranchHandle fillCachedTransition(final BranchCacheEntry entry,
                                                               final double effectiveBranchLength) {
        if (preparedTransition == null) {
            processModel.fillCanonicalTransition(effectiveBranchLength, entry.transition);
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

    private void ensureStationaryMeanSnapshot() {
        if (stationaryMeanSnapshotValid) {
            return;
        }
        processModel.getInitialMean(stationaryMeanScratch);
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

    private static final class BranchCacheEntry {
        private CanonicalGaussianTransition transition;
        private CanonicalPreparedBranchSnapshot snapshot;
        private CanonicalPreparedBranchHandle preparedBranchHandle;
        private double effectiveBranchLength;
        private boolean valid;
        private boolean preparedBasisValid;
        private long transitionSelectionVersion;
        private long transitionDiffusionVersion;
        private long transitionStationaryMeanVersion;
        private long transitionBranchLengthVersion;
        private long preparedSelectionVersion;
        private long preparedStationaryMeanVersion;
        private long preparedBranchLengthVersion;

        private boolean isTransitionCurrent(final double effectiveBranchLength,
                                            final CacheEpoch epoch) {
            return valid
                    && Double.doubleToLongBits(this.effectiveBranchLength)
                    == Double.doubleToLongBits(effectiveBranchLength)
                    && transitionSelectionVersion == epoch.selectionVersion
                    && transitionDiffusionVersion == epoch.diffusionVersion
                    && transitionStationaryMeanVersion == epoch.stationaryMeanVersion
                    && transitionBranchLengthVersion == epoch.branchLengthVersion;
        }

        private void markTransitionCurrent(final double effectiveBranchLength,
                                           final CacheEpoch epoch) {
            this.effectiveBranchLength = effectiveBranchLength;
            this.transitionSelectionVersion = epoch.selectionVersion;
            this.transitionDiffusionVersion = epoch.diffusionVersion;
            this.transitionStationaryMeanVersion = epoch.stationaryMeanVersion;
            this.transitionBranchLengthVersion = epoch.branchLengthVersion;
            this.valid = true;
        }

        private boolean isPreparedBasisCurrent(final CacheEpoch epoch) {
            return preparedBasisValid
                    && preparedSelectionVersion == epoch.selectionVersion
                    && preparedStationaryMeanVersion == epoch.stationaryMeanVersion
                    && preparedBranchLengthVersion == epoch.branchLengthVersion;
        }

        private void markPreparedBasisCurrent(final CacheEpoch epoch) {
            this.preparedSelectionVersion = epoch.selectionVersion;
            this.preparedStationaryMeanVersion = epoch.stationaryMeanVersion;
            this.preparedBranchLengthVersion = epoch.branchLengthVersion;
            this.preparedBasisValid = true;
        }

        private boolean invalidate(final CanonicalTransitionCacheInvalidationReason reason) {
            valid = false;
            if (preparedBranchHandle == null) {
                preparedBasisValid = false;
                return false;
            }
            if (reason == CanonicalTransitionCacheInvalidationReason.DIFFUSION_CHANGED) {
                preparedBranchHandle.invalidateCovariance();
                return true;
            }
            preparedBasisValid = false;
            preparedBranchHandle.invalidateCovariance();
            return true;
        }
    }

    private static final class CacheEpoch {
        private long selectionVersion;
        private long diffusionVersion;
        private long stationaryMeanVersion;
        private long branchLengthVersion;

        private void advance(final CanonicalTransitionCacheInvalidationReason reason) {
            if (reason == CanonicalTransitionCacheInvalidationReason.DIFFUSION_CHANGED) {
                diffusionVersion++;
            } else if (reason == CanonicalTransitionCacheInvalidationReason.STATIONARY_MEAN_CHANGED) {
                stationaryMeanVersion++;
            } else if (reason == CanonicalTransitionCacheInvalidationReason.SELECTION_CHANGED) {
                selectionVersion++;
            } else if (reason == CanonicalTransitionCacheInvalidationReason.BRANCH_LENGTH_CHANGED) {
                branchLengthVersion++;
            } else {
                selectionVersion++;
                diffusionVersion++;
                stationaryMeanVersion++;
                branchLengthVersion++;
            }
        }
    }

}
