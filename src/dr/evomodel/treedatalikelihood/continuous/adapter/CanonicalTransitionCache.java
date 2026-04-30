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
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.SpecializedCanonicalSelectionParameterization;
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
    private final SpecializedCanonicalSelectionParameterization specializedSelection;
    private final ThreadLocal<CanonicalBranchWorkspace> specializedWorkspace;
    private final double[] stationaryMeanScratch;
    private final CanonicalTransitionCacheDiagnosticsRecorder diagnostics;
    private boolean stationaryMeanSnapshotValid;

    CanonicalTransitionCache(final int dimension,
                             final int nodeCount,
                             final OUProcessModel processModel,
                             final SpecializedCanonicalSelectionParameterization specializedSelection,
                             final BranchLengthProvider branchLengthProvider,
                             final CanonicalTransitionCacheOptions options) {
        this.dimension = dimension;
        this.processModel = processModel;
        this.branchLengthProvider = branchLengthProvider;
        this.entries = new BranchCacheEntry[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            this.entries[i] = new BranchCacheEntry();
        }
        this.specializedSelection = specializedSelection;
        this.specializedWorkspace =
                specializedSelection == null
                        ? null
                        : ThreadLocal.withInitial(specializedSelection::createBranchWorkspace);
        this.stationaryMeanScratch = specializedSelection == null ? null : new double[dimension];
        this.diagnostics = new CanonicalTransitionCacheDiagnosticsRecorder(options.isDiagnosticsEnabled());
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
        stationaryMeanSnapshotValid = false;
        for (int i = 0; i < entries.length; i++) {
            entries[i].invalidate(reason);
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
        if (!entry.valid
                || Double.doubleToLongBits(entry.effectiveBranchLength)
                != Double.doubleToLongBits(effectiveBranchLength)) {
            recordMiss();
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
            entry.effectiveBranchLength = effectiveBranchLength;
            entry.valid = true;
        } else {
            recordHit();
        }
        return entry.transition;
    }

    private CanonicalPreparedBranchHandle fillCachedTransition(final BranchCacheEntry entry,
                                                               final double effectiveBranchLength) {
        if (specializedSelection == null) {
            processModel.fillCanonicalTransition(effectiveBranchLength, entry.transition);
            return null;
        }

        CanonicalPreparedBranchHandle prepared =
                entry.preparedBranchHandle;
        if (prepared == null) {
            prepared = specializedSelection.createPreparedBranchHandle();
            entry.preparedBranchHandle = prepared;
            entry.preparedBasisValid = false;
        }
        ensureStationaryMeanSnapshot();
        if (!entry.preparedBasisValid) {
            specializedSelection.prepareBranch(effectiveBranchLength, stationaryMeanScratch, prepared);
            entry.preparedBasisValid = true;
        }
        specializedSelection.fillCanonicalTransitionPrepared(
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

        private void invalidate(final CanonicalTransitionCacheInvalidationReason reason) {
            valid = false;
            if (preparedBranchHandle == null) {
                preparedBasisValid = false;
                return;
            }
            if (reason == CanonicalTransitionCacheInvalidationReason.DIFFUSION_CHANGED) {
                preparedBranchHandle.invalidateCovariance();
                return;
            }
            preparedBasisValid = false;
            preparedBranchHandle.invalidateCovariance();
        }
    }

}
