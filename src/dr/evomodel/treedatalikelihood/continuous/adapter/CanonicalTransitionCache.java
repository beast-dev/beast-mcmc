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

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockBranchGradientWorkspace;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockPreparedBranchBasis;
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
    private final OrthogonalBlockCanonicalParameterization orthogonalSelection;
    private final ThreadLocal<OrthogonalBlockBranchGradientWorkspace>
            orthogonalWorkspace;
    private final double[] stationaryMeanScratch;
    private final CanonicalTransitionCacheDiagnosticsRecorder diagnostics;

    CanonicalTransitionCache(final int dimension,
                             final int nodeCount,
                             final OUProcessModel processModel,
                             final OrthogonalBlockCanonicalParameterization orthogonalSelection,
                             final BranchLengthProvider branchLengthProvider,
                             final CanonicalTransitionCacheOptions options) {
        this.dimension = dimension;
        this.processModel = processModel;
        this.branchLengthProvider = branchLengthProvider;
        this.entries = new BranchCacheEntry[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            this.entries[i] = new BranchCacheEntry();
        }
        this.orthogonalSelection = orthogonalSelection;
        this.orthogonalWorkspace =
                orthogonalSelection == null
                        ? null
                        : ThreadLocal.withInitial(orthogonalSelection::createBranchGradientWorkspace);
        this.stationaryMeanScratch = orthogonalSelection == null ? null : new double[dimension];
        this.diagnostics = new CanonicalTransitionCacheDiagnosticsRecorder(options.isDiagnosticsEnabled());
    }

    void fillTransition(final int childNodeIndex, final CanonicalGaussianTransition out) {
        copyTransition(ensureTransition(childNodeIndex), out);
    }

    OrthogonalBlockPreparedBranchBasis
    getOrthogonalPreparedBranchBasis(final int childNodeIndex) {
        final CanonicalPreparedBranchSnapshot snapshot = getPreparedBranchSnapshot(childNodeIndex);
        return snapshot == null ? null : snapshot.getOrthogonalPreparedBasis();
    }

    CanonicalPreparedBranchSnapshot getPreparedBranchSnapshot(final int childNodeIndex) {
        ensureTransition(childNodeIndex);
        return entries[childNodeIndex].snapshot;
    }

    void clear(final CanonicalTransitionCacheInvalidationReason reason) {
        diagnostics.recordClear(reason);
        for (int i = 0; i < entries.length; i++) {
            entries[i].valid = false;
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
            final OrthogonalBlockPreparedBranchBasis prepared =
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

    private OrthogonalBlockPreparedBranchBasis fillCachedTransition(final BranchCacheEntry entry,
                                                                    final double effectiveBranchLength) {
        if (orthogonalSelection == null) {
            processModel.fillCanonicalTransition(effectiveBranchLength, entry.transition);
            return null;
        }

        OrthogonalBlockPreparedBranchBasis prepared =
                entry.orthogonalPreparedBasis;
        if (prepared == null) {
            prepared = orthogonalSelection.createPreparedBranchBasis();
            entry.orthogonalPreparedBasis = prepared;
        }
        processModel.getInitialMean(stationaryMeanScratch);
        orthogonalSelection.prepareBranchBasis(effectiveBranchLength, stationaryMeanScratch, prepared);
        orthogonalSelection.fillCanonicalTransitionPrepared(
                prepared,
                processModel.getDiffusionMatrix(),
                orthogonalWorkspace.get(),
                entry.transition);
        return prepared;
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

    private static void copyTransition(final CanonicalGaussianTransition source,
                                       final CanonicalGaussianTransition target) {
        final int dimension = source.getDimension();
        System.arraycopy(source.informationX, 0, target.informationX, 0, dimension);
        System.arraycopy(source.informationY, 0, target.informationY, 0, dimension);
        System.arraycopy(source.precisionXX, 0, target.precisionXX, 0, dimension * dimension);
        System.arraycopy(source.precisionXY, 0, target.precisionXY, 0, dimension * dimension);
        System.arraycopy(source.precisionYX, 0, target.precisionYX, 0, dimension * dimension);
        System.arraycopy(source.precisionYY, 0, target.precisionYY, 0, dimension * dimension);
        target.logNormalizer = source.logNormalizer;
    }

    private static final class BranchCacheEntry {
        private CanonicalGaussianTransition transition;
        private CanonicalPreparedBranchSnapshot snapshot;
        private OrthogonalBlockPreparedBranchBasis orthogonalPreparedBasis;
        private double effectiveBranchLength;
        private boolean valid;
    }

}
