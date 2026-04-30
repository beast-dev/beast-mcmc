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

import java.util.concurrent.atomic.AtomicLong;

/**
 * Per-branch canonical transition cache for homogeneous OU providers.
 */
final class CanonicalTransitionCache {

    interface BranchLengthProvider {
        double getEffectiveBranchLength(int childNodeIndex);
    }

    private static final String PHASE_OTHER = "other";

    private final int dimension;
    private final OUProcessModel processModel;
    private final BranchLengthProvider branchLengthProvider;
    private final CanonicalGaussianTransition[] transitions;
    private final CanonicalPreparedBranchSnapshot[] snapshots;
    private final double[] cachedEffectiveBranchLength;
    private final boolean[] valid;
    private final OrthogonalBlockCanonicalParameterization orthogonalSelection;
    private final OrthogonalBlockPreparedBranchBasis[] orthogonalPreparedCache;
    private final ThreadLocal<OrthogonalBlockBranchGradientWorkspace>
            orthogonalWorkspace;
    private final double[] stationaryMeanScratch;
    private final CacheDiagnostics diagnostics;

    CanonicalTransitionCache(final int dimension,
                             final int nodeCount,
                             final OUProcessModel processModel,
                             final OrthogonalBlockCanonicalParameterization orthogonalSelection,
                             final BranchLengthProvider branchLengthProvider,
                             final boolean diagnosticsEnabled) {
        this.dimension = dimension;
        this.processModel = processModel;
        this.branchLengthProvider = branchLengthProvider;
        this.transitions = new CanonicalGaussianTransition[nodeCount];
        this.snapshots = new CanonicalPreparedBranchSnapshot[nodeCount];
        this.cachedEffectiveBranchLength = new double[nodeCount];
        this.valid = new boolean[nodeCount];
        this.orthogonalSelection = orthogonalSelection;
        this.orthogonalPreparedCache =
                orthogonalSelection == null
                        ? null
                        : new OrthogonalBlockPreparedBranchBasis[nodeCount];
        this.orthogonalWorkspace =
                orthogonalSelection == null
                        ? null
                        : ThreadLocal.withInitial(orthogonalSelection::createBranchGradientWorkspace);
        this.stationaryMeanScratch = orthogonalSelection == null ? null : new double[dimension];
        this.diagnostics = diagnosticsEnabled ? new CacheDiagnostics() : null;
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
        if (orthogonalSelection == null) {
            ensureTransition(childNodeIndex);
            return snapshots[childNodeIndex];
        }
        ensureTransition(childNodeIndex);
        return snapshots[childNodeIndex];
    }

    void clear() {
        recordClear();
        for (int i = 0; i < valid.length; i++) {
            valid[i] = false;
        }
    }

    String pushDiagnosticPhase(final String phase) {
        if (diagnostics == null) {
            return null;
        }
        final String previous = diagnostics.phase.get();
        diagnostics.phase.set(phase);
        return previous;
    }

    void popDiagnosticPhase(final String previous) {
        if (diagnostics == null) {
            return;
        }
        if (previous == null) {
            diagnostics.phase.remove();
        } else {
            diagnostics.phase.set(previous);
        }
    }

    void recordSnapshotRefresh() {
        if (diagnostics != null) {
            diagnostics.recordSnapshotRefresh();
        }
    }

    void recordStore() {
        if (diagnostics != null) {
            diagnostics.recordStore();
        }
    }

    void recordRestore() {
        if (diagnostics != null) {
            diagnostics.recordRestore();
        }
    }

    void recordAccept() {
        if (diagnostics != null) {
            diagnostics.recordAccept();
        }
    }

    void report(final String label) {
        if (diagnostics != null) {
            diagnostics.report(label);
        }
    }

    private CanonicalGaussianTransition ensureTransition(final int childNodeIndex) {
        recordRequest();
        final double effectiveBranchLength = branchLengthProvider.getEffectiveBranchLength(childNodeIndex);
        if (!valid[childNodeIndex]
                || Double.doubleToLongBits(cachedEffectiveBranchLength[childNodeIndex])
                != Double.doubleToLongBits(effectiveBranchLength)) {
            recordMiss();
            CanonicalGaussianTransition cached = transitions[childNodeIndex];
            if (cached == null) {
                cached = new CanonicalGaussianTransition(dimension);
                transitions[childNodeIndex] = cached;
            }
            final OrthogonalBlockPreparedBranchBasis prepared =
                    fillCachedTransition(childNodeIndex, effectiveBranchLength, cached);
            snapshots[childNodeIndex] = new CanonicalPreparedBranchSnapshot(
                    childNodeIndex,
                    effectiveBranchLength,
                    cached,
                    prepared);
            cachedEffectiveBranchLength[childNodeIndex] = effectiveBranchLength;
            valid[childNodeIndex] = true;
        } else {
            recordHit();
        }
        return transitions[childNodeIndex];
    }

    private OrthogonalBlockPreparedBranchBasis fillCachedTransition(final int childNodeIndex,
                                                                    final double effectiveBranchLength,
                                                                    final CanonicalGaussianTransition cached) {
        if (orthogonalSelection == null) {
            processModel.fillCanonicalTransition(effectiveBranchLength, cached);
            return null;
        }

        OrthogonalBlockPreparedBranchBasis prepared =
                orthogonalPreparedCache[childNodeIndex];
        if (prepared == null) {
            prepared = orthogonalSelection.createPreparedBranchBasis();
            orthogonalPreparedCache[childNodeIndex] = prepared;
        }
        processModel.getInitialMean(stationaryMeanScratch);
        orthogonalSelection.prepareBranchBasis(effectiveBranchLength, stationaryMeanScratch, prepared);
        orthogonalSelection.fillCanonicalTransitionPrepared(
                prepared,
                processModel.getDiffusionMatrix(),
                orthogonalWorkspace.get(),
                cached);
        return prepared;
    }

    private void recordRequest() {
        if (diagnostics != null) {
            diagnostics.recordRequest(currentDiagnosticPhase());
        }
    }

    private void recordHit() {
        if (diagnostics != null) {
            diagnostics.recordHit();
        }
    }

    private void recordMiss() {
        if (diagnostics != null) {
            diagnostics.recordMiss(currentDiagnosticPhase());
        }
    }

    private void recordClear() {
        if (diagnostics != null) {
            diagnostics.recordClear();
        }
    }

    private String currentDiagnosticPhase() {
        final String phase = diagnostics.phase.get();
        return phase == null ? PHASE_OTHER : phase;
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

    private static final class CacheDiagnostics {
        private final ThreadLocal<String> phase = new ThreadLocal<>();
        private final AtomicLong requests = new AtomicLong();
        private final AtomicLong hits = new AtomicLong();
        private final AtomicLong misses = new AtomicLong();
        private final AtomicLong clears = new AtomicLong();
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

        private void recordRequest(final String phase) {
            requests.incrementAndGet();
            phaseRequestCounter(phase).incrementAndGet();
        }

        private void recordHit() {
            hits.incrementAndGet();
        }

        private void recordMiss(final String phase) {
            misses.incrementAndGet();
            phaseMissCounter(phase).incrementAndGet();
        }

        private void recordClear() {
            clears.incrementAndGet();
        }

        private void recordSnapshotRefresh() {
            snapshotRefreshes.incrementAndGet();
        }

        private void recordStore() {
            stores.incrementAndGet();
        }

        private void recordRestore() {
            restores.incrementAndGet();
        }

        private void recordAccept() {
            accepts.incrementAndGet();
        }

        private AtomicLong phaseRequestCounter(final String phase) {
            if ("postorder".equals(phase)) {
                return postOrderRequests;
            }
            if ("preorder".equals(phase)) {
                return preOrderRequests;
            }
            if ("gradientPrep".equals(phase)) {
                return gradientPrepRequests;
            }
            if ("branchLengthGradient".equals(phase)) {
                return branchLengthRequests;
            }
            return otherRequests;
        }

        private AtomicLong phaseMissCounter(final String phase) {
            if ("postorder".equals(phase)) {
                return postOrderMisses;
            }
            if ("preorder".equals(phase)) {
                return preOrderMisses;
            }
            if ("gradientPrep".equals(phase)) {
                return gradientPrepMisses;
            }
            if ("branchLengthGradient".equals(phase)) {
                return branchLengthMisses;
            }
            return otherMisses;
        }

        private void report(final String label) {
            System.err.println("[canonical-transition-cache] " + label
                    + " requests=" + requests.get()
                    + " hits=" + hits.get()
                    + " misses=" + misses.get()
                    + " clears=" + clears.get()
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
    }
}
