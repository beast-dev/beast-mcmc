/*
 * HomogeneousCanonicalOUBranchTransitionProvider.java
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

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalNativeBranchMessageProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalOUTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalPreparedBranchSnapshot;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalPreparedBranchSnapshotProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTransitionCacheDiagnostics;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTransitionPreloader;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.TipObservationMode;
import dr.util.TaskPool;

/**
 * Homogeneous OU branch-transition provider for the canonical tree pathway.
 */
public class HomogeneousCanonicalOUBranchTransitionProvider extends AbstractModel
        implements CanonicalBranchTransitionProvider,
        CanonicalOUTransitionProvider,
        CanonicalPreparedBranchSnapshotProvider,
        CanonicalTransitionCacheDiagnostics,
        CanonicalTransitionPreloader,
        CanonicalNativeBranchMessageProvider {

    private final int dimension;
    private final CanonicalOUBranchProcessProvider branchProcess;
    private final CanonicalTransitionCache transitionCache;
    private final CanonicalNativeBranchMessageProvider nativeBranchMessageProvider;

    private boolean dirty = false;
    private CanonicalTransitionCacheInvalidationReason dirtyReason =
            CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED;

    public HomogeneousCanonicalOUBranchTransitionProvider(final Tree tree,
                                                          final MultivariateElasticModel elasticModel,
                                                          final MultivariateDiffusionModel diffusionModel,
                                                          final Parameter stationaryMean,
                                                          final BranchRateModel rateModel) {
        this(tree, elasticModel, diffusionModel, stationaryMean, rateModel, null);
    }

    public HomogeneousCanonicalOUBranchTransitionProvider(final Tree tree,
                                                          final MultivariateElasticModel elasticModel,
                                                          final MultivariateDiffusionModel diffusionModel,
                                                          final Parameter stationaryMean,
                                                          final BranchRateModel rateModel,
                                                          final ContinuousRateTransformation rateTransformation) {
        this(new HomogeneousCanonicalOUBranchProcessProvider(
                tree,
                elasticModel,
                diffusionModel,
                stationaryMean,
                rateModel,
                rateTransformation));
    }

    private HomogeneousCanonicalOUBranchTransitionProvider(
            final CanonicalOUBranchProcessProvider branchProcess) {
        super("homogeneousCanonicalOUTransitionProvider");
        this.branchProcess = branchProcess;
        this.branchProcess.registerWith(this);
        this.dimension = branchProcess.getDimension();

        this.transitionCache = new CanonicalTransitionCache(
                branchProcess,
                CanonicalTransitionCacheOptions.fromSystemProperties());
        this.nativeBranchMessageProvider = createNativeBranchMessageProvider(branchProcess.getNodeCount());

        refreshSnapshot();
    }

    protected CanonicalNativeBranchMessageProvider createNativeBranchMessageProvider(final int nodeCount) {
        return branchProcess.createNativeBranchMessageProvider(nodeCount);
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public void fillCanonicalTransition(final int childNodeIndex, final CanonicalGaussianTransition out) {
        ensureCurrentSnapshot();
        transitionCache.fillTransition(childNodeIndex, out);
    }

    @Override
    public CanonicalGaussianTransition getCanonicalTransitionView(final int childNodeIndex) {
        ensureCurrentSnapshot();
        return transitionCache.getTransitionView(childNodeIndex);
    }

    @Override
    public void preloadCanonicalTransitions(final int rootNodeIndex,
                                            final TaskPool taskPool,
                                            final int chunkSize) {
        ensureCurrentSnapshot();
        transitionCache.preloadTransitions(rootNodeIndex, taskPool, chunkSize);
    }

    @Override
    public void fillCanonicalTransitionForLength(final double branchLength,
                                                 final CanonicalGaussianTransition out) {
        ensureCurrentSnapshot();
        branchProcess.fillCanonicalTransitionForLength(branchLength, out);
    }

    @Override
    public void fillTransitionMatrixFlat(final double branchLength, final double[] out) {
        ensureCurrentSnapshot();
        branchProcess.fillTransitionMatrixFlat(branchLength, out);
    }

    @Override
    public void fillTransitionOffset(final double branchLength, final double[] out) {
        ensureCurrentSnapshot();
        branchProcess.fillTransitionOffset(branchLength, out);
    }

    @Override
    public void fillTransitionCovarianceFlat(final double branchLength, final double[] out) {
        ensureCurrentSnapshot();
        branchProcess.fillTransitionCovarianceFlat(branchLength, out);
    }

    @Override
    public void fillTransitionMomentsFlat(final double branchLength,
                                          final double[] transitionMatrixOut,
                                          final double[] transitionOffsetOut,
                                          final double[] transitionCovarianceOut) {
        ensureCurrentSnapshot();
        branchProcess.fillTransitionMomentsFlat(
                branchLength,
                transitionMatrixOut,
                transitionOffsetOut,
                transitionCovarianceOut);
    }

    public CanonicalPreparedBranchHandle getPreparedBranchHandle(final int childNodeIndex) {
        ensureCurrentSnapshot();
        return transitionCache.getPreparedBranchHandle(childNodeIndex);
    }

    @Override
    public CanonicalPreparedBranchSnapshot getPreparedBranchSnapshot(final int childNodeIndex) {
        ensureCurrentSnapshot();
        return transitionCache.getPreparedBranchSnapshot(childNodeIndex);
    }

    @Override
    public String pushDiagnosticPhase(final String phase) {
        ensureCurrentSnapshot();
        branchProcess.beginPhase();
        return transitionCache.pushDiagnosticPhase(phase);
    }

    @Override
    public void popDiagnosticPhase(final String previous) {
        transitionCache.popDiagnosticPhase(previous);
        branchProcess.endPhase();
    }

    @Override
    public double getEffectiveBranchLength(final int childNodeIndex) {
        return branchProcess.getEffectiveBranchLength(childNodeIndex);
    }

    @Override
    public void fillTraitCovarianceFlat(final double[] out) {
        ensureCurrentSnapshot();
        branchProcess.fillTraitCovarianceFlat(out);
    }

    public OUProcessModel getProcessModel() {
        ensureCurrentSnapshot();
        return branchProcess.getProcessModel();
    }

    @Override
    public boolean fillBackwardMessageForState(final int childNodeIndex,
                                               final CanonicalGaussianState childState,
                                               final BranchGradientWorkspace workspace,
                                               final CanonicalGaussianState out) {
        if (nativeBranchMessageProvider == null) {
            return false;
        }
        ensureCurrentSnapshot();
        return nativeBranchMessageProvider.fillBackwardMessageForState(
                childNodeIndex,
                childState,
                workspace,
                out);
    }

    @Override
    public boolean fillBackwardMessageForObservedTip(final int childNodeIndex,
                                                     final double[] observedValue,
                                                     final CanonicalGaussianState out) {
        if (nativeBranchMessageProvider == null) {
            return false;
        }
        ensureCurrentSnapshot();
        return nativeBranchMessageProvider.fillBackwardMessageForObservedTip(
                childNodeIndex,
                observedValue,
                out);
    }

    @Override
    public boolean fillForwardMessage(final int childNodeIndex,
                                      final CanonicalGaussianState parentState,
                                      final BranchGradientWorkspace workspace,
                                      final CanonicalGaussianState out) {
        if (nativeBranchMessageProvider == null) {
            return false;
        }
        ensureCurrentSnapshot();
        return nativeBranchMessageProvider.fillForwardMessage(
                childNodeIndex,
                parentState,
                workspace,
                out);
    }

    @Override
    public boolean hasNativeBranchMessages() {
        return nativeBranchMessageProvider != null
                && nativeBranchMessageProvider.hasNativeBranchMessages();
    }

    @Override
    public boolean supportsPostorderInternalMessage() {
        return nativeBranchMessageProvider != null
                && nativeBranchMessageProvider.supportsPostorderInternalMessage();
    }

    @Override
    public boolean supportsPostorderTipMessage(final TipObservationMode mode) {
        return nativeBranchMessageProvider != null
                && nativeBranchMessageProvider.supportsPostorderTipMessage(mode);
    }

    @Override
    public boolean supportsPreorderForwardMessage() {
        return nativeBranchMessageProvider != null
                && nativeBranchMessageProvider.supportsPreorderForwardMessage();
    }

    @Override
    public boolean supportsFixedRootForwardMessage() {
        return nativeBranchMessageProvider != null
                && nativeBranchMessageProvider.supportsFixedRootForwardMessage();
    }

    public boolean usesNativeBranchMessages() {
        return hasNativeBranchMessages();
    }

    private void ensureCurrentSnapshot() {
        if (!dirty) {
            return;
        }
        refreshSnapshot();
        dirty = false;
        transitionCache.clear(dirtyReason);
        clearNativeBranchMessageCache();
        dirtyReason = CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED;
    }

    private void refreshSnapshot() {
        transitionCache.recordSnapshotRefresh();
        branchProcess.refreshSnapshot();
    }

    @Override
    protected void handleModelChangedEvent(final Model model, final Object object, final int index) {
        markDirty(branchProcess.classifyModelChange(model));
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable,
                                              final int index,
                                              final Parameter.ChangeType type) {
        markDirty(branchProcess.classifyVariableChange(variable, type));
    }

    private void markDirty(final CanonicalTransitionCacheInvalidationReason reason) {
        dirty = true;
        if (reason == CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED
                && dirtyReason != CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED) {
            return;
        }
        dirtyReason = reason;
    }

    @Override
    public void storeState() {
        transitionCache.recordStore();
    }

    @Override
    public void restoreState() {
        transitionCache.recordRestore();
        transitionCache.clear(CanonicalTransitionCacheInvalidationReason.RESTORE_STATE);
        clearNativeBranchMessageCache();
        dirty = true;
        dirtyReason = CanonicalTransitionCacheInvalidationReason.RESTORE_STATE;
    }

    @Override
    public void acceptState() {
        transitionCache.recordAccept();
    }

    public void reportTransitionCacheDiagnostics(final String label) {
        transitionCache.report(label);
    }

    @Override
    public long getTransitionCacheMissCount(final String phase) {
        return transitionCache.getMissCount(phase);
    }

    @Override
    public long getTransitionCacheRequestCount(final String phase) {
        return transitionCache.getRequestCount(phase);
    }

    @Override
    public long getTransitionCacheRebuildCount() {
        return transitionCache.getTransitionRebuildCount();
    }

    @Override
    public long getPreparedBranchRebuildCount() {
        return transitionCache.getPreparedBranchRebuildCount();
    }

    @Override
    public long getTransitionCacheClearCount(final String reason) {
        return transitionCache.getClearCount(reason);
    }

    private void clearNativeBranchMessageCache() {
        if (nativeBranchMessageProvider != null) {
            nativeBranchMessageProvider.clearBranchMessageCache();
        }
    }
}
