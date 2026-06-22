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
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedTransitionSupport;
import dr.inference.model.AbstractModel;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumerics;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumericsOptions;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservationModel;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.IdentityCanonicalTipObservationModel;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.TipObservationMode;
import dr.util.TaskPool;

/**
 * Homogeneous OU branch-transition provider for the canonical tree pathway.
 */
public final class HomogeneousCanonicalOUBranchTransitionProvider extends AbstractModel
        implements CanonicalBranchTransitionProvider,
        CanonicalOUTransitionProvider,
        CanonicalPreparedBranchSnapshotProvider,
        CanonicalTransitionCacheDiagnostics,
        CanonicalTransitionPreloader,
        CanonicalNativeBranchMessageProvider {

    private static final String USE_NATIVE_DIAGONAL_MESSAGES_PROPERTY =
            "beast.experimental.canonicalUseNativeDiagonalOUMessages";

    private final int dimension;
    private final MultivariateDiffusionModel diffusionModel;
    private final MultivariateElasticModel elasticModel;
    private final MatrixParameterInterface driftMatrix;
    private final Parameter stationaryMean;
    private final MatrixParameter diffusionCovariance;
    private final CanonicalOUDiffusionSnapshot diffusionSnapshot;
    private final CanonicalBranchLengthSnapshot branchLengths;
    private final OUProcessModel processModel;
    private final CanonicalTransitionCache transitionCache;
    private final DiagonalBranchMessageCacheEntry[] diagonalMessageCache;

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
        super("homogeneousCanonicalOUTransitionProvider");
        this.elasticModel = elasticModel;
        this.diffusionModel = diffusionModel;
        addModel(diffusionModel);

        this.driftMatrix = elasticModel.getStrengthOfSelectionMatrixParameter();
        this.stationaryMean = stationaryMean;
        addVariable(driftMatrix);
        addVariable(stationaryMean);
        this.dimension = driftMatrix.getRowDimension();
        this.diffusionCovariance = new MatrixParameter("canonicalOuProvider.diffusion", dimension, dimension);
        this.diffusionSnapshot = new CanonicalOUDiffusionSnapshot(
                dimension,
                diffusionModel,
                diffusionCovariance);
        this.branchLengths = new CanonicalBranchLengthSnapshot(
                tree,
                rateModel,
                rateTransformation);

        final MatrixParameter initialCovariance = new MatrixParameter("canonicalOuProvider.initial", dimension, dimension);
        setIdentity(initialCovariance);

        this.processModel = new OUProcessModel(
                "canonicalOuProvider.process",
                dimension,
                driftMatrix,
                diffusionCovariance,
                stationaryMean,
                initialCovariance);
        addModel(processModel);

        this.transitionCache = new CanonicalTransitionCache(
                dimension,
                tree.getNodeCount(),
                processModel,
                CanonicalPreparedTransitionSupport.create(processModel),
                this::getEffectiveBranchLength,
                CanonicalTransitionCacheOptions.fromSystemProperties());
        this.diagonalMessageCache = driftMatrix instanceof DiagonalMatrix
                ? new DiagonalBranchMessageCacheEntry[tree.getNodeCount()]
                : null;

        refreshSnapshot();
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
        if (nativeDiagonalMessagesEnabled()) {
            return;
        }
        ensureCurrentSnapshot();
        transitionCache.preloadTransitions(rootNodeIndex, taskPool, chunkSize);
    }

    @Override
    public void fillCanonicalTransitionForLength(final double branchLength,
                                                 final CanonicalGaussianTransition out) {
        ensureCurrentSnapshot();
        processModel.fillCanonicalTransition(branchLength, out);
    }

    @Override
    public void fillTransitionMatrixFlat(final double branchLength, final double[] out) {
        ensureCurrentSnapshot();
        processModel.fillTransitionMatrixFlat(branchLength, out);
    }

    @Override
    public void fillTransitionOffset(final double branchLength, final double[] out) {
        ensureCurrentSnapshot();
        processModel.fillTransitionOffset(branchLength, out);
    }

    @Override
    public void fillTransitionCovarianceFlat(final double branchLength, final double[] out) {
        ensureCurrentSnapshot();
        processModel.fillTransitionCovarianceFlat(branchLength, out);
    }

    @Override
    public void fillTransitionMomentsFlat(final double branchLength,
                                          final double[] transitionMatrixOut,
                                          final double[] transitionOffsetOut,
                                          final double[] transitionCovarianceOut) {
        ensureCurrentSnapshot();
        processModel.fillTransitionMomentsFlat(
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
        branchLengths.beginPhase();
        return transitionCache.pushDiagnosticPhase(phase);
    }

    @Override
    public void popDiagnosticPhase(final String previous) {
        transitionCache.popDiagnosticPhase(previous);
        branchLengths.endPhase();
    }

    @Override
    public double getEffectiveBranchLength(final int childNodeIndex) {
        return branchLengths.getEffectiveBranchLength(childNodeIndex);
    }

    @Override
    public void fillTraitCovarianceFlat(final double[] out) {
        ensureCurrentSnapshot();
        diffusionSnapshot.fillCovarianceFlat(out);
    }

    public OUProcessModel getProcessModel() {
        ensureCurrentSnapshot();
        return processModel;
    }

    @Override
    public boolean fillBackwardMessageForState(final int childNodeIndex,
                                               final CanonicalGaussianState childState,
                                               final BranchGradientWorkspace workspace,
                                               final CanonicalGaussianState out) {
        final DiagonalBranchMessageCacheEntry entry = diagonalMessageEntry(childNodeIndex, workspace);
        if (entry == null) {
            return false;
        }
        CanonicalGaussianMessageOps.pushBackwardDiagonalTransition(
                childState,
                entry.actualizationDiagonal,
                entry.displacement,
                entry.precision,
                entry.logNormalizer,
                workspace.gaussianWorkspace,
                out);
        return true;
    }

    @Override
    public boolean fillBackwardMessageForTip(final int childNodeIndex,
                                             final CanonicalTipObservationModel observationModel,
                                             final BranchGradientWorkspace workspace,
                                             final CanonicalGaussianState out) {
        if (!(observationModel instanceof IdentityCanonicalTipObservationModel)) {
            return false;
        }
        if (observationModel.getMode() == TipObservationMode.MISSING) {
            CanonicalGaussianMessageOps.clearState(out);
            return true;
        }
        if (observationModel.getMode() != TipObservationMode.EXACT_IDENTITY) {
            return false;
        }
        final DiagonalBranchMessageCacheEntry entry = diagonalMessageEntry(childNodeIndex, workspace);
        if (entry == null) {
            return false;
        }
        final IdentityCanonicalTipObservationModel identity =
                (IdentityCanonicalTipObservationModel) observationModel;
        for (int i = 0; i < dimension; ++i) {
            workspace.diagonalObserved[i] = identity.valueAt(i);
        }
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlockDiagonalTransition(
                entry.actualizationDiagonal,
                entry.displacement,
                entry.precision,
                entry.logNormalizer,
                workspace.diagonalObserved,
                out);
        return true;
    }

    @Override
    public boolean fillForwardMessage(final int childNodeIndex,
                                      final CanonicalGaussianState parentState,
                                      final BranchGradientWorkspace workspace,
                                      final CanonicalGaussianState out) {
        final DiagonalBranchMessageCacheEntry entry = diagonalMessageEntry(childNodeIndex, workspace);
        if (entry == null) {
            return false;
        }
        CanonicalGaussianMessageOps.pushForwardDiagonalTransition(
                parentState,
                entry.actualizationDiagonal,
                entry.displacement,
                entry.precision,
                entry.logNormalizer,
                workspace.gaussianWorkspace,
                out);
        return true;
    }

    @Override
    public boolean bypassCanonicalTransitionPreload() {
        return nativeDiagonalMessagesEnabled();
    }

    public boolean usesNativeDiagonalBranchMessages() {
        return nativeDiagonalMessagesEnabled();
    }

    private void ensureCurrentSnapshot() {
        if (!dirty) {
            return;
        }
        refreshSnapshot();
        dirty = false;
        transitionCache.clear(dirtyReason);
        clearDiagonalMessageCache();
        dirtyReason = CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED;
    }

    private void refreshSnapshot() {
        transitionCache.recordSnapshotRefresh();
        diffusionSnapshot.refresh();
        processModel.fireModelChanged();
    }

    private static void setIdentity(final MatrixParameter matrix) {
        for (int i = 0; i < matrix.getRowDimension(); i++) {
            for (int j = 0; j < matrix.getColumnDimension(); j++) {
                matrix.setParameterValueQuietly(i, j, i == j ? 1.0 : 0.0);
            }
        }
        matrix.fireParameterChangedEvent();
    }

    @Override
    protected void handleModelChangedEvent(final Model model, final Object object, final int index) {
        if (model == diffusionModel) {
            markDirty(CanonicalTransitionCacheInvalidationReason.DIFFUSION_CHANGED);
        } else if (model == elasticModel) {
            markDirty(CanonicalTransitionCacheInvalidationReason.SELECTION_CHANGED);
        } else {
            markDirty(CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED);
        }
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable,
                                              final int index,
                                              final Parameter.ChangeType type) {
        if (variable == driftMatrix) {
            markDirty(CanonicalTransitionCacheInvalidationReason.SELECTION_CHANGED);
        } else if (variable == stationaryMean) {
            markDirty(CanonicalTransitionCacheInvalidationReason.STATIONARY_MEAN_CHANGED);
        } else {
            markDirty(CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED);
        }
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
        clearDiagonalMessageCache();
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

    private boolean nativeDiagonalMessagesEnabled() {
        return diagonalMessageCache != null
                && Boolean.parseBoolean(System.getProperty(USE_NATIVE_DIAGONAL_MESSAGES_PROPERTY, "true"));
    }

    private DiagonalBranchMessageCacheEntry diagonalMessageEntry(final int childNodeIndex,
                                                                final BranchGradientWorkspace workspace) {
        if (!nativeDiagonalMessagesEnabled()) {
            return null;
        }
        ensureCurrentSnapshot();
        final double branchLength = getEffectiveBranchLength(childNodeIndex);
        if (!(branchLength > 0.0) || !Double.isFinite(branchLength)) {
            return null;
        }
        DiagonalBranchMessageCacheEntry entry = diagonalMessageCache[childNodeIndex];
        if (entry == null) {
            entry = new DiagonalBranchMessageCacheEntry(dimension);
            diagonalMessageCache[childNodeIndex] = entry;
        }
        if (!entry.isCurrent(branchLength)) {
            rebuildDiagonalMessageEntry(branchLength, entry, workspace);
        }
        return entry;
    }

    private void rebuildDiagonalMessageEntry(final double branchLength,
                                             final DiagonalBranchMessageCacheEntry entry,
                                             final BranchGradientWorkspace workspace) {
        diffusionSnapshot.fillCovarianceFlat(workspace.diagonalBranchCovarianceScratch);
        final double[] oneMinusActualization = workspace.diagonalOneMinusActualization;
        for (int i = 0; i < dimension; ++i) {
            final double alpha = driftMatrix.getParameterValue(i, i);
            final double scaled = -alpha * branchLength;
            entry.actualizationDiagonal[i] = Math.exp(scaled);
            oneMinusActualization[i] = -Math.expm1(scaled);
            entry.displacement[i] = oneMinusActualization[i] * stationaryMeanValue(i);
        }

        for (int i = 0; i < dimension; ++i) {
            final double alphaI = driftMatrix.getParameterValue(i, i);
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                final double alphaJ = driftMatrix.getParameterValue(j, j);
                final double denominator = alphaI + alphaJ;
                final double diffusionCovariance =
                        workspace.diagonalBranchCovarianceScratch[iOff + j];
                final double covariance;
                if (denominator == 0.0) {
                    covariance = branchLength * diffusionCovariance;
                } else {
                    final double factor = oneMinusActualization[i]
                            + oneMinusActualization[j]
                            - oneMinusActualization[i] * oneMinusActualization[j];
                    covariance = diffusionCovariance * factor / denominator;
                }
                workspace.diagonalBranchCovariance[iOff + j] = covariance;
            }
        }

        final double logDet = CanonicalNumerics.copyAndInvertPositiveDefiniteFlat(
                workspace.diagonalBranchCovariance,
                dimension,
                workspace.diagonalBranchCovarianceScratch,
                entry.precision,
                workspace.diagonalBranchCholeskyScratch,
                workspace.diagonalBranchLowerInverseScratch,
                CanonicalNumericsOptions.OU_TREE);
        entry.logNormalizer = 0.5 * (dimension * MultivariateElasticModel.LOG2PI
                + logDet
                + quadratic(entry.precision, entry.displacement));
        entry.effectiveBranchLength = branchLength;
        entry.valid = true;
    }

    private double stationaryMeanValue(final int index) {
        return stationaryMean.getDimension() == 1
                ? stationaryMean.getParameterValue(0)
                : stationaryMean.getParameterValue(index);
    }

    private double quadratic(final double[] precision, final double[] value) {
        double out = 0.0;
        for (int i = 0; i < dimension; ++i) {
            double row = 0.0;
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                row += precision[iOff + j] * value[j];
            }
            out += value[i] * row;
        }
        return out;
    }

    private void clearDiagonalMessageCache() {
        if (diagonalMessageCache == null) {
            return;
        }
        for (DiagonalBranchMessageCacheEntry entry : diagonalMessageCache) {
            if (entry != null) {
                entry.valid = false;
            }
        }
    }

    private static final class DiagonalBranchMessageCacheEntry {
        private final double[] actualizationDiagonal;
        private final double[] displacement;
        private final double[] precision;
        private double effectiveBranchLength;
        private double logNormalizer;
        private boolean valid;

        private DiagonalBranchMessageCacheEntry(final int dimension) {
            this.actualizationDiagonal = new double[dimension];
            this.displacement = new double[dimension];
            this.precision = new double[dimension * dimension];
            this.effectiveBranchLength = Double.NaN;
            this.logNormalizer = Double.NaN;
            this.valid = false;
        }

        private boolean isCurrent(final double branchLength) {
            return valid
                    && Double.doubleToLongBits(effectiveBranchLength)
                    == Double.doubleToLongBits(branchLength);
        }
    }
}
