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

package dr.evomodel.treedatalikelihood.continuous.adapter;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalOUTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalPreparedBranchSnapshot;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalPreparedBranchSnapshotProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalPreparedBranchBasisProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTransitionCacheDiagnostics;
import dr.evomodel.continuous.ou.CanonicalPreparedBranchHandle;
import dr.evomodel.continuous.ou.CanonicalPreparedTransitionCapability;
import dr.inference.model.AbstractModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockPreparedBranchBasis;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;

/**
 * Homogeneous OU branch-transition provider for the canonical tree pathway.
 */
public final class HomogeneousCanonicalOUBranchTransitionProvider extends AbstractModel
        implements CanonicalBranchTransitionProvider,
        CanonicalOUTransitionProvider,
        CanonicalPreparedBranchBasisProvider,
        CanonicalPreparedBranchSnapshotProvider,
        CanonicalTransitionCacheDiagnostics {

    private final int dimension;
    private final MultivariateDiffusionModel diffusionModel;
    private final MatrixParameter diffusionCovariance;
    private final CanonicalOUDiffusionSnapshot diffusionSnapshot;
    private final CanonicalBranchLengthSnapshot branchLengths;
    private final OUProcessModel processModel;
    private final CanonicalTransitionCache transitionCache;

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
        this.diffusionModel = diffusionModel;
        addModel(diffusionModel);

        final MatrixParameterInterface driftMatrix = elasticModel.getStrengthOfSelectionMatrixParameter();
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

        final CanonicalPreparedTransitionCapability preparedTransition =
                processModel.getSelectionMatrixParameterization()
                        instanceof CanonicalPreparedTransitionCapability
                        ? (CanonicalPreparedTransitionCapability)
                        processModel.getSelectionMatrixParameterization()
                        : null;
        this.transitionCache = new CanonicalTransitionCache(
                dimension,
                tree.getNodeCount(),
                processModel,
                preparedTransition,
                this::getEffectiveBranchLength,
                CanonicalTransitionCacheOptions.fromSystemProperties());

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
    public void fillCanonicalTransitionForLength(final double branchLength,
                                                 final CanonicalGaussianTransition out) {
        ensureCurrentSnapshot();
        processModel.fillCanonicalTransition(branchLength, out);
    }

    @Override
    public void fillTransitionMatrix(final double branchLength, final double[][] out) {
        ensureCurrentSnapshot();
        processModel.fillTransitionMatrix(branchLength, out);
    }

    @Override
    public void fillTransitionOffset(final double branchLength, final double[] out) {
        ensureCurrentSnapshot();
        processModel.fillTransitionOffset(branchLength, out);
    }

    @Override
    public void fillTransitionCovariance(final double branchLength, final double[][] out) {
        ensureCurrentSnapshot();
        processModel.fillTransitionCovariance(branchLength, out);
    }

    @Override
    public CanonicalPreparedBranchHandle getPreparedBranchHandle(final int childNodeIndex) {
        ensureCurrentSnapshot();
        return transitionCache.getPreparedBranchHandle(childNodeIndex);
    }

    @Override
    public OrthogonalBlockPreparedBranchBasis
    getOrthogonalPreparedBranchBasis(final int childNodeIndex) {
        ensureCurrentSnapshot();
        return transitionCache.getOrthogonalPreparedBranchBasis(childNodeIndex);
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
    public void fillTraitCovariance(final double[][] out) {
        ensureCurrentSnapshot();
        diffusionSnapshot.fillCovariance(out);
    }

    public OUProcessModel getProcessModel() {
        ensureCurrentSnapshot();
        return processModel;
    }

    private void ensureCurrentSnapshot() {
        if (!dirty) {
            return;
        }
        refreshSnapshot();
        dirty = false;
        transitionCache.clear(dirtyReason);
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
        dirty = true;
        if (model == diffusionModel && dirtyReason == CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED) {
            dirtyReason = CanonicalTransitionCacheInvalidationReason.DIFFUSION_CHANGED;
        } else if (model != diffusionModel) {
            dirtyReason = CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED;
        }
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable,
                                              final int index,
                                              final Parameter.ChangeType type) {
        dirty = true;
        dirtyReason = CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED;
    }

    @Override
    public void storeState() {
        transitionCache.recordStore();
    }

    @Override
    public void restoreState() {
        transitionCache.recordRestore();
        transitionCache.clear(CanonicalTransitionCacheInvalidationReason.RESTORE_STATE);
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
}
