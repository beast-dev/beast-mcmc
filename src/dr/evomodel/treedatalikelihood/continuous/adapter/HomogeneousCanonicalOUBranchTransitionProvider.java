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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalBranchTransitionProvider;
import dr.inference.model.AbstractModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.OrthogonalBlockDiagonalSelectionMatrixParameterization;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.concurrent.atomic.AtomicLong;

/**
 * Homogeneous OU branch-transition provider for the canonical tree pathway.
 */
public final class HomogeneousCanonicalOUBranchTransitionProvider extends AbstractModel
        implements CanonicalBranchTransitionProvider {

    private static final String DEBUG_CACHE_PROPERTY = "beast.debug.canonicalTransitionCache";
    private static final boolean DEBUG_CACHE = Boolean.getBoolean(DEBUG_CACHE_PROPERTY);
    private static final String PHASE_OTHER = "other";

    private final Tree tree;
    private final int dimension;
    private final BranchRateModel rateModel;
    private final ContinuousRateTransformation rateTransformation;
    private final MultivariateDiffusionModel diffusionModel;
    private final MatrixParameter diffusionCovariance;
    private final OUProcessModel processModel;
    private final CanonicalGaussianTransition[] transitionCache;
    private final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalSelection;
    private final OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis[] orthogonalPreparedCache;
    private final ThreadLocal<OrthogonalBlockDiagonalSelectionMatrixParameterization.BranchGradientWorkspace>
            orthogonalWorkspace;
    private final double[] stationaryMeanScratch;
    private final double[] cachedEffectiveBranchLength;
    private final boolean[] transitionCacheValid;

    private final DenseMatrix64F ejmlPrecision;
    private final DenseMatrix64F ejmlCovariance;
    private final CacheDiagnostics cacheDiagnostics;

    private boolean dirty = false;

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
        this.tree = tree;
        this.rateModel = rateModel;
        this.rateTransformation = rateTransformation;
        this.diffusionModel = diffusionModel;
        addModel(diffusionModel);

        final MatrixParameterInterface driftMatrix = elasticModel.getStrengthOfSelectionMatrixParameter();
        this.dimension = driftMatrix.getRowDimension();
        this.diffusionCovariance = new MatrixParameter("canonicalOuProvider.diffusion", dimension, dimension);

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

        this.transitionCache = new CanonicalGaussianTransition[tree.getNodeCount()];
        this.orthogonalSelection =
                processModel.getSelectionMatrixParameterization()
                        instanceof OrthogonalBlockDiagonalSelectionMatrixParameterization
                        ? (OrthogonalBlockDiagonalSelectionMatrixParameterization)
                        processModel.getSelectionMatrixParameterization()
                        : null;
        this.orthogonalPreparedCache =
                orthogonalSelection == null
                        ? null
                        : new OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis[tree.getNodeCount()];
        this.orthogonalWorkspace =
                orthogonalSelection == null
                        ? null
                        : ThreadLocal.withInitial(orthogonalSelection::createBranchGradientWorkspace);
        this.stationaryMeanScratch = orthogonalSelection == null ? null : new double[dimension];
        this.cachedEffectiveBranchLength = new double[tree.getNodeCount()];
        this.transitionCacheValid = new boolean[tree.getNodeCount()];

        this.ejmlPrecision = new DenseMatrix64F(dimension, dimension);
        this.ejmlCovariance = new DenseMatrix64F(dimension, dimension);
        this.cacheDiagnostics = DEBUG_CACHE ? new CacheDiagnostics() : null;
        refreshSnapshot();
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public void fillCanonicalTransition(final int childNodeIndex, final CanonicalGaussianTransition out) {
        ensureCurrentSnapshot();
        copyTransition(ensureCachedTransition(childNodeIndex), out);
    }

    public OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis
    getOrthogonalPreparedBranchBasis(final int childNodeIndex) {
        if (orthogonalSelection == null) {
            return null;
        }
        ensureCurrentSnapshot();
        ensureCachedTransition(childNodeIndex);
        return orthogonalPreparedCache[childNodeIndex];
    }

    public String pushDiagnosticPhase(final String phase) {
        if (cacheDiagnostics == null) {
            return null;
        }
        final String previous = cacheDiagnostics.phase.get();
        cacheDiagnostics.phase.set(phase);
        return previous;
    }

    public void popDiagnosticPhase(final String previous) {
        if (cacheDiagnostics == null) {
            return;
        }
        if (previous == null) {
            cacheDiagnostics.phase.remove();
        } else {
            cacheDiagnostics.phase.set(previous);
        }
    }

    private CanonicalGaussianTransition ensureCachedTransition(final int childNodeIndex) {
        recordTransitionRequest();
        final double effectiveBranchLength = getEffectiveBranchLength(childNodeIndex);
        if (!transitionCacheValid[childNodeIndex]
                || Double.doubleToLongBits(cachedEffectiveBranchLength[childNodeIndex])
                != Double.doubleToLongBits(effectiveBranchLength)) {
            recordTransitionMiss();
            CanonicalGaussianTransition cached = transitionCache[childNodeIndex];
            if (cached == null) {
                cached = new CanonicalGaussianTransition(dimension);
                transitionCache[childNodeIndex] = cached;
            }
            if (orthogonalSelection == null) {
                processModel.fillCanonicalTransition(effectiveBranchLength, cached);
            } else {
                OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis prepared =
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
            }
            cachedEffectiveBranchLength[childNodeIndex] = effectiveBranchLength;
            transitionCacheValid[childNodeIndex] = true;
        } else {
            recordTransitionHit();
        }
        return transitionCache[childNodeIndex];
    }

    @Override
    public double getEffectiveBranchLength(final int childNodeIndex) {
        final NodeRef node = tree.getNode(childNodeIndex);
        final double rawLength = tree.getBranchLength(node);
        final double normalization = rateTransformation == null ? 1.0 : rateTransformation.getNormalization();
        if (rateModel == null) {
            return rawLength * normalization;
        }
        return rawLength * rateModel.getBranchRate(tree, node) * normalization;
    }

    @Override
    public void fillTraitCovariance(final double[][] out) {
        ensureCurrentSnapshot();
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                out[i][j] = diffusionCovariance.getParameterValue(i, j);
            }
        }
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
        clearTransitionCache();
    }

    private void refreshSnapshot() {
        recordSnapshotRefresh();
        final double[][] precision = diffusionModel.getPrecisionMatrix();
        final double[] precData = ejmlPrecision.data;
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                precData[i * dimension + j] = precision[i][j];
            }
        }
        ejmlCovariance.set(ejmlPrecision);
        CommonOps.invert(ejmlCovariance);

        final double[] covData = ejmlCovariance.data;
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                diffusionCovariance.setParameterValueQuietly(i, j, covData[i * dimension + j]);
            }
        }
        diffusionCovariance.fireParameterChangedEvent();
        processModel.fireModelChanged();
    }

    private void clearTransitionCache() {
        recordCacheClear();
        for (int i = 0; i < transitionCacheValid.length; i++) {
            transitionCacheValid[i] = false;
        }
    }

    private void recordTransitionRequest() {
        if (cacheDiagnostics != null) {
            cacheDiagnostics.recordRequest(currentDiagnosticPhase());
        }
    }

    private void recordTransitionHit() {
        if (cacheDiagnostics != null) {
            cacheDiagnostics.recordHit(currentDiagnosticPhase());
        }
    }

    private void recordTransitionMiss() {
        if (cacheDiagnostics != null) {
            cacheDiagnostics.recordMiss(currentDiagnosticPhase());
        }
    }

    private void recordCacheClear() {
        if (cacheDiagnostics != null) {
            cacheDiagnostics.recordClear();
        }
    }

    private void recordSnapshotRefresh() {
        if (cacheDiagnostics != null) {
            cacheDiagnostics.recordSnapshotRefresh();
        }
    }

    private String currentDiagnosticPhase() {
        final String phase = cacheDiagnostics.phase.get();
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
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable,
                                              final int index,
                                              final Parameter.ChangeType type) {
        dirty = true;
    }

    @Override
    public void storeState() {
        if (cacheDiagnostics != null) {
            cacheDiagnostics.recordStore();
        }
    }

    @Override
    public void restoreState() {
        if (cacheDiagnostics != null) {
            cacheDiagnostics.recordRestore();
        }
        dirty = true;
    }

    @Override
    public void acceptState() {
        if (cacheDiagnostics != null) {
            cacheDiagnostics.recordAccept();
        }
    }

    public void reportTransitionCacheDiagnostics(final String label) {
        if (cacheDiagnostics != null) {
            cacheDiagnostics.report(label);
        }
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

        private void recordHit(final String phase) {
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
