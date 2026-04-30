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
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalOUTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalPreparedBranchBasisProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTransitionCacheDiagnostics;
import dr.inference.model.AbstractModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockPreparedBranchBasis;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Homogeneous OU branch-transition provider for the canonical tree pathway.
 */
public final class HomogeneousCanonicalOUBranchTransitionProvider extends AbstractModel
        implements CanonicalBranchTransitionProvider,
        CanonicalOUTransitionProvider,
        CanonicalPreparedBranchBasisProvider,
        CanonicalTransitionCacheDiagnostics {

    private static final String DEBUG_CACHE_PROPERTY = "beast.debug.canonicalTransitionCache";
    private static final boolean DEBUG_CACHE = Boolean.getBoolean(DEBUG_CACHE_PROPERTY);

    private final Tree tree;
    private final int dimension;
    private final BranchRateModel rateModel;
    private final ContinuousRateTransformation rateTransformation;
    private final MultivariateDiffusionModel diffusionModel;
    private final MatrixParameter diffusionCovariance;
    private final OUProcessModel processModel;
    private final CanonicalTransitionCache transitionCache;

    private final DenseMatrix64F ejmlPrecision;
    private final DenseMatrix64F ejmlCovariance;

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

        final OrthogonalBlockCanonicalParameterization orthogonalSelection =
                processModel.getSelectionMatrixParameterization()
                        instanceof OrthogonalBlockCanonicalParameterization
                        ? (OrthogonalBlockCanonicalParameterization)
                        processModel.getSelectionMatrixParameterization()
                        : null;
        this.transitionCache = new CanonicalTransitionCache(
                dimension,
                tree.getNodeCount(),
                processModel,
                orthogonalSelection,
                this::getEffectiveBranchLength,
                DEBUG_CACHE);

        this.ejmlPrecision = new DenseMatrix64F(dimension, dimension);
        this.ejmlCovariance = new DenseMatrix64F(dimension, dimension);
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
    public OrthogonalBlockPreparedBranchBasis
    getOrthogonalPreparedBranchBasis(final int childNodeIndex) {
        ensureCurrentSnapshot();
        return transitionCache.getOrthogonalPreparedBranchBasis(childNodeIndex);
    }

    @Override
    public String pushDiagnosticPhase(final String phase) {
        return transitionCache.pushDiagnosticPhase(phase);
    }

    @Override
    public void popDiagnosticPhase(final String previous) {
        transitionCache.popDiagnosticPhase(previous);
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
        transitionCache.clear();
    }

    private void refreshSnapshot() {
        transitionCache.recordSnapshotRefresh();
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
        transitionCache.recordStore();
    }

    @Override
    public void restoreState() {
        transitionCache.recordRestore();
        dirty = true;
    }

    @Override
    public void acceptState() {
        transitionCache.recordAccept();
    }

    public void reportTransitionCacheDiagnostics(final String label) {
        transitionCache.report(label);
    }
}
