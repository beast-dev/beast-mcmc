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
import dr.inference.timeseries.gaussian.OUProcessModel;
import dr.inference.timeseries.representation.CanonicalGaussianTransition;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Homogeneous OU branch-transition provider for the canonical tree pathway.
 *
 * <p>This mirrors {@link HomogeneousOUBranchKernelProvider}, but exposes the OU
 * branch factor directly in canonical form through
 * {@link OUProcessModel#fillCanonicalTransition(double, CanonicalGaussianTransition)}.
 */
public final class HomogeneousCanonicalOUBranchTransitionProvider extends AbstractModel
        implements CanonicalBranchTransitionProvider {

    private final Tree tree;
    private final int dimension;
    private final BranchRateModel rateModel;
    private final ContinuousRateTransformation rateTransformation;
    private final MultivariateDiffusionModel diffusionModel;
    private final MatrixParameter diffusionCovariance;
    private final OUProcessModel processModel;
    private final CanonicalGaussianTransition[] transitionCache;
    private final double[] cachedEffectiveBranchLength;
    private final boolean[] transitionCacheValid;

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

        this.transitionCache = new CanonicalGaussianTransition[tree.getNodeCount()];
        this.cachedEffectiveBranchLength = new double[tree.getNodeCount()];
        this.transitionCacheValid = new boolean[tree.getNodeCount()];

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

        final double effectiveBranchLength = getEffectiveBranchLength(childNodeIndex);
        if (!transitionCacheValid[childNodeIndex]
                || Double.doubleToLongBits(cachedEffectiveBranchLength[childNodeIndex])
                != Double.doubleToLongBits(effectiveBranchLength)) {
            CanonicalGaussianTransition cached = transitionCache[childNodeIndex];
            if (cached == null) {
                cached = new CanonicalGaussianTransition(dimension);
                transitionCache[childNodeIndex] = cached;
            }
            processModel.fillCanonicalTransition(effectiveBranchLength, cached);
            cachedEffectiveBranchLength[childNodeIndex] = effectiveBranchLength;
            transitionCacheValid[childNodeIndex] = true;
        }
        copyTransition(transitionCache[childNodeIndex], out);
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
        for (int i = 0; i < transitionCacheValid.length; i++) {
            transitionCacheValid[i] = false;
        }
    }

    private static void copyTransition(final CanonicalGaussianTransition source,
                                       final CanonicalGaussianTransition target) {
        final int dimension = source.getDimension();
        for (int i = 0; i < dimension; i++) {
            target.informationX[i] = source.informationX[i];
            target.informationY[i] = source.informationY[i];
            for (int j = 0; j < dimension; j++) {
                target.precisionXX[i][j] = source.precisionXX[i][j];
                target.precisionXY[i][j] = source.precisionXY[i][j];
                target.precisionYX[i][j] = source.precisionYX[i][j];
                target.precisionYY[i][j] = source.precisionYY[i][j];
            }
        }
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
    public void storeState() { }

    @Override
    public void restoreState() {
        dirty = true;
    }

    @Override
    public void acceptState() { }
}
