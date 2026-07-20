/*
 * HomogeneousCanonicalBMBranchTransitionProvider.java
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

package dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuousNew.BMProcessModel;
import dr.evomodel.treedatalikelihood.continuousNew.ContinuousBranchTimeModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuousNew.DefaultContinuousBranchTimeModel;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalTransitionMomentProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.Arrays;

/**
 * Homogeneous Brownian-motion branch-transition provider for canonical traversal.
 *
 * <p>Unlike the OU provider, BM does not allocate or cache a transition object per
 * branch. Each request fills the caller-owned {@link CanonicalGaussianTransition}
 * from a cached diffusion snapshot and scalar branch-time arrays.
 */
public final class HomogeneousCanonicalBMBranchTransitionProvider extends AbstractModel
        implements CanonicalBranchTransitionProvider, CanonicalTransitionMomentProvider {

    private final Tree tree;
    private final BranchRateModel rateModel;
    private final BMProcessModel processModel;
    private final ContinuousBranchTimeModel branchTimeModel;
    private final int dimension;
    private final double[] effectiveBranchTime;
    private final double[] inverseEffectiveBranchTime;
    private final double[] logEffectiveBranchTime;

    private boolean branchTimesDirty = true;

    public HomogeneousCanonicalBMBranchTransitionProvider(final Tree tree,
                                                          final MultivariateDiffusionModel diffusionModel,
                                                          final BranchRateModel rateModel,
                                                          final ContinuousRateTransformation rateTransformation) {
        this(tree,
                new BMProcessModel("canonicalBmProvider.process", diffusionModel),
                rateModel,
                rateTransformation);
    }

    public HomogeneousCanonicalBMBranchTransitionProvider(final Tree tree,
                                                          final BMProcessModel processModel,
                                                          final BranchRateModel rateModel,
                                                          final ContinuousRateTransformation rateTransformation) {
        super("homogeneousCanonicalBMTransitionProvider");
        this.tree = tree;
        this.rateModel = rateModel;
        this.processModel = processModel;
        this.branchTimeModel = new DefaultContinuousBranchTimeModel(tree, rateModel, rateTransformation);
        this.dimension = processModel.getDimension();
        this.effectiveBranchTime = new double[tree.getNodeCount()];
        this.inverseEffectiveBranchTime = new double[tree.getNodeCount()];
        this.logEffectiveBranchTime = new double[tree.getNodeCount()];

        addModel(processModel);
        if (rateModel != null) {
            addModel(rateModel);
        }
        if (tree instanceof Model) {
            addModel((Model) tree);
        }
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public void fillCanonicalTransition(final int childNodeIndex, final CanonicalGaussianTransition out) {
        ensureBranchTimes();
        processModel.fillTransitionPrecision(effectiveBranchTime[childNodeIndex], out.precisionYY);
        final int dim2 = dimension * dimension;
        for (int i = 0; i < dim2; i++) {
            final double value = out.precisionYY[i];
            out.precisionXX[i] = value;
            out.precisionXY[i] = -value;
            out.precisionYX[i] = -value;
        }
        Arrays.fill(out.informationX, 0.0);
        Arrays.fill(out.informationY, 0.0);
        out.logNormalizer = processModel.getTransitionLogNormalizer(effectiveBranchTime[childNodeIndex]);
    }

    @Override
    public double getEffectiveBranchLength(final int childNodeIndex) {
        ensureBranchTimes();
        return effectiveBranchTime[childNodeIndex];
    }

    public double getInverseEffectiveBranchLength(final int childNodeIndex) {
        ensureBranchTimes();
        return inverseEffectiveBranchTime[childNodeIndex];
    }

    public double getLogEffectiveBranchLength(final int childNodeIndex) {
        ensureBranchTimes();
        return logEffectiveBranchTime[childNodeIndex];
    }

    @Override
    public void fillTraitCovarianceFlat(final double[] out) {
        processModel.fillTraitCovarianceFlat(out);
    }

    @Override
    public void fillTransitionMatrixFlat(final double branchLength, final double[] out) {
        processModel.fillTransitionMatrixFlat(out);
    }

    @Override
    public void fillTransitionOffset(final double branchLength, final double[] out) {
        processModel.fillTransitionOffset(out);
    }

    @Override
    public void fillTransitionCovarianceFlat(final double branchLength, final double[] out) {
        processModel.fillTransitionCovarianceFlat(branchLength, out);
    }

    @Override
    public void fillTransitionMomentsFlat(final double branchLength,
                                          final double[] transitionMatrixOut,
                                          final double[] transitionOffsetOut,
                                          final double[] transitionCovarianceOut) {
        processModel.fillTransitionMoments(
                branchLength,
                transitionMatrixOut,
                transitionOffsetOut,
                transitionCovarianceOut);
    }

    public BMProcessModel getProcessModel() {
        return processModel;
    }

    public void makeDirty() {
        branchTimesDirty = true;
        processModel.makeDirty();
    }

    private void ensureBranchTimes() {
        if (branchTimesDirty) {
            refreshBranchTimes();
            branchTimesDirty = false;
        }
    }

    private void refreshBranchTimes() {
        branchTimeModel.updateNormalization();
        final int nodeCount = tree.getNodeCount();
        for (int childNodeIndex = 0; childNodeIndex < nodeCount; childNodeIndex++) {
            final double time = branchTimeModel.getEffectiveBranchTime(childNodeIndex);
            effectiveBranchTime[childNodeIndex] = time;
            if (time > 0.0) {
                inverseEffectiveBranchTime[childNodeIndex] = 1.0 / time;
                logEffectiveBranchTime[childNodeIndex] = Math.log(time);
            } else {
                inverseEffectiveBranchTime[childNodeIndex] = Double.POSITIVE_INFINITY;
                logEffectiveBranchTime[childNodeIndex] = Double.NEGATIVE_INFINITY;
            }
        }
    }

    @Override
    protected void handleModelChangedEvent(final Model model, final Object object, final int index) {
        if (model == rateModel || model == tree) {
            branchTimesDirty = true;
        }
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable,
                                              final int index,
                                              final Parameter.ChangeType type) {
    }

    @Override
    public void storeState() {
    }

    @Override
    public void restoreState() {
        makeDirty();
    }

    @Override
    public void acceptState() {
    }
}
