/*
 * TreeTraitLogger.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.inference.model;

import dr.evolution.tree.Tree;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.DiffusionProcessDelegate;
import dr.evomodel.treedatalikelihood.continuous.MultivariateTraitDebugUtilities;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Using the "population" variances, computed from the model.
 *
 * @author Gabriel Hassler
 * @author Paul Bastide
 */

public class VarianceProportionStatisticPopulation extends AbstractVarianceProportionStatistic implements VariableListener, ModelListener {

    private final MultivariateDiffusionModel diffusionModel;
    private VarianceProportionStatistic.TreeVarianceSums treeSums;

    private DenseMatrix64F diffusionVariance;
    private DenseMatrix64F samplingVariance;

    private boolean treeKnown = false;
    private boolean varianceKnown = false;

    private final DiffusionProcessDelegate diffusionProcessDelegate;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;

    private double[] treeDepths;

    public VarianceProportionStatisticPopulation(Tree tree, TreeDataLikelihood treeLikelihood,
                                                 RepeatedMeasuresTraitDataModel dataModel,
                                                 MultivariateDiffusionModel diffusionModel,
                                                 MatrixRatios ratio) {

        super(tree, treeLikelihood, dataModel, ratio);

        this.likelihoodDelegate = (ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate();

        this.diffusionModel = diffusionModel;
        this.diffusionProcessDelegate = likelihoodDelegate.getDiffusionProcessDelegate();

        this.treeDepths = new double[tree.getExternalNodeCount()];

        this.diffusionVariance = null;
        this.samplingVariance = null;

        if (isTreeRandom) {
            ((AbstractModel) tree).addModelListener(this);
        }
        this.diffusionModel.getPrecisionParameter().addParameterListener(this);
        this.dataModel.getExtensionPrecision().addParameterListener(this);
        diffusionProcessDelegate.addModelListener(this);
    }


    protected void updateVarianceComponents() {

        diffusionProcessDelegate.getMeanTipVariances(
                likelihoodDelegate.getRootProcessDelegate().getPseudoObservations(), treeDepths, diffusionVariance, diffusionComponent);

        dataModel.getMeanTipVariances(samplingVariance, samplingComponent);
    }

    @Override
    protected boolean needToUpdate(int dim) {
        boolean needToUpdate = false;
        if (!treeKnown) {

            updateTreeDepths();
            treeKnown = true;
            needToUpdate = true;

        }

        if (!varianceKnown) {

            samplingVariance = MissingOps.wrap(dataModel.getSamplingVariance().toArrayComponents(), 0, dimTrait, dimTrait);
            diffusionVariance = MissingOps.wrap(diffusionModel.getPrecisionmatrixAsVector(), 0, dimTrait, dimTrait);
            CommonOps.invert(diffusionVariance);

            varianceKnown = true;

            needToUpdate = true;

        }
        return needToUpdate;
    }

    private void updateTreeDepths() {
        double normalization = likelihoodDelegate.getRateTransformationNormalization();

        treeDepths = MultivariateTraitDebugUtilities.getTreeDepths(
                tree, treeLikelihood.getBranchRateModel(), normalization);
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        assert (variable == dataModel.getExtensionPrecision() || variable == diffusionModel.getPrecisionParameter());

        varianceKnown = false;
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        assert (model == tree || model == diffusionProcessDelegate);

        if (model == tree) {
            if (!isTreeRandom) throw new IllegalStateException("Attempting to change a fixed tree");
            treeKnown = false;
        }
        if (model == diffusionProcessDelegate) {
            varianceKnown = false;
        }
    }

    @Override
    public void modelRestored(Model model) {
        // Do nothing
    }
}
