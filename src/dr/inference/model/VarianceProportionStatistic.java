/*
 * TreeTraitLogger.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evomodel.treedatalikelihood.RateRescalingScheme;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.MultivariateTraitDebugUtilities;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.math.matrixAlgebra.Matrix;

/**
 * A Statistic class that computes the expected proportion of the variance in the data due to diffusion on the tree
 * versus sampling error.
 *
 * @author Gabriel Hassler
 */

public class VarianceProportionStatistic extends AbstractVarianceProportionStatistic implements VariableListener, ModelListener {

    private final MultivariateDiffusionModel diffusionModel;
    private TreeVarianceSums treeSums;

    private Matrix diffusionVariance;
    private Matrix samplingVariance;

    private boolean treeKnown = false;
    private boolean varianceKnown = false;

    public VarianceProportionStatistic(Tree tree, TreeDataLikelihood treeLikelihood,
                                       RepeatedMeasuresTraitDataModel dataModel,
                                       MultivariateDiffusionModel diffusionModel,
                                       MatrixRatios ratio) {

        super(tree, treeLikelihood, dataModel, ratio);

        this.diffusionModel = diffusionModel;
        this.treeSums = new TreeVarianceSums(0, 0);

        this.diffusionVariance = null;
        this.samplingVariance = null;

        if (isTreeRandom) {
            ((AbstractModel) tree).addModelListener(this);
        }
        diffusionModel.getPrecisionParameter().addParameterListener(this);
        dataModel.getExtensionPrecision().addParameterListener(this);
    }


    /**
     * a class that stores the sum of the diagonal elements and all elements of a matrix
     */
    protected class TreeVarianceSums {

        private double diagonalSum;
        private double totalSum;

        private TreeVarianceSums(double diagonalSum, double totalSum) {

            this.diagonalSum = diagonalSum;
            this.totalSum = totalSum;
        }

        public double getDiagonalSum() {
            return diagonalSum;
        }

        public double getTotalSum() {
            return totalSum;
        }

        public void setDiagonalSum(double diagonalSum) {
            this.diagonalSum = diagonalSum;
        }

        public void setTotalSum(double totalSum) {
            this.totalSum = totalSum;
        }
    }

    protected void updateVarianceComponents() {

        double n = tree.getExternalNodeCount();

        double diffusionScale = (treeSums.getDiagonalSum() / n - treeSums.getTotalSum() / (n * n));
        double samplingScale = (n - 1) / n;

        for (int i = 0; i < dimTrait; i++) {

            diffusionComponent.set(i, i, diffusionScale * diffusionVariance.component(i, i));
            samplingComponent.set(i, i, samplingScale * samplingVariance.component(i, i));

            for (int j = i + 1; j < dimTrait; j++) {

                double diffValue = diffusionScale * diffusionVariance.component(i, j);
                double sampValue = samplingScale * samplingVariance.component(i, j);

                diffusionComponent.set(i, j, diffValue);
                samplingComponent.set(i, j, sampValue);

                diffusionComponent.set(j, i, diffValue);
                samplingComponent.set(j, i, sampValue);

            }
        }
    }

    @Override
    protected boolean needToUpdate(int dim) {
        boolean needToUpdate = false;
        if (!treeKnown) {

            updateTreeSums();
            treeKnown = true;
            needToUpdate = true;

        }

        if (!varianceKnown) {

            samplingVariance = dataModel.getSamplingVariance();
            diffusionVariance = new Matrix(diffusionModel.getPrecisionmatrix()).inverse();

            varianceKnown = true;

            needToUpdate = true;

        }
        return needToUpdate;
    }

    /**
     * recalculates the the sum of the diagonal elements and sum of all the elements of the tree variance
     * matrix statistic based on current parameters
     */
    private void updateTreeSums() {
        double diagonalSum = MultivariateTraitDebugUtilities.getVarianceDiagonalSum(tree,
                treeLikelihood.getBranchRateModel(), 1.0);

        double offDiagonalSum = MultivariateTraitDebugUtilities.getVarianceOffDiagonalSum(tree,
                treeLikelihood.getBranchRateModel(), 1.0);

        RateRescalingScheme rescalingScheme = treeLikelihood.getDataLikelihoodDelegate().getRateRescalingScheme();

        double normalization = 1.0;
        if (rescalingScheme == RateRescalingScheme.TREE_HEIGHT) {
            normalization = tree.getNodeHeight(tree.getRoot());
        } else if (rescalingScheme == RateRescalingScheme.TREE_LENGTH) {
            //TODO: find function that returns tree length
            throw new RuntimeException("VarianceProportionStatistic not yet implemented for " +
                    "traitDataLikelihood argument useTreeLength='true'.");
        } else if (rescalingScheme != RateRescalingScheme.NONE) {
            throw new RuntimeException("VarianceProportionStatistic not yet implemented for RateRescalingShceme" +
                    rescalingScheme.getText() + ".");
        }

        treeSums.setDiagonalSum(diagonalSum / normalization);
        treeSums.setTotalSum((diagonalSum + offDiagonalSum) / normalization);
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        assert (variable == dataModel.getExtensionPrecision() || variable == diffusionModel.getPrecisionParameter());

        varianceKnown = false;
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        assert (model == tree);

        if (!isTreeRandom) throw new IllegalStateException("Attempting to change a fixed tree");

        treeKnown = false;
    }

    @Override
    public void modelRestored(Model model) {
        // Do nothing
    }
}
