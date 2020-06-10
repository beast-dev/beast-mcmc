/*
 * VarianceProportionStatisticEmpirical.java
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

import dr.evolution.tree.TreeTrait;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.evomodel.treedatalikelihood.preorder.ContinuousExtensionDelegate;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;

/**
 * Using Empirical variance.
 *
 * @author Gabriel Hassler
 */

public class VarianceProportionStatisticEmpirical extends AbstractVarianceProportionStatistic {

    private final ContinuousExtensionDelegate extensionDelegate;
    private final TreeTrait treeTrait;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;

    private final boolean forceResample;

    public VarianceProportionStatisticEmpirical(TreeModel tree, TreeDataLikelihood treeLikelihood,
                                                RepeatedMeasuresTraitDataModel dataModel,
                                                MultivariateDiffusionModel diffusionModel,
                                                MatrixRatios ratio,
                                                boolean forceResample) {

        super(tree, treeLikelihood, dataModel, ratio);

        this.treeTrait = treeLikelihood.getTreeTrait(REALIZED_TIP_TRAIT + "." + dataModel.getTraitName());
        this.likelihoodDelegate = (ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate();
        this.extensionDelegate = dataModel.getExtensionDelegate(likelihoodDelegate, treeTrait,
                treeLikelihood.getTree());

        this.forceResample = forceResample;

    }

    @Override
    protected void updateVarianceComponents() {

        if (forceResample) {
            likelihoodDelegate.fireModelChanged(); //Forces new sample
        }
        double[] tipTraits = (double[]) treeTrait.getTrait(treeLikelihood.getTree(), null);
        double[] data = extensionDelegate.getExtendedValues(tipTraits);
        treeLikelihood.getLogLikelihood(); //Needed to avoid assertion error

        int nTaxa = tree.getExternalNodeCount();

        computeVariance(diffusionComponent, tipTraits, nTaxa, dimTrait);
        computeVariance(samplingComponent, data, nTaxa, dimTrait);

        CommonOps.addEquals(samplingComponent, -1, diffusionComponent);
    }

    //TODO: move to difference class
    private void computeVariance(DenseMatrix64F matrix, double[] data, int numRows, int numCols) {

        double[] buffer = new double[numRows];
        DenseMatrix64F sumVec = new DenseMatrix64F(numCols, 1);
        DenseMatrix64F matrixBuffer = new DenseMatrix64F(numCols, numCols);

        Arrays.fill(matrix.getData(), 0);

        for (int i = 0; i < numRows; i++) {
            int offset = numCols * i;

            DenseMatrix64F wrapper = MissingOps.wrap(data, offset, numCols, 1, buffer);

            CommonOps.multTransB(wrapper, wrapper, matrixBuffer);
            CommonOps.addEquals(matrix, matrixBuffer);
            CommonOps.addEquals(sumVec, wrapper);

        }

        CommonOps.multTransB(sumVec, sumVec, matrixBuffer);
        CommonOps.addEquals(matrix, -1.0 / numRows, matrixBuffer);
        CommonOps.scale(1.0 / numRows, matrix);
    }

    @Override
    protected boolean needToUpdate(int dim) {
        return dim == 0;
    }
}
