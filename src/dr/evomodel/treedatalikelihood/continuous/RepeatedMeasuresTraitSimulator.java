/*
 * ContinuousTraitData.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.continuous;


import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Matrix;
import dr.math.distributions.MultivariateNormalDistribution;

import java.util.ArrayList;

/**
 * @author Gabriel Hassler
 */


public class RepeatedMeasuresTraitSimulator {

    private final RepeatedMeasuresTraitDataModel dataModel;
    private final TreeDataLikelihood dataLikelihood;
    private final Tree tree;
    private final int dimTrait;
    private final int nTaxa;
    private final Parameter dataParameter;
    private final MatrixParameterInterface samplingPrecision;

    RepeatedMeasuresTraitSimulator(RepeatedMeasuresTraitDataModel dataModel, TreeDataLikelihood dataLikelihood) {

        this.dataModel = dataModel;
        this.dataLikelihood = dataLikelihood;
        this.tree = dataLikelihood.getTree();
        this.dimTrait = dataModel.dimTrait;
        this.nTaxa = tree.getExternalNodeCount();
        this.dataParameter = dataModel.getParameter();
        this.samplingPrecision = dataModel.getPrecisionMatrix();


    }

    public void simulateMissingData(double[] tipTraits) {

        //TODO: simulate diffusion on tree
        for (int i = 0; i < nTaxa; i++) {
            simulateMissingTaxonData(tipTraits, i);
        }

        dataParameter.fireParameterChangedEvent();

    }

    private void simulateMissingTaxonData(double[] tipTraits, int taxon) {

        double[] tipTrait = new double[dimTrait];
        for (int i = 0; i < dimTrait; i++) {
            tipTrait[i] = tipTraits[taxon * dimTrait + i];
        }


        double[] wholeData = new double[dimTrait];
        int nMissing = 0;
        ArrayList<Integer> observedArray = new ArrayList<Integer>();

        ArrayList<Integer> missingArray = new ArrayList<Integer>();


        for (int i = 0; i < dimTrait; i++) {

            if (dataModel.getMissingVector()[taxon * dimTrait + i]) {
                missingArray.add(i);
            } else {
                observedArray.add(i);
            }
            wholeData[i] = dataParameter.getParameterValue(taxon * dimTrait + i);

        }


        nMissing = missingArray.size();

        if (nMissing == 0) {
            return;
        }

        double[] observedData = new double[dimTrait - nMissing];
        double[] observedTip = new double[dimTrait - nMissing];
        double[] missingTip = new double[nMissing];


        for (int i = 0; i < dimTrait - nMissing; i++) {
            observedData[i] = wholeData[observedArray.get(i)];
            observedTip[i] = tipTrait[observedArray.get(i)];

        }

        for (int i = 0; i < nMissing; i++) {
            missingTip[i] = tipTrait[missingArray.get(i)];

        }


        Matrix missingPrecisionBlock = new Matrix(nMissing, nMissing);

        for (int i = 0; i < nMissing; i++) {
            for (int j = 0; j < nMissing; j++) {
                missingPrecisionBlock.set(i, j, samplingPrecision.getParameterValue(i, j));
            }

        }

        Matrix missingVarianceBlock = missingPrecisionBlock.inverse();

        double[] draw = null;

        if (nMissing == dimTrait){
            draw = MultivariateNormalDistribution.nextMultivariateNormalVariance(tipTrait,
                    missingVarianceBlock.toComponents());
        } else {

            Matrix missingObservedPrecisionBlock = new Matrix(nMissing, dimTrait - nMissing);

            for (int i = 0; i < nMissing; i++) {
                for (int j = nMissing; j < dimTrait; j++) {
                    missingObservedPrecisionBlock.set(i, j - nMissing, samplingPrecision.getParameterValue(i, j));
                }

            }

            double[] adjustedMean = computeAdjustedMean(missingVarianceBlock, missingObservedPrecisionBlock, observedData, observedTip, missingTip);


            draw = MultivariateNormalDistribution.nextMultivariateNormalVariance(adjustedMean,
                    missingVarianceBlock.toComponents());
        }


        for (int i = 0; i < nMissing; i++) {

            dataParameter.setParameterValueQuietly(taxon * dimTrait + missingArray.get(i), draw[i]);


        }

    }

    private double[] computeAdjustedMean(Matrix missingVarianceBlock, Matrix missingObservedPrecisionBlock,
                                         double[] observedData, double[] observedTip, double[] missingTip) {

        int nMissing = missingVarianceBlock.rows();
        double[] adjustedMean = new double[nMissing];
        double[] storageMis = new double[nMissing];

        for (int i = 0; i < nMissing; i++) {
            for (int j = 0; j < dimTrait - nMissing; j++) {
                storageMis[i] += missingObservedPrecisionBlock.component(i, j) * (observedData[j] - observedTip[j]);
            }
        }

        for (int i = 0; i < nMissing; i++) {

            for (int j = 0; j < nMissing; j++) {
                adjustedMean[i] += storageMis[j] * missingVarianceBlock.component(i, j);
            }

            adjustedMean[i] = missingTip[i] - adjustedMean[i];

        }


        return adjustedMean;

    }
}