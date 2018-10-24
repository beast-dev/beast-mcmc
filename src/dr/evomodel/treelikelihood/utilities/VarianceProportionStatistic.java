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

package dr.evomodel.treelikelihood.utilities;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import dr.math.matrixAlgebra.Matrix;

/**
 * A utility class to format traits on a tree in a LogColumn-based file.  This class takes an array of TraitProvider
 * and its package location is bound to move
 *
 * @author Marc A. Suchard
 */

public class VarianceProportionStatistic extends Statistic.Abstract {


    public static final String VARIANCE_PROPORTION_STAT = "varianceProportionStatistic";

    private Tree tree;
    private RepeatedMeasuresTraitDataModel dataModel;
    private MultivariateDiffusionModel diffusionModel;
    private TreeDataLikelihood treeLikelihood;
    private double[] diffusionProportion;
    private boolean proportionKnown = false;


    public VarianceProportionStatistic(Tree tree, TreeDataLikelihood treeLikelihood, RepeatedMeasuresTraitDataModel dataModel, MultivariateDiffusionModel diffusionModel) {
        this.tree = tree;
        this.dataModel = dataModel;
        this.treeLikelihood = treeLikelihood;
        this.diffusionModel = diffusionModel;
        this.diffusionProportion = getDiffusionProportion();
//        diffusionProportion.addParameterListener(this);

    }


    private double[] getDiffusionProportion(){
        double[] diffusionVariance = getDiffusionVariance();
        double[] sampleVariance = getSampleVariance();
        int dim = dataModel.getTraitDimension();
        double[] diffusionProportion = new double[dim];
        for (int i = 0; i <= dim - 1; i++){
            diffusionProportion[i] = diffusionVariance[i] / (sampleVariance[i] + diffusionVariance[i]);
        }
        return diffusionProportion;
    }

    private double[] getDiffusionVariance() {
        double[][] treeVariance = MultivariateTraitDebugUtilities.getTreeVariance(tree, treeLikelihood.getBranchRateModel(),
                1.0, Double.POSITIVE_INFINITY);
        Matrix diffusivityMatrix = new Matrix(diffusionModel.getPrecisionmatrix()).inverse();
        int n = treeVariance.length;
        int dim = diffusivityMatrix.rows();
        double diagonalSum = 0;
        double offDiagonalSum= 0;
        for (int i = 0; i <= n - 1; i++){
            diagonalSum = diagonalSum + treeVariance[i][i];
        }
        for (int i = 0; i <= n - 2; i++){
            for (int j = i + 1; j <= n - 1; j++)
                offDiagonalSum = offDiagonalSum + treeVariance[i][j];
        }
        offDiagonalSum = offDiagonalSum * 2;
        double diffusionScalar =  diagonalSum / n - (diagonalSum + offDiagonalSum) / (n * n);

        double[] diffusionVariance = new double[dim];
        for (int i = 0; i <= dim - 1; i++){
            diffusionVariance[i] = diffusionScalar * diffusivityMatrix.component(i, i);
        }
        return diffusionVariance;
    }
    
    private double[] getSampleVariance(){
        int n = tree.getExternalNodeCount();
        int dim = dataModel.getTraitDimension();
        double[] samplingPrecision = dataModel.getSamplingPrecision().getParameterValues();
        double[] sampleVariance = new double[dim];
        for (int i = 0; i <= dim - 1; i++){
            sampleVariance[i] = (n - 1) / (n * samplingPrecision[i]);
        }
        return sampleVariance;
    }

//    private double getSamplingVariance(double[][] treeVariance){
//
//    }


    @Override
    public int getDimension() {
        return diffusionProportion.length;
    }

    @Override
    public double getStatisticValue(int dim) {
        if (!proportionKnown){
            diffusionProportion = getDiffusionProportion();
//            proportionKnown = true;
        }
        return diffusionProportion[dim];
    }

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        proportionKnown = false;
    }
}
