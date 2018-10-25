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
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.inference.model.*;
import dr.math.matrixAlgebra.Matrix;
import dr.xml.*;

import java.util.List;

/**
 * A Statistic class that computes the expected proportion of the variance in the data due to diffusion on the tree
 * versus sampling error.
 *
 * @author Gabriel Hassler
 */

public class VarianceProportionStatistic extends Statistic.Abstract implements VariableListener, ModelListener {

    public static final String PARSER_NAME = "varianceProportionStatistic";
    public static final String SCALE_BY_HEIGHT = "scaleByTreeHeight";

    private TreeModel tree;
    private MultivariateDiffusionModel diffusionModel;
    private TreeDataLikelihood treeLikelihood;
    private Parameter samplingPrecision;
    private Parameter diffusionPrecision;
    private double[] diffusionProportion;
    private boolean scaleByHeight;
    private double[] treeSums = new double[2]; // Array that stores [diagonalSum, totalSum] of tree variance matrix
    private double[] diffusionVariance;
    private double[] samplingVariance;
    private int[] observedCounts;

    private boolean treeKnown = false;
    private boolean varianceKnown = false;


    public VarianceProportionStatistic(TreeModel tree, TreeDataLikelihood treeLikelihood,
                                       RepeatedMeasuresTraitDataModel dataModel,
                                       MultivariateDiffusionModel diffusionModel,
                                       boolean scaleByHeight) {
        this.tree = tree;
        this.treeLikelihood = treeLikelihood;
        this.diffusionModel = diffusionModel;
        this.samplingPrecision = dataModel.getSamplingPrecision();
        this.diffusionPrecision = diffusionModel.getPrecisionParameter();
        this.scaleByHeight = scaleByHeight;

        this.observedCounts = getObservedCounts(dataModel);

        int dim = samplingPrecision.getDimension();
        this.diffusionVariance = new double[dim];
        this.samplingVariance = new double[dim];
        this.diffusionProportion = new double[dim];

        updateTreeSums();
        updateDiffusionVariance();
        updateSamplingVariance();
        updateDiffusionProportion();

        tree.addModelListener(this);
        samplingPrecision.addParameterListener(this);
        diffusionPrecision.addParameterListener(this);

    }

    /**
     * @return an array with the number of taxa with observed data for each trait
     */
    private int[] getObservedCounts(RepeatedMeasuresTraitDataModel dataModel){

        List<Integer> missingInds = dataModel.getMissingIndices();
        int n = tree.getExternalNodeCount();
        int dim = dataModel.getTraitDimension();
        int[] observedCounts = new int[dim];

        for (int i = 0; i < dim; i++){
            observedCounts[i] = n;
        }

        int threshold = n;
        int currentDim = 0;

        for (int index : missingInds){

            if (index >= threshold){
                threshold += n;
                currentDim += 1;
            }

            observedCounts[currentDim] -= 1;
        }

        return observedCounts;

    }

    /**
     * recalculates the diffusionProportion statistic based on current parameters
     */
    private void updateDiffusionProportion() {
        int dim = samplingPrecision.getDimension();
        double[] diffusionProportion = new double[dim];
        for (int i = 0; i < dim; i++) {
            double diffusionComponent = diffusionVariance[i] * (treeSums[0] / observedCounts[i]
                    + treeSums[1] / (observedCounts[i] * observedCounts[i]));
            double samplingComponent = samplingVariance[i] * (observedCounts[i] - 1) / observedCounts[i];
            diffusionProportion[i] = diffusionComponent / (diffusionComponent + samplingComponent);
        }
    }

    /**
     * recalculates the the sum of the diagonal elements and sum of all the elements of the tree variance
     * matrix statistic based on current parameters
     */
    private void updateTreeSums(){
        assert (treeSums.length == 2);

        double normalization = 1.0;
        if (scaleByHeight){
            normalization = 1 / tree.getNodeHeight(tree.getRoot());
        }
        //TODO: implement more efficient method for computing treeVariance
        double[][] treeVariance = MultivariateTraitDebugUtilities.getTreeVariance(tree,
                treeLikelihood.getBranchRateModel(),
                normalization, Double.POSITIVE_INFINITY);

        int n = treeVariance.length;

        double diagonalSum = 0;
        double offDiagonalSum = 0;

        for (int i = 0; i < n; i++) {
            diagonalSum = diagonalSum + treeVariance[i][i];
        }

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                offDiagonalSum = offDiagonalSum + treeVariance[i][j];
            }
        }
        offDiagonalSum = offDiagonalSum * 2;
        treeSums[0] = diagonalSum;
        treeSums[1] = diagonalSum + offDiagonalSum;
    }

    /**
     * recalculates the diffusionVariance variable, which stores the diagonal elements of the diffusion variance matrix,
     * by inverting the current diffusion precision matrix
     */
    private void updateDiffusionVariance() {

        Matrix diffusivityMatrix = new Matrix(diffusionModel.getPrecisionmatrix()).inverse();
        int dim = diffusivityMatrix.rows();

        for (int i = 0; i < dim; i++) {
            diffusionVariance[i] = diffusivityMatrix.component(i, i);
        }
    }

    /**
     * recalculates the sampling variance for each trait based on the current sampling precision
     */
    private void updateSamplingVariance() {

        int dim = samplingPrecision.getDimension();
        double[] samplingPrecisionVals = samplingPrecision.getParameterValues();

        for (int i = 0; i < dim; i++) {
            samplingVariance[i] = 1 / samplingPrecisionVals[i];
        }
    }


    @Override
    public int getDimension() {
        return diffusionProportion.length;
    }

    @Override
    public double getStatisticValue(int dim) {
        boolean needToUpdate = false;

        if (!treeKnown) {

            updateTreeSums();
            treeKnown = true;
            needToUpdate = true;

        }

        if (!varianceKnown){

            updateDiffusionVariance();
            updateSamplingVariance();
            varianceKnown = true;
            needToUpdate = true;

        }

        if (needToUpdate){

            updateDiffusionProportion();
        }

        return diffusionProportion[dim];
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        varianceKnown = false;
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        assert (model == tree);

        treeKnown = false;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
            RepeatedMeasuresTraitDataModel dataModel = (RepeatedMeasuresTraitDataModel)
                    xo.getChild(RepeatedMeasuresTraitDataModel.class);
            MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel)
                    xo.getChild(MultivariateDiffusionModel.class);

//        MatrixInverseStatistic diffusionVariance = (MatrixInverseStatistic) xo.getChild(MatrixInverseStatistic.class);
//        MatrixParameter diffusionPrecision = (MatrixParameter) xo.getChild(MatrixParameter.class);
//        Parameter samplingPrecision = (Parameter) xo.getChild(Parameter.class);

            TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
            final boolean scaleByHeight;
            if (xo.hasAttribute(SCALE_BY_HEIGHT)){
                scaleByHeight = xo.getBooleanAttribute(SCALE_BY_HEIGHT);
            } else{
                scaleByHeight = false;
            }
            return new VarianceProportionStatistic(tree, treeLikelihood, dataModel, diffusionModel, scaleByHeight);
        }

        private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newStringRule(SCALE_BY_HEIGHT, true),
                new ElementRule(TreeModel.class),
                new ElementRule(TreeDataLikelihood.class),
                new ElementRule(RepeatedMeasuresTraitDataModel.class),
                new ElementRule(MultivariateDiffusionModel.class)
        };

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return "This element returns a statistic that computes proportion of variance due to diffusion on the tree";
        }

        @Override
        public Class getReturnType() {
            return VarianceProportionStatistic.class;
        }

        @Override
        public String getParserName() {
            return PARSER_NAME;
        }
    };


    @Override
    public void modelRestored(Model model) {
        // Do nothing
    }
}
