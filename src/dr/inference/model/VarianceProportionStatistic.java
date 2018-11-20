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

import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.MultivariateTraitDebugUtilities;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.math.matrixAlgebra.Matrix;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * A Statistic class that computes the expected proportion of the variance in the data due to diffusion on the tree
 * versus sampling error.
 *
 * @author Gabriel Hassler
 */

public class VarianceProportionStatistic extends Statistic.Abstract implements VariableListener, ModelListener {


    public static final String PARSER_NAME = "varianceProportionStatistic";
    private static final String SCALE_BY_HEIGHT = "scaleByTreeHeight";

    private TreeModel tree;
    private MultivariateDiffusionModel diffusionModel;
    private RepeatedMeasuresTraitDataModel dataModel;
    private TreeDataLikelihood treeLikelihood;
    private Parameter samplingPrecision;
    private double[] diffusionProportion;
    private boolean scaleByHeight;
    private TreeVarianceSums[] treeSums;
    private double[] diffusionVariance;
    private double[] samplingVariance;
    private int[] observedCounts;
    private ArrayList<Integer>[] perTraitMissingIndices;

    private boolean treeKnown = false;
    private boolean varianceKnown = false;


    public VarianceProportionStatistic(TreeModel tree, TreeDataLikelihood treeLikelihood,
                                       RepeatedMeasuresTraitDataModel dataModel,
                                       MultivariateDiffusionModel diffusionModel,
                                       boolean scaleByHeight) {
        this.tree = tree;
        this.treeLikelihood = treeLikelihood;
        this.diffusionModel = diffusionModel;
        this.dataModel = dataModel;
        this.samplingPrecision = dataModel.getSamplingPrecision();
        this.scaleByHeight = scaleByHeight;


//        int dim = samplingPrecision.getDimension();
        int dim = dataModel.getTraitDimension();
        this.diffusionVariance = new double[dim];
        this.samplingVariance = new double[dim];
        this.diffusionProportion = new double[dim];
        this.perTraitMissingIndices = getPerTraitMissingIndices();
        this.observedCounts = getObservedCounts();
        this.treeSums = new TreeVarianceSums[dim];
        for (int i = 0; i < dim; i++){
            treeSums[i] = new TreeVarianceSums(0, 0);
        }


        tree.addModelListener(this);
//        diffusionModel.addModelListener(this);

        diffusionModel.getPrecisionParameter().addParameterListener(this);
        samplingPrecision.addParameterListener(this);
    }

    /**
     * a class that stores the sum of the diagonal elements and all elements of a matrix
     */
    private class TreeVarianceSums {

        private double diagonalSum;
        private double totalSum;

        private TreeVarianceSums(double diagonalSum, double totalSum) {

            this.diagonalSum = diagonalSum;
            this.totalSum = totalSum;
        }

        private double getDiagonalSum() {
            return this.diagonalSum;
        }


        private double getTotalSum() {
            return this.totalSum;
        }


    }

    /**
     * @return an array with the number of taxa with observed data for each trait
     */
    private ArrayList<Integer>[] getPerTraitMissingIndices() {

        List<Integer> missingIndices = dataModel.getMissingIndices();
        int n = tree.getExternalNodeCount();
        int dim = dataModel.getTraitDimension();

        ArrayList<Integer>[] missingArrays = (ArrayList<Integer>[]) new ArrayList[dim];

        for(int i = 0; i < dim; i++){
            missingArrays[i] = new ArrayList<Integer>();
        }


        int threshold = n;
        int currentDim = 0;
        int currentTaxon = 0;

        for (int index : missingIndices) {

            currentTaxon = index / dim;
            currentDim = index - currentTaxon * dim;


            missingArrays[currentDim].add(currentTaxon);
        }

        return missingArrays;

    }

    private int[] getObservedCounts(){
        int dim = dataModel.getTraitDimension();
        int n = tree.getExternalNodeCount();
        int[] counts = new int[dim];

        for(int i = 0; i < dim; i++){

            counts[i] = n - perTraitMissingIndices[i].size();

        }

        return counts;

    }


    /**
     * recalculates the diffusionProportion statistic based on current parameters
     */
    private void updateDiffusionProportion() {
        int dim = dataModel.getTraitDimension();

        for (int i = 0; i < dim; i++) {

            double diffusionComponent = diffusionVariance[i] * (treeSums[i].getDiagonalSum() / observedCounts[i]
                    - treeSums[i].getTotalSum() / (observedCounts[i] * observedCounts[i]));

            double samplingComponent = samplingVariance[i] * (observedCounts[i] - 1) / observedCounts[i];
            diffusionProportion[i] = diffusionComponent / (diffusionComponent + samplingComponent);
        }
    }

    /**
     * recalculates the the sum of the diagonal elements and sum of all the elements of the tree variance
     * matrix statistic based on current parameters
     */
    private void updateTreeSums() {

        int dim = treeSums.length;

        double normalization = 1.0;
        if (scaleByHeight) {
            normalization = 1 / tree.getNodeHeight(tree.getRoot());
        }

        for (int i = 0; i < dim; i ++){
            double diagonalSum = MultivariateTraitDebugUtilities.getVarianceDiagonalSum(tree,
                    treeLikelihood.getBranchRateModel(), normalization, perTraitMissingIndices[i]);

            double offDiagonalSum = MultivariateTraitDebugUtilities.getVarianceOffDiagonalSum(tree,
                    treeLikelihood.getBranchRateModel(), normalization, perTraitMissingIndices[i]);

            treeSums[i].diagonalSum = diagonalSum;
            treeSums[i].totalSum = diagonalSum + offDiagonalSum;
        }


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

        int dim = dataModel.getTraitDimension();
        Matrix samplingMatrix = new Matrix(dataModel.getPrecisionMatrix().getParameterAsMatrix()).inverse();
        
        for (int i = 0; i < dim; i++) {
            samplingVariance[i] = samplingMatrix.component(i, i);
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

        if (!varianceKnown) {

            updateSamplingVariance();
            updateDiffusionVariance();
            varianceKnown = true;
            needToUpdate = true;

        }

        if (needToUpdate) {

            updateDiffusionProportion();
        }

        return diffusionProportion[dim];
    }


    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        assert  (variable == samplingPrecision || variable == diffusionModel.getPrecisionParameter());

        varianceKnown = false;
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        assert (model == tree);

        treeKnown = false;
    }

    //TODO: make its own class in evomodelxml

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
            RepeatedMeasuresTraitDataModel dataModel = (RepeatedMeasuresTraitDataModel)
                    xo.getChild(RepeatedMeasuresTraitDataModel.class);

            MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel)
                    xo.getChild(MultivariateDiffusionModel.class);

            TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
            boolean scaleByHeight = xo.getAttribute(SCALE_BY_HEIGHT, false);

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
