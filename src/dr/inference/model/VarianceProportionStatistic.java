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
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.inference.model.*;
import dr.math.matrixAlgebra.Matrix;
import dr.xml.*;

/**
 * A Statistic class that computes the expected proportion of the variance in the data due to diffusion on the tree
 * versus sampling error.
 *
 * @author Gabriel Hassler
 */

public class VarianceProportionStatistic extends Statistic.Abstract implements VariableListener {

    public static final String PARSER_NAME = "varianceProportionStatistic";
    public static final String SCALE_BY_HEIGHT = "scaleByTreeHeight";

    private Tree tree;
    private MultivariateDiffusionModel diffusionModel;
    private TreeDataLikelihood treeLikelihood;
    private Parameter samplingPrecision;
    private Parameter diffusionPrecision;
    private double[] diffusionProportion;
    private boolean scaleByHeight;


    public VarianceProportionStatistic(Tree tree, TreeDataLikelihood treeLikelihood,
                                       RepeatedMeasuresTraitDataModel dataModel,
                                       MultivariateDiffusionModel diffusionModel,
                                       boolean scaleByHeight) {
        this.tree = tree;
        this.treeLikelihood = treeLikelihood;
        this.diffusionModel = diffusionModel;
        this.samplingPrecision = dataModel.getSamplingPrecision();
        this.diffusionPrecision = diffusionModel.getPrecisionParameter();
        this.diffusionProportion = getDiffusionProportion(scaleByHeight);
        this.scaleByHeight = scaleByHeight;
        samplingPrecision.addParameterListener(this);
        diffusionPrecision.addParameterListener(this);
        //TODO: find parameters for tree branch lengths and add them to parameter listener
        //TODO: implement proportionKnown boolean variable that changes to false when one of the parameters changes

    }


    private double[] getDiffusionProportion(boolean scaleByHeight) {
        double[] diffusionVariance = getDiffusionVariance(scaleByHeight);
        double[] sampleVariance = getSampleVariance();
        int dim = samplingPrecision.getDimension();
        double[] diffusionProportion = new double[dim];
        for (int i = 0; i < dim; i++) {
            diffusionProportion[i] = diffusionVariance[i] / (sampleVariance[i] + diffusionVariance[i]);
        }
        return diffusionProportion;
    }

    private double[] getDiffusionVariance(boolean scaleByHeight) {
        double normalization = 1.0;
        if (scaleByHeight){
            normalization = 1 / tree.getNodeHeight(tree.getRoot());
        }
        //TODO: implement more efficient method for computing treeVariance
        double[][] treeVariance = MultivariateTraitDebugUtilities.getTreeVariance(tree,
                treeLikelihood.getBranchRateModel(),
                normalization, Double.POSITIVE_INFINITY);
        Matrix diffusivityMatrix = new Matrix(diffusionModel.getPrecisionmatrix()).inverse();
        int n = treeVariance.length;
        int dim = diffusivityMatrix.rows();
        double diagonalSum = 0;
        double offDiagonalSum = 0;
        for (int i = 0; i <= n - 1; i++) {
            diagonalSum = diagonalSum + treeVariance[i][i];
        }
        for (int i = 0; i <= n - 2; i++) {
            for (int j = i + 1; j <= n - 1; j++)
                offDiagonalSum = offDiagonalSum + treeVariance[i][j];
        }
        offDiagonalSum = offDiagonalSum * 2;
        double diffusionScalar = diagonalSum / n - (diagonalSum + offDiagonalSum) / (n * n);

        double[] diffusionVariance = new double[dim];
        for (int i = 0; i <= dim - 1; i++) {
            diffusionVariance[i] = diffusionScalar * diffusivityMatrix.component(i, i);
        }
        return diffusionVariance;
    }

    private double[] getSampleVariance() {
        int n = tree.getExternalNodeCount();
        int dim = samplingPrecision.getDimension();
        double[] samplingPrecisionVals = samplingPrecision.getParameterValues();
        double[] sampleVariance = new double[dim];
        for (int i = 0; i <= dim - 1; i++) {
            sampleVariance[i] = (n - 1) / (n * samplingPrecisionVals[i]);
        }
        return sampleVariance;
    }


    @Override
    public int getDimension() {
        return diffusionProportion.length;
    }

    @Override
    public double getStatisticValue(int dim) {
        if (dim == 0) {
            diffusionProportion = getDiffusionProportion(scaleByHeight);
        }
        return diffusionProportion[dim];
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {

    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Tree tree = (Tree) xo.getChild(Tree.class);
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
                new ElementRule(Tree.class),
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
}
