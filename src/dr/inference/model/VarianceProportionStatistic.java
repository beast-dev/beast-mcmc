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

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.RateRescalingScheme;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.MultivariateTraitDebugUtilities;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitSimulator;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.RobustEigenDecomposition;
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
    private static final String OFF_DIAGONAL = "computeOffDiagonal";
    private static final String REMOVE_MISSING = "removeMissingBranches";

    private final TreeModel tree;
    private final MultivariateDiffusionModel diffusionModel;
    private final RepeatedMeasuresTraitDataModel dataModel;
    private final TreeDataLikelihood treeLikelihood;
    private final Parameter samplingPrecision;
    private final boolean offDiagonal;
    private final boolean removeMissing;
    private final RateRescalingScheme rescalingScheme;
    private Matrix diffusionProportion;
    private TreeVarianceSums[][] treeSums;
    private Matrix diffusionVariance;
    private Matrix samplingVariance;
    private int[][] observedCounts;
    private ArrayList<Integer>[][] perTraitMissingIndices;

    private Matrix diffusionComponent;
    private Matrix samplingComponent;

    private boolean treeKnown = false;
    private boolean varianceKnown = false;


    public VarianceProportionStatistic(TreeModel tree, TreeDataLikelihood treeLikelihood,
                                       RepeatedMeasuresTraitDataModel dataModel,
                                       MultivariateDiffusionModel diffusionModel,
                                       boolean removeMissing,
                                       boolean offDiagonal) {
        this.tree = tree;
        this.treeLikelihood = treeLikelihood;
        this.diffusionModel = diffusionModel;
        this.dataModel = dataModel;
        this.samplingPrecision = dataModel.getSamplingPrecision();
        this.rescalingScheme = treeLikelihood.getDataLikelihoodDelegate().getRateRescalingScheme();
        this.offDiagonal = offDiagonal;
        this.removeMissing = removeMissing;


        int dim = dataModel.getTraitDimension();
        this.diffusionVariance = new Matrix(dim, dim);
        this.samplingVariance = new Matrix(dim, dim);
        this.diffusionProportion = new Matrix(dim, dim);
        this.diffusionComponent = new Matrix(dim, dim);
        this.samplingComponent = new Matrix(dim, dim);


        this.perTraitMissingIndices = getPerTraitMissingIndices();
        this.observedCounts = getObservedCounts();
        this.treeSums = new TreeVarianceSums[dim][dim];
        for (int i = 0; i < dim; i++) {
            treeSums[i][i] = new TreeVarianceSums(0, 0);
            if (offDiagonal) {
                for (int j = 0; j < dim; j++) {
                    treeSums[i][j] = new TreeVarianceSums(0, 0);
                    treeSums[j][i] = treeSums[i][j];
                }
            }
        }


        tree.addModelListener(this);

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
    private ArrayList<Integer>[][] getPerTraitMissingIndices() {

        int n = tree.getExternalNodeCount();
        int dim = dataModel.getTraitDimension();

        ArrayList<Integer>[][] missingArrays = (ArrayList<Integer>[][]) new ArrayList[dim][dim];

        for (int i = 0; i < dim; i++) {
            missingArrays[i][i] = new ArrayList<Integer>();

            if (offDiagonal) {
                for (int j = i + 1; j < dim; j++) {
                    missingArrays[i][j] = new ArrayList<Integer>();
                    missingArrays[j][i] = new ArrayList<Integer>();
                }
            }
        }

        if (removeMissing) {

            boolean[] missingVector = dataModel.getMissingVector();

            for (int taxon = 0; taxon < n; taxon++) {
                ArrayList<Integer> added = new ArrayList<Integer>();

                for (int trait = 0; trait < dim; trait++) {
                    if (missingVector[taxon * dim + trait]) {
                        missingArrays[trait][trait].add(taxon);
                        added.add(trait);

                    }
                }

                for (int i : added) {
                    for (int j = 0; j < dim; j++) {
                        if (i != j) {
                            int lastDim = missingArrays[i][j].size() - 1;
                            if (lastDim == -1 || missingArrays[i][j].get(lastDim) != taxon) {
                                missingArrays[i][j].add(taxon);
                                missingArrays[j][i].add(taxon);
                            }
                        }
                    }
                }
            }
        }

        return missingArrays;

    }

    private int[][] getObservedCounts() {
        int dim = dataModel.getTraitDimension();
        int n = tree.getExternalNodeCount();
        int[][] counts = new int[dim][dim];

        for (int i = 0; i < dim; i++) {

            counts[i][i] = n - perTraitMissingIndices[i][i].size();

            if (offDiagonal) {

                for (int j = i + 1; j < dim; j++) {

                    counts[i][j] = n - perTraitMissingIndices[i][j].size();
                    counts[j][i] = counts[i][j];
                }
            }

        }

        return counts;

    }


    /**
     * recalculates the diffusionProportion statistic based on current parameters
     */
    private void updateDiffusionProportion() throws IllegalDimension {
        int dim = dataModel.getTraitDimension();


        for (int i = 0; i < dim; i++) {

            setDiffusionComponent(i, i);
            setSamplingComponent(i, i);

            if (offDiagonal) {
                for (int j = i + 1; j < dim; j++) {
                    setDiffusionComponent(i, j);
                    setSamplingComponent(i, j);
                }
            }
        }
        Matrix totalVariance = diffusionComponent.add(samplingComponent);
        Matrix totalPrecision = diffusionComponent.add(samplingComponent).inverse();
        Matrix sqrt = getMatrixSqrt(totalPrecision);


        diffusionProportion = sqrt.product(diffusionComponent.product(sqrt));
    }

    private Matrix getMatrixSqrt(Matrix M) {
        DoubleMatrix2D S = new DenseDoubleMatrix2D(M.toComponents());
        RobustEigenDecomposition eigenDecomp = new RobustEigenDecomposition(S, 100);
        DoubleMatrix1D eigenValues = eigenDecomp.getRealEigenvalues();
        int dim = eigenValues.size();
        for (int i = 0; i < dim; i++) {
            eigenValues.set(i, Math.sqrt(eigenValues.get(i)));
        }

        DoubleMatrix2D eigenVectors = eigenDecomp.getV();
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                eigenVectors.set(i, j, eigenVectors.get(i, j) * eigenValues.get(j));

            }
        }
        DoubleMatrix2D storageMatrix = new DenseDoubleMatrix2D(dim, dim);
        eigenVectors.zMult(eigenDecomp.getV(), storageMatrix, 1, 0, false, true);


        return new Matrix(storageMatrix.toArray());

    }

    private void setDiffusionComponent(int i, int j) {

        double value = diffusionVariance.component(i, j) * (treeSums[i][j].getDiagonalSum() / observedCounts[i][j]
                - treeSums[i][j].getTotalSum() / (observedCounts[i][j] * observedCounts[i][j]));

        diffusionComponent.set(i, j, value);
        diffusionComponent.set(j, i, value);

    }

    private void setSamplingComponent(int i, int j) {

        double value = samplingVariance.component(i, j) * (observedCounts[i][j] - 1) / observedCounts[i][j];

        samplingComponent.set(i, j, value);
        samplingComponent.set(j, i, value);

    }


    /**
     * recalculates the the sum of the diagonal elements and sum of all the elements of the tree variance
     * matrix statistic based on current parameters
     */
    private void updateTreeSums() {

        int dim = treeSums.length;


        for (int i = 0; i < dim; i++) {
            updateTreeSum(i, i);

            if (offDiagonal) {
                for (int j = i + 1; j < dim; j++) {
                    updateTreeSum(i, j);
                }
            }
        }


    }

    private void updateTreeSum(int i, int j) {
        double diagonalSum = MultivariateTraitDebugUtilities.getVarianceDiagonalSum(tree,
                treeLikelihood.getBranchRateModel(), 1.0, perTraitMissingIndices[i][j]);

        double offDiagonalSum = MultivariateTraitDebugUtilities.getVarianceOffDiagonalSum(tree,
                treeLikelihood.getBranchRateModel(), 1.0, perTraitMissingIndices[i][j]);

        double normalization = 1.0;
        if (rescalingScheme == RateRescalingScheme.TREE_HEIGHT) {
            normalization = tree.getNodeHeight(tree.getRoot());
        } else if (rescalingScheme == RateRescalingScheme.TREE_LENGTH) {
            //TODO: find function that returns tree length
            System.err.println("VarianceProportionStatistic not yet implemented for " +
                    "traitDataLikelihood argument useTreeLength='true'.");
        } else if (rescalingScheme != RateRescalingScheme.NONE) {
            System.err.println("VarianceProportionStatistic not yet implemented for RateRescalingShceme" +
                    rescalingScheme.getText() + ".");
        }

        treeSums[i][j].diagonalSum = diagonalSum / normalization;
        treeSums[i][j].totalSum = (diagonalSum + offDiagonalSum) / normalization;
    }

    /**
     * recalculates the diffusionVariance variable, which stores the diagonal elements of the diffusion variance matrix,
     * by inverting the current diffusion precision matrix
     */
    private void updateDiffusionVariance() {

        Matrix diffusivityMatrix = new Matrix(diffusionModel.getPrecisionmatrix()).inverse();

        if (offDiagonal) {

            diffusionVariance = diffusivityMatrix;

        } else {

            int dim = dataModel.getTraitDimension();

            for (int i = 0; i < dim; i++) {

                diffusionVariance.set(i, i, diffusivityMatrix.component(i, i));

            }
        }
    }


    /**
     * recalculates the sampling variance for each trait based on the current sampling precision
     */
    private void updateSamplingVariance() {

        Matrix samplingMatrix = new Matrix(dataModel.getPrecisionMatrix().getParameterAsMatrix()).inverse();

        if (offDiagonal) {
            samplingVariance = samplingMatrix;
        } else {

            int dim = dataModel.getTraitDimension();

            for (int i = 0; i < dim; i++) {

                samplingVariance.set(i, i, samplingMatrix.component(i, i));
            }
        }
    }


    @Override
    public int getDimension() {
        int dim = dataModel.getTraitDimension();
        if (offDiagonal) {
            return dim * dim;
        } else {
            return dim;
        }
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

            try {
                updateDiffusionProportion();
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }
        }

        if (offDiagonal) {
            int p = dataModel.getTraitDimension();
            int d1 = dim / p;
            int d2 = dim - d1 * p;
            return diffusionProportion.component(d1, d2);

        } else {
            return diffusionProportion.component(dim, dim);
        }

    }


    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        assert (variable == samplingPrecision || variable == diffusionModel.getPrecisionParameter());

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

            boolean offDiagonal = xo.getAttribute(OFF_DIAGONAL, false);
            boolean removeMissing = xo.getAttribute(REMOVE_MISSING, false);

            return new VarianceProportionStatistic(tree, treeLikelihood, dataModel, diffusionModel,
                    removeMissing, offDiagonal);
        }

        private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newStringRule(OFF_DIAGONAL, true),
                AttributeRule.newStringRule(REMOVE_MISSING, true),
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
