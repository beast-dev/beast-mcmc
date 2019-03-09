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
import dr.evolution.tree.TreeTrait;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.RateRescalingScheme;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.MultivariateTraitDebugUtilities;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.RobustEigenDecomposition;
import dr.math.matrixAlgebra.missingData.MissingOps;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;

/**
 * A Statistic class that computes the expected proportion of the variance in the data due to diffusion on the tree
 * versus sampling error.
 *
 * @author Gabriel Hassler
 */

public class VarianceProportionStatistic extends Statistic.Abstract implements VariableListener, ModelListener {


    public static final String PARSER_NAME = "varianceProportionStatistic";
    private static final String MATRIX_RATIO = "matrixRatio";
    private static final String ELEMENTWISE = "elementWise";
    private static final String SYMMETRIC_DIVISION = "symmetricDivision";

    private final TreeModel tree;
    private final MultivariateDiffusionModel diffusionModel;
    private final RepeatedMeasuresTraitDataModel dataModel;
    private final TreeDataLikelihood treeLikelihood;
    private DenseMatrix64F diffusionProportion;
    private TreeVarianceSums treeSums;
    private Matrix diffusionVariance;
    private Matrix samplingVariance;

    private DenseMatrix64F diffusionComponent;
    private DenseMatrix64F samplingComponent;

    private boolean treeKnown = false;
    private boolean varianceKnown = false;

    private final MatrixRatios ratio;

    private final int dimTrait;

    private final static boolean useEmpiricalVariance = true;

    private final static boolean forceResample = true;


    public VarianceProportionStatistic(TreeModel tree, TreeDataLikelihood treeLikelihood,
                                       RepeatedMeasuresTraitDataModel dataModel,
                                       MultivariateDiffusionModel diffusionModel,
                                       MatrixRatios ratio) {
        this.tree = tree;
        this.treeLikelihood = treeLikelihood;
        this.diffusionModel = diffusionModel;
        this.dataModel = dataModel;
        this.dimTrait = dataModel.getTraitDimension();
        this.diffusionVariance = new Matrix(dimTrait, dimTrait);
        this.samplingVariance = new Matrix(dimTrait, dimTrait);
        this.diffusionProportion = new DenseMatrix64F(dimTrait, dimTrait);
        this.diffusionComponent = new DenseMatrix64F(dimTrait, dimTrait);
        this.samplingComponent = new DenseMatrix64F(dimTrait, dimTrait);


        this.treeSums = new TreeVarianceSums(0, 0);

        tree.addModelListener(this);
        if (useEmpiricalVariance && !forceResample) {
            dataModel.getParameter().addParameterListener(this);
        } else {
            diffusionModel.getPrecisionParameter().addParameterListener(this);
            dataModel.getPrecisionMatrix().addParameterListener(this);
        }

        this.ratio = ratio;
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

    }

    private enum MatrixRatios {
        ELEMENT_WISE {
            @Override
            void setMatrixRatio(DenseMatrix64F numeratorMatrix, DenseMatrix64F otherMatrix,
                                DenseMatrix64F destination) {
                int dim = destination.numRows;

                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {

                        double n = Math.abs(numeratorMatrix.get(i, j));
                        double d = Math.abs(otherMatrix.get(i, j));

                        if (n == 0 && d == 0) {
                            destination.set(i, j, 0);
                        } else {
                            destination.set(i, j, n / (n + d));
                        }

                    }
                }

            }
        },
        SYMMETRIC_DIVISION {
            @Override
            void setMatrixRatio(DenseMatrix64F numeratorMatrix, DenseMatrix64F otherMatrix, DenseMatrix64F destination)
                    throws IllegalDimension {


                //TODO: implement for eigendecomposition with DensMatrix64F

                System.err.println(this.name() + " not yet implemented for " + SYMMETRIC_DIVISION + ".");


//                int dim = destination.numRows;

//
//                Matrix M1 = numeratorMatrix.add(otherMatrix); //M1 = numeratorMatrix + otherMatrix
//                Matrix M2 = getMatrixSqrt(M1, true); //M2 = inv(sqrt(numeratorMatrix + otherMatrix))
//                Matrix M3 = M2.product(numeratorMatrix.product(M2));//M3 = inv(sqrt(numeratorMatrix + otherMatrix)) *
//                //                                            numeratorMatrix * inv(sqrt(numeratorMatrix + otherMatrix))
//                for (int i = 0; i < dim; i++) {
//                    for (int j = 0; j < dim; j++) {
//                        destination.set(i, j, M3.component(i, j));
//                    }
//                }

            }
        };

        abstract void setMatrixRatio(DenseMatrix64F numeratorMatrix, DenseMatrix64F otherMatrix,
                                     DenseMatrix64F destination)
                throws IllegalDimension;
    }

    private void updateDiffsionProportion() throws IllegalDimension {
        updateVarianceComponents();
        ratio.setMatrixRatio(diffusionComponent, samplingComponent, diffusionProportion);
    }


    private void updateVarianceComponents() {

        if (useEmpiricalVariance) {

            String key = REALIZED_TIP_TRAIT + "." + dataModel.getTraitName();
            TreeTrait trait = treeLikelihood.getTreeTrait(key);
            double[] tipTraits = (double[]) trait.getTrait(treeLikelihood.getTree(), null);
            double[] data = dataModel.getParameter().getParameterValues();

            int nTaxa = tree.getExternalNodeCount();

            if (forceResample) {
                //TODO: use repeatedMeasuresTraitSimulator to force sampling
            }

            computeVariance(diffusionComponent, tipTraits, nTaxa, dimTrait);
            computeVariance(samplingComponent, data, nTaxa, dimTrait);

            CommonOps.addEquals(samplingComponent, -1, diffusionComponent);


        } else {
            double n = tree.getExternalNodeCount();

            double diffusionScale = (treeSums.diagonalSum / n - treeSums.totalSum / (n * n));
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
            System.err.println("VarianceProportionStatistic not yet implemented for " +
                    "traitDataLikelihood argument useTreeLength='true'.");
        } else if (rescalingScheme != RateRescalingScheme.NONE) {
            System.err.println("VarianceProportionStatistic not yet implemented for RateRescalingShceme" +
                    rescalingScheme.getText() + ".");
        }

        treeSums.diagonalSum = diagonalSum / normalization;
        treeSums.totalSum = (diagonalSum + offDiagonalSum) / normalization;
    }

    private void updateSamplingVariance() {

        samplingVariance = dataModel.getSamplingVariance();

    }

    private void updateDiffusionVariance() {

        diffusionVariance = new Matrix(diffusionModel.getPrecisionmatrix()).inverse();

    }


    //TODO: move to difference class
    private void computeVariance(DenseMatrix64F matrix, double[] data, int numRows, int numCols) {

        double[] buffer = new double[numRows];
        DenseMatrix64F sumVec = new DenseMatrix64F(numRows, 1);
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
        CommonOps.addEquals(matrix, -1 / numRows, matrixBuffer);


    }

    //TODO: Move method below to a different class
    //TODO: implement this for DenseMatrix64F rather than Matrix
    private static Matrix getMatrixSqrt(Matrix M, Boolean invert) {
        DoubleMatrix2D S = new DenseDoubleMatrix2D(M.toComponents());
        RobustEigenDecomposition eigenDecomp = new RobustEigenDecomposition(S, 100);
        DoubleMatrix1D eigenValues = eigenDecomp.getRealEigenvalues();
        int dim = eigenValues.size();
        for (int i = 0; i < dim; i++) {
            double value = Math.sqrt(eigenValues.get(i));
            if (invert) {
                value = 1 / value;
            }
            eigenValues.set(i, value);
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


    @Override
    public int getDimension() {
        return dimTrait * dimTrait;
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

            if (useEmpiricalVariance) {

            } else {
                updateSamplingVariance();
                updateDiffusionVariance();
            }
            if (!forceResample) {
                varianceKnown = true;
            }
            needToUpdate = true;

        }

        if (needToUpdate) {

            try {
                updateDiffsionProportion();
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }

        }


        int d1 = dim / dimTrait;
        int d2 = dim - d1 * dimTrait;
        return diffusionProportion.get(d1, d2);

    }


    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        assert (variable == dataModel.getSamplingPrecision() || variable == diffusionModel.getPrecisionParameter());

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

            String ratioString = xo.getStringAttribute(MATRIX_RATIO);

            MatrixRatios ratio = null;

            if (ratioString.equalsIgnoreCase(ELEMENTWISE)) {
                ratio = MatrixRatios.ELEMENT_WISE;
            } else if (ratioString.equalsIgnoreCase(SYMMETRIC_DIVISION)) {
                ratio = MatrixRatios.SYMMETRIC_DIVISION;
            } else {
                throw new RuntimeException(PARSER_NAME + " must have attibute " + MATRIX_RATIO +
                        " with one of the following values: " + MatrixRatios.values());
            }

            return new VarianceProportionStatistic(tree, treeLikelihood, dataModel, diffusionModel,
                    ratio);
        }

        private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newStringRule(MATRIX_RATIO, false),
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
