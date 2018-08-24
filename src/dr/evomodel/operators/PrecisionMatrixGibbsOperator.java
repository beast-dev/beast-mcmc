/*
 * PrecisionMatrixGibbsOperator.java
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

package dr.evomodel.operators;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.continuous.SampledMultivariateTraitLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.distribution.MultivariateNormalDistributionModel;
import dr.inference.distribution.WishartGammalDistributionModel;
import dr.inference.model.DiagonalConstrainedMatrixView;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.distributions.WishartDistribution;
import dr.math.distributions.WishartStatistics;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.math.matrixAlgebra.Vector;
import dr.util.Attribute;
import dr.xml.*;

import java.util.List;

//import dr.math.matrixAlgebra.Matrix;

/**
 * @author Marc Suchard
 */
public class PrecisionMatrixGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String VARIANCE_OPERATOR = "precisionGibbsOperator";
    //    public static final String PRECISION_MATRIX = "precisionMatrix";
    public static final String TREE_MODEL = "treeModel";
    //    public static final String OUTCOME = "outcome";
//    public static final String MEAN = "mean";
    public static final String DISTRIBUTION = "distribution";
    public static final String PRIOR = "prior";
//    public static final String TRAIT_MODEL = "traitModel";

    private final AbstractMultivariateTraitLikelihood traitModel;
    private AbstractMultivariateTraitLikelihood debugModel = null;
    private final ConjugateWishartStatisticsProvider conjugateWishartProvider;
    private final MultivariateDistributionLikelihood multivariateLikelihood;
    private final Parameter meanParam;
    private final MatrixParameterInterface precisionParam;

    private double priorDf;
    private SymmetricMatrix priorInverseScaleMatrix;
    private final MutableTreeModel treeModel;
    private final int dim;
    private double numberObservations;
    private final String traitName;
    private final boolean isSampledTraitLikelihood;
    private double pathWeight = 1.0;

    private boolean wishartIsModel = false;
    private WishartGammalDistributionModel priorModel = null;

    public PrecisionMatrixGibbsOperator(
            MultivariateDistributionLikelihood likelihood,
            WishartStatistics priorDistribution,
            double weight) {
        super();

        // Unnecessary variables
        this.traitModel = null;
        this.treeModel = null;
        this.traitName = null;
        this.conjugateWishartProvider = null;
        this.isSampledTraitLikelihood = false;

        this.multivariateLikelihood = likelihood;
        MultivariateNormalDistributionModel density = (MultivariateNormalDistributionModel) likelihood.getDistribution();
        this.meanParam = density.getMeanParameter();
        this.precisionParam = density.getPrecisionMatrixParameter();
        this.dim = meanParam.getDimension();

        setupWishartStatistics(priorDistribution);
        if (priorDistribution instanceof WishartGammalDistributionModel) {
            wishartIsModel = true;
            priorModel = (WishartGammalDistributionModel) priorDistribution;
        }

        setWeight(weight);
    }

    private void setupWishartStatistics(WishartStatistics priorDistribution) {
        this.priorDf = priorDistribution.getDF();
        this.priorInverseScaleMatrix = null;
        double[][] scale = priorDistribution.getScaleMatrix();
        if (scale != null)
            this.priorInverseScaleMatrix =
                    (SymmetricMatrix) (new SymmetricMatrix(scale)).inverse();
    }

    @Deprecated
    public PrecisionMatrixGibbsOperator(
            MatrixParameterInterface precisionParam,
            AbstractMultivariateTraitLikelihood traitModel,
            WishartStatistics priorDistribution,
            double weight) {
        super();
        this.traitModel = traitModel;
        this.conjugateWishartProvider = null;
        this.meanParam = null;
        this.precisionParam = precisionParam;
//        this.precisionParam = traitModel.getDiffusionModel().getPrecisionParameter();

        setupWishartStatistics(priorDistribution);
        if (priorDistribution instanceof WishartGammalDistributionModel) {
            wishartIsModel = true;
            priorModel = (WishartGammalDistributionModel) priorDistribution;
        }

        setWeight(weight);

        this.treeModel = traitModel.getTreeModel();
        traitName = traitModel.getTraitName();
        dim = precisionParam.getRowDimension(); // assumed to be square

        isSampledTraitLikelihood = (traitModel instanceof SampledMultivariateTraitLikelihood);

        if (!isSampledTraitLikelihood &&
                !(traitModel instanceof ConjugateWishartStatisticsProvider)) {
            throw new RuntimeException("Only implemented for a SampledMultivariateTraitLikelihood or " +
                    "ConjugateWishartStatisticsProvider");
        }

        multivariateLikelihood = null;
    }

    public PrecisionMatrixGibbsOperator(
            ConjugateWishartStatisticsProvider wishartStatisticsProvider,
            MatrixParameterInterface extraPrecisionParam,
            WishartStatistics priorDistribution,
            double weight,
            AbstractMultivariateTraitLikelihood debugModel) {
        super();
        this.traitModel = null;
        this.debugModel = debugModel;
        this.conjugateWishartProvider = wishartStatisticsProvider;
        this.meanParam = null;
        this.precisionParam = (extraPrecisionParam != null ? extraPrecisionParam :
                conjugateWishartProvider.getPrecisionParameter());
        isSampledTraitLikelihood = false;
        this.treeModel = null;
        this.traitName = null;

        setupWishartStatistics(priorDistribution);
        if (priorDistribution instanceof WishartGammalDistributionModel) {
            wishartIsModel = true;
            priorModel = (WishartGammalDistributionModel) priorDistribution;
        }

        setWeight(weight);

//        this.treeModel = traitModel.getTreeModel();
//        traitName = traitModel.getTraitName();
        dim = precisionParam.getRowDimension(); // assumed to be square

//        isSampledTraitLikelihood = (traitModel instanceof SampledMultivariateTraitLikelihood);

//        if (!isSampledTraitLikelihood &&
//                !(traitModel instanceof ConjugateWishartStatisticsProvider)) {
//            throw new RuntimeException("Only implemented for a SampledMultivariateTraitLikelihood or " +
//                    "ConjugateWishartStatisticsProvider");
//        }

        multivariateLikelihood = null;
    }

    public void setPathParameter(double beta) {
        if (beta < 0 || beta > 1) {
            throw new IllegalArgumentException("Illegal path weight of " + beta);
        }
        pathWeight = beta;
    }

    public int getStepCount() {
        return 1;
    }

//    private void incrementScaledSquareMatrix(double[][] out, double[][] in, double scalar, int dim) {
//        for (int i = 0; i < dim; i++) {
//            for (int j = 0; j < dim; j++) {
//                out[i][j] += scalar * in[i][j];
//            }
//        }
//    }

//    private void zeroSquareMatrix(double[][] out, int dim) {
//        for (int i = 0; i < dim; i++) {
//            for (int j = 0; j < dim; j++) {
//                out[i][j] = 0.0;
//            }
//        }
//    }

    private void incrementOuterProduct(double[][] S,
                                       MultivariateDistributionLikelihood likelihood) {

        double[] mean = likelihood.getDistribution().getMean();

        numberObservations = 0;
        List<Attribute<double[]>> dataList = likelihood.getDataList();
        int count = 0;
        for (Attribute<double[]> d : dataList) {

            double[] data = d.getAttributeValue();
            for (int i = 0; i < dim; i++) {
                data[i] -= mean[i];
            }

            for (int i = 0; i < dim; i++) {  // symmetric matrix,
                for (int j = i; j < dim; j++) {
                    S[j][i] = S[i][j] += data[i] * data[j];
                }
            }
            numberObservations += 1;
        }
    }

    private void incrementOuterProduct(double[][] S,
                                       ConjugateWishartStatisticsProvider integratedLikelihood) {


        final WishartSufficientStatistics sufficientStatistics = integratedLikelihood.getWishartStatistics();
        final double[] outerProducts = sufficientStatistics.getScaleMatrix();

        final double df = sufficientStatistics.getDf();
        if (DEBUG) {
            System.err.println("OP df = " + df);
            System.err.println("OP    = " + new Vector(outerProducts));
        }
//        System.exit(-1);

        if (debugModel != null) {
            final WishartSufficientStatistics debug = ((ConjugateWishartStatisticsProvider) debugModel).getWishartStatistics();
            System.err.println(df + " ?= " + debug.getDf());
            System.err.println(new Vector(outerProducts));
            System.err.println("");
            System.err.println(new Vector(debug.getScaleMatrix()));
            System.exit(-1);
        }


//        final double df = 2;

//        final double df = integratedLikelihood.getTotalTreePrecision();

//        System.err.println("OuterProducts = \n" + new Matrix(outerProducts));
//        System.err.println("Total tree DF  = " + df);
//        System.exit(-1);

        final int dim = S.length;
        for (int i = 0; i < dim; i++) {
            System.arraycopy(outerProducts, i * dim, S[i], 0, dim);
        }
        numberObservations = df;


//        checkDiagonals(outerProducts);


    }

//    private void checkDiagonals(double[][] S) {
//        for (int i = 0; i < S.length; ++i) {
//            if (S[i][i] < 0.0) {
//                System.err.println("ERROR diag(S)\n" + new Matrix(S));
//                System.exit(-1);
//            }
//        }
//    }

    private void incrementOuterProduct(double[][] S, NodeRef node) {

        if (!treeModel.isRoot(node)) {

            NodeRef parent = treeModel.getParent(node);
            double[] parentTrait = treeModel.getMultivariateNodeTrait(parent, traitName);
            double[] childTrait = treeModel.getMultivariateNodeTrait(node, traitName);
            double time = traitModel.getRescaledBranchLengthForPrecision(node);

            if (time > 0) {

                double sqrtTime = Math.sqrt(time);

                double[] delta = new double[dim];

                for (int i = 0; i < dim; i++)
                    delta[i] = (childTrait[i] - parentTrait[i]) / sqrtTime;

                for (int i = 0; i < dim; i++) {            // symmetric matrix,
                    for (int j = i; j < dim; j++)
                        S[j][i] = S[i][j] += delta[i] * delta[j];
                }
                numberObservations += 1; // This assumes a *single* observation per tip
            }
        }
        // recurse down tree
        for (int i = 0; i < treeModel.getChildCount(node); i++)
            incrementOuterProduct(S, treeModel.getChild(node, i));
    }

    public double[][] getOperationScaleMatrixAndSetObservationCount() {

        // calculate sum-of-the-weighted-squares matrix over tree
        double[][] S = new double[dim][dim];
        SymmetricMatrix S2;
        SymmetricMatrix inverseS2 = null;
        numberObservations = 0; // Need to reset, as incrementOuterProduct can be recursive

        if (isSampledTraitLikelihood) {
            incrementOuterProduct(S, treeModel.getRoot());
        } else { // IntegratedTraitLikelihood
            if (traitModel != null) { // is a tree
                incrementOuterProduct(S, (ConjugateWishartStatisticsProvider) traitModel); // TODO deprecate usage
            } else if (conjugateWishartProvider != null) {
                incrementOuterProduct(S, conjugateWishartProvider);
            } else { // is a normal-normal-wishart model
                incrementOuterProduct(S, multivariateLikelihood);
            }
        }

        try {
            S2 = new SymmetricMatrix(S);
            if (pathWeight != 1.0) {
                S2 = (SymmetricMatrix) S2.product(pathWeight);
            }
            if (priorInverseScaleMatrix != null)
                S2 = priorInverseScaleMatrix.add(S2);
            inverseS2 = (SymmetricMatrix) S2.inverse();

//            if (S[0][0] < 0.0) {
//                 System.err.println("ERROR A");
//                 System.err.println(new Matrix(S));
//             }
//
//            if (S2.component(0, 0) < 0.0) {
//                 System.err.println("ERROR B");
//                 System.err.println(S2);
//             }
//
//            if (inverseS2.component(0, 0) < 0.0) {
//                 System.err.println("ERROR C");
//                 System.err.println("S:\n" + new Matrix(S));
//                 System.err.println("S2:\n" + S2);
//                 System.err.println(inverseS2);
//             }

        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        assert inverseS2 != null;

        return inverseS2.toComponents();
    }

    public double doOperation() {

        if (wishartIsModel) {
            setupWishartStatistics(priorModel);
        }

        final double[][] scaleMatrix = getOperationScaleMatrixAndSetObservationCount();
        final double treeDf = numberObservations;
        final double df = priorDf + treeDf * pathWeight;

//        if (scaleMatrix[0][0] < 0.0) {
//             System.err.println("ERROR");
//             System.err.println(new Matrix(scaleMatrix));
//         }

        double[][] draw = WishartDistribution.nextWishart(df, scaleMatrix);
//        int tries  = 0;
//        int limit = 100;
//        boolean success = false;
//
//        double[][] draw = null;
//
//        while (!success && tries < limit) {

        if (DEBUG) {
            System.err.println("draw = " + new Matrix(draw));
        }

//
//
//        draw = WishartDistribution.nextWishart(df, scaleMatrix);
//
//        Matrix m = new Matrix(draw);
//        try {
//            double logDet = m.logDeterminant();
//            if (Double.isNaN(logDet)) {
//                System.err.println("Bad proposal!");
//                System.err.println("df = " + df);
//
////                System.err.println(m);
////                System.exit(-1);
//            } else {
//                success = true;
//            }
//
//        } catch (IllegalDimension illegalDimension) {
//            illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            System.exit(-1);
//        }
//
//            tries++;
//        }
//
//        if (tries >= limit) {
//            System.err.println("Too many attempts!");
//            System.exit(-1);
//        }

        for (int i = 0; i < dim; i++) {
            Parameter column = precisionParam.getParameter(i);
            for (int j = 0; j < dim; j++)
                column.setParameterValueQuietly(j, draw[j][i]);
        }
        precisionParam.fireParameterChangedEvent();

        return 0;
    }

    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return VARIANCE_OPERATOR;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return VARIANCE_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            AbstractMultivariateTraitLikelihood traitModel = (AbstractMultivariateTraitLikelihood) xo.getChild(AbstractMultivariateTraitLikelihood.class);
            ConjugateWishartStatisticsProvider ws = (ConjugateWishartStatisticsProvider) xo.getChild(ConjugateWishartStatisticsProvider.class);
            if (ws == traitModel) {
                ws = null;
            }

            MultivariateDistributionLikelihood prior = null;
            MatrixParameterInterface precMatrix = null;
            MultivariateDistributionLikelihood likelihood = null;

            if (traitModel != null) {
                precMatrix = traitModel.getDiffusionModel().getPrecisionParameter();
                prior = (MultivariateDistributionLikelihood) xo.getChild(MultivariateDistributionLikelihood.class);
            }

            if (ws != null) {
                precMatrix = ws.getPrecisionParameter();
                prior = (MultivariateDistributionLikelihood) xo.getChild(MultivariateDistributionLikelihood.class);
            }

            if (traitModel == null && ws == null) { // generic likelihood and prior
                for (int i = 0; i < xo.getChildCount(); ++i) {
                    MultivariateDistributionLikelihood density = (MultivariateDistributionLikelihood) xo.getChild(i);
                    if (density.getDistribution() instanceof WishartStatistics) {
                        prior = density;
                    } else if (density.getDistribution() instanceof MultivariateNormalDistributionModel) {
                        likelihood = density;
                        precMatrix = ((MultivariateNormalDistributionModel)
                                density.getDistribution()).getPrecisionMatrixParameter();
                    }
                }

                if (prior == null || likelihood == null) {
                    throw new XMLParseException(
                            "Must provide a multivariate normal likelihood and Wishart prior in element '" +
                                    xo.getName() + "'\n"
                    );
                }
            }

            if (!(prior.getDistribution() instanceof WishartStatistics)) {
                throw new XMLParseException("Only a Wishart distribution is conjugate for Gibbs sampling");
            }

            // Make sure precMatrix is square and dim(precMatrix) = dim(parameter)
            if (precMatrix.getColumnDimension() != precMatrix.getRowDimension()) {
                throw new XMLParseException("The variance matrix is not square or of wrong dimension");
            }

            if (traitModel != null && ws == null) {

                if (precMatrix instanceof DiagonalConstrainedMatrixView) {
                    precMatrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

                    if (precMatrix == null) {
                        throw new XMLParseException("Must provide unconstrained precision matrix");
                    }
                }

                return new PrecisionMatrixGibbsOperator(
                        precMatrix,
                        traitModel, (WishartStatistics) prior.getDistribution(), weight
                );
            } else if (ws != null) {

                if (precMatrix instanceof DiagonalConstrainedMatrixView) {
                    precMatrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

                    if (precMatrix == null) {
                        throw new XMLParseException("Must provide unconstrained precision matrix");
                    }
                } else {
                    precMatrix = null;
                }
                
                return new PrecisionMatrixGibbsOperator(
                        ws, precMatrix, (WishartStatistics) prior.getDistribution(), weight, traitModel
                );

            } else {
                return new PrecisionMatrixGibbsOperator(likelihood, (WishartStatistics) prior.getDistribution(), weight);
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a multivariate normal random walk operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(AbstractMultivariateTraitLikelihood.class, true),
                new ElementRule(ConjugateWishartStatisticsProvider.class, true),
                new ElementRule(MultivariateDistributionLikelihood.class, 1, 2),
                new ElementRule(MatrixParameterInterface.class, true),
        };
    };

    private static final boolean DEBUG = false;
}
