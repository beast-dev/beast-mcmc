/*
 * PrecisionMatrixByCorrelationGibbsOperator.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.math.distributions.InverseGammaDistribution;
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
import sun.jvm.hotspot.debugger.cdbg.Sym;

import java.util.List;

/**
 * @author Marc Suchard
 */
public class PrecisionMatrixByCorrelationGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private static final String VARIANCE_OPERATOR = "precisionGibbsOperator1";
    public static final String TREE_MODEL = "treeModel";
    public static final String DISTRIBUTION = "distribution";
    public static final String PRIOR = "prior";
    private static final String WORKING = "workingDistribution";

    private final AbstractMultivariateTraitLikelihood traitModel;
    private AbstractMultivariateTraitLikelihood debugModel = null;
    private final ConjugateWishartStatisticsProvider conjugateWishartProvider;
    private final MultivariateDistributionLikelihood multivariateLikelihood;
    private final Parameter meanParam;
    private final MatrixParameterInterface precisionParam;

    private Statistics priorStatistics;
    private Statistics workingStatistics;

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

    public PrecisionMatrixByCorrelationGibbsOperator(
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

        setupWishartStatistics(priorDistribution); // TODO Deprecate
        priorStatistics = setupStatistics(priorDistribution);

        if (priorDistribution instanceof WishartGammalDistributionModel) {
            wishartIsModel = true;
            priorModel = (WishartGammalDistributionModel) priorDistribution;
        }

        setWeight(weight);
    }

    private class Statistics {
        final double degreesOfFreedom;
        final double[][] rateMatrix;

        Statistics(double degreesOfFreedom, double[][] rateMatrix) {
            this.degreesOfFreedom = degreesOfFreedom;
            this.rateMatrix = rateMatrix;
        }
    }

    private Statistics setupStatistics(WishartStatistics distribution) {

        double[][] scale = distribution.getScaleMatrix();
        double[][] rate = null;

        if (scale != null) {
            rate = (new SymmetricMatrix(scale)).inverse().toComponents();
        }

        return new Statistics(distribution.getDF(), rate);
    }

    private void setupWishartStatistics(WishartStatistics priorDistribution) {
        this.priorDf = priorDistribution.getDF();
        this.priorInverseScaleMatrix = null;
        double[][] scale = priorDistribution.getScaleMatrix();
        if (scale != null)
            this.priorInverseScaleMatrix =
                    (SymmetricMatrix) (new SymmetricMatrix(scale)).inverse();
    }

    private void renormalizeGetInvCorr(double[][] draw) {

//        SymmetricMatrix(draw).
        double[][] cov = new SymmetricMatrix(draw).inverse().toComponents();
        double[] scaleD = new double[dim];

        for (int i = 0; i < dim; i++) {
            scaleD[i] = Math.sqrt(cov[i][i]);
        }


        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                draw[i][j] = scaleD[i] * scaleD[j] * draw[i][j];
            }
        }
    }

    @Deprecated
    public PrecisionMatrixByCorrelationGibbsOperator(
            MatrixParameterInterface precisionParam,
            AbstractMultivariateTraitLikelihood traitModel,
            WishartStatistics priorDistribution,
            double weight) {
        super();
        this.traitModel = traitModel;
        this.conjugateWishartProvider = null;
        this.meanParam = null;
        this.precisionParam = precisionParam;

        setupWishartStatistics(priorDistribution); // TODO Deprecate
        priorStatistics = setupStatistics(priorDistribution);

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

    public PrecisionMatrixByCorrelationGibbsOperator(
            ConjugateWishartStatisticsProvider wishartStatisticsProvider,
            MatrixParameterInterface extraPrecisionParam,
            WishartStatistics priorDistribution,
            WishartStatistics workingDistribution,
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

        setupWishartStatistics(priorDistribution); // TODO Deprecate
        priorStatistics = setupStatistics(priorDistribution);

        if (priorDistribution instanceof WishartGammalDistributionModel) {
            wishartIsModel = true;
            priorModel = (WishartGammalDistributionModel) priorDistribution;
        }

        if (workingDistribution != null) {
            workingStatistics = setupStatistics(workingDistribution);
        }

        setWeight(weight);

        dim = precisionParam.getRowDimension(); // assumed to be square



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

    private double[] getscalarMatrix() {

        double[] scalarMatrix = new double[dim];
        for (int i = 0; i<dim;i++){
            scalarMatrix[i] = Math.sqrt(InverseGammaDistribution.nextInverseGamma((dim + 1) / 2.0,
                    precisionParam.getParameterValue(i, i) / 2.0));
        }



        return scalarMatrix;
    }



    private void rescaleOuterProduct(double[] outerProduct) {

        double[] scalarMatrix = getscalarMatrix();

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                //outerProduct[i * dim + j] = getscalar(i) * getscalar(j) * outerProduct[i * dim + j];
                System.err.println(outerProduct[j * dim + i]);
                System.out.println("scalar 1 = " + scalarMatrix[i] + "scalar 2 = " + scalarMatrix[j]);
                outerProduct[j * dim + i] = scalarMatrix[i] * scalarMatrix[j] * outerProduct[j * dim + i];
                System.err.println(outerProduct[j * dim + i]);

            }
        }
    }

    private void incrementOuterProduct2(double[][] S,
                                       ConjugateWishartStatisticsProvider integratedLikelihood) {


        final WishartSufficientStatistics sufficientStatistics = integratedLikelihood.getWishartStatistics();
        final double[] outerProducts = sufficientStatistics.getScaleMatrix();
        if (DEBUG) {

            System.err.println("OP used to be symme   = " + new Vector(outerProducts));
        }

        rescaleOuterProduct(outerProducts);

        if (DEBUG) {

            System.err.println("OP not symme   = " + new Vector(outerProducts));
        }


        final double df = sufficientStatistics.getDf();

        if (DEBUG) {
            System.err.println("OP df = " + df);
            System.err.println("OP    = " + new Vector(outerProducts));
        }

        if (debugModel != null) {
            final WishartSufficientStatistics debug = ((ConjugateWishartStatisticsProvider) debugModel).getWishartStatistics();
            System.err.println(df + " ?= " + debug.getDf());
            System.err.println(new Vector(outerProducts));
            System.err.println("");
            System.err.println(new Vector(debug.getScaleMatrix()));
            System.exit(-1);
        }

        final int dim = S.length;
        for (int i = 0; i < dim; i++) {
            System.arraycopy(outerProducts, i * dim, S[i], 0, dim);
        }
        numberObservations = df;


//        checkDiagonals(outerProducts);


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

        if (debugModel != null) {
            final WishartSufficientStatistics debug = ((ConjugateWishartStatisticsProvider) debugModel).getWishartStatistics();
            System.err.println(df + " ?= " + debug.getDf());
            System.err.println(new Vector(outerProducts));
            System.err.println("");
            System.err.println(new Vector(debug.getScaleMatrix()));
            System.exit(-1);
        }

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

    private double[][] getOperationScaleMatrixAndSetObservationCount2() {

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
                incrementOuterProduct2(S, conjugateWishartProvider); //todo zy: outer product
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
                S2 = priorInverseScaleMatrix.add(S2); //todo zy: new inverse scale matrix = prior inverse scale + outer product. S2 in add(S2) is the outer produt
            inverseS2 = (SymmetricMatrix) S2.inverse();

        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        assert inverseS2 != null;

        return inverseS2.toComponents();
    }

    private double[][] getOperationScaleMatrixAndSetObservationCount() {

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
                incrementOuterProduct(S, conjugateWishartProvider); //todo zy: outer product
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
                S2 = priorInverseScaleMatrix.add(S2); //todo zy: new inverse scale matrix = prior inverse scale + outer product
            inverseS2 = (SymmetricMatrix) S2.inverse();

        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        assert inverseS2 != null;

        return inverseS2.toComponents();
    }

    public double doOperation() {

        if (wishartIsModel) {
            setupWishartStatistics(priorModel); // TODO Deprecate
            priorStatistics = setupStatistics(priorModel);
        }

        final double[][] scaleMatrix = getOperationScaleMatrixAndSetObservationCount2();
        final double treeDf = numberObservations;

        final double df = priorDf + treeDf * pathWeight;

        double[][] draw = WishartDistribution.nextWishart(df, scaleMatrix);

        if (DEBUG) {
            System.err.println("scale matrix for draw = " + new Matrix(scaleMatrix));
        }

//        if (DEBUG) {
//            System.err.println("draw1 = " + new Matrix(draw));
//        }
        //todo : before setting R, renomailzie SIGMA TO R
        renormalizeGetInvCorr(draw);


            System.err.println("draw2 = " + new Matrix(draw));
            System.err.println("corr = " + new Matrix(draw).inverse());



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

                return new PrecisionMatrixByCorrelationGibbsOperator(
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

                WishartStatistics workingDistribution = null;
                if (xo.hasChildNamed(WORKING)) {
                    workingDistribution = (WishartStatistics) xo.getElementFirstChild(WORKING);
                }

                return new PrecisionMatrixByCorrelationGibbsOperator(
                        ws, precMatrix, (WishartStatistics) prior.getDistribution(),
                        workingDistribution,
                        weight, traitModel
                );

            } else {
                return new PrecisionMatrixByCorrelationGibbsOperator(likelihood, (WishartStatistics) prior.getDistribution(), weight);
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
