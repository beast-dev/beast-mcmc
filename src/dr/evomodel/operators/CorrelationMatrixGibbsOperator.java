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

import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.distribution.WishartGammalDistributionModel;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.distributions.*;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;


/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */
public class CorrelationMatrixGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private static final String CORRELATION_OPERATOR = "correlationGibbsOperator";
    public static final String TREE_MODEL = "treeModel";
    public static final String DISTRIBUTION = "distribution";
    public static final String PRIOR = "prior";

    private final ConjugateWishartStatisticsProvider conjugateWishartProvider;
    private final MatrixParameterInterface inverseCorrelation;

    private Statistics priorStatistics;
    private Statistics workingStatistics;

    private double priorDf;
    private SymmetricMatrix priorInverseScaleMatrix;

    private final int dim;
    private double numberObservations;
    private double pathWeight = 1.0;

    private boolean wishartIsModel = false;
    private WishartGammalDistributionModel priorModel = null;

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

    private void normalizeToGetInverseCorrelation(double[][] precisionMatrix) {

        double[][] covarianceMatrix = new SymmetricMatrix(precisionMatrix).inverse().toComponents();
        double[] covarianceSqrtDiagonals = new double[dim];

        for (int i = 0; i < dim; i++) {
            covarianceSqrtDiagonals[i] = Math.sqrt(covarianceMatrix[i][i]);
        }


        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                precisionMatrix[i][j] = covarianceSqrtDiagonals[i] * covarianceSqrtDiagonals[j] * precisionMatrix[i][j];
            }
        }
    }

    public CorrelationMatrixGibbsOperator(
            ConjugateWishartStatisticsProvider wishartStatisticsProvider,
            MatrixParameterInterface extraPrecisionParam,
            WishartStatistics priorDistribution,
            WishartStatistics workingDistribution,
            double weight) {
        super();

        this.conjugateWishartProvider = wishartStatisticsProvider;
        this.inverseCorrelation = (extraPrecisionParam != null ? extraPrecisionParam :
                conjugateWishartProvider.getPrecisionParameter());

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

        dim = inverseCorrelation.getRowDimension(); // assumed to be square
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

    private double[] getDiagRescaleMatrix() {

        double[] scalarMatrix = new double[dim];

        for (int i = 0; i < dim; i++) {
            double g = GammaDistribution.nextGamma((dim + 1) / 2.0,
                    1.0);
            scalarMatrix[i] = Math.sqrt(inverseCorrelation.getParameterValue(i, i) / (2.0 * g));
        }

        return scalarMatrix;
    }


    private void rescaleOuterProduct(double[] outerProduct) {

        double[] scalarMatrix = getDiagRescaleMatrix();

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                outerProduct[j * dim + i] = scalarMatrix[i] * scalarMatrix[j] * outerProduct[j * dim + i];
            }
        }
    }

    private void incrementOuterProductWithRescale(double[][] S,
                                                  ConjugateWishartStatisticsProvider integratedLikelihood) {


        final WishartSufficientStatistics sufficientStatistics = integratedLikelihood.getWishartStatistics();
        final double[] outerProducts = sufficientStatistics.getScaleMatrix();

        rescaleOuterProduct(outerProducts);

        final double df = sufficientStatistics.getDf();

        if (DEBUG) {
            System.err.println("OP df = " + df);
            System.err.println("OP    = " + new Vector(outerProducts));
        }

        final int dim = S.length;
        for (int i = 0; i < dim; i++) {
            System.arraycopy(outerProducts, i * dim, S[i], 0, dim);
        }
        numberObservations = df;

    }

    private double[][] getOperationScaleMatrixAndSetObservationCount2() {

        // calculate sum-of-the-weighted-squares matrix over tree
        double[][] S = new double[dim][dim];
        SymmetricMatrix S2;
        SymmetricMatrix inverseS2 = null;
        numberObservations = 0; // Need to reset, as incrementOuterProduct can be recursive

        incrementOuterProductWithRescale(S, conjugateWishartProvider); //todo zy: outer product

        try {
            S2 = new SymmetricMatrix(S);
            if (priorInverseScaleMatrix != null) {
                S2 = priorInverseScaleMatrix.add(S2); //todo zy: new inverse scale matrix = prior inverse scale +
                // outer product. S2 in add(S2) is the outer produt

            }
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

        normalizeToGetInverseCorrelation(draw);

        for (int i = 0; i < dim; i++) {
            Parameter column = inverseCorrelation.getParameter(i);
            for (int j = 0; j < dim; j++)
                column.setParameterValueQuietly(j, draw[j][i]);
        }

        inverseCorrelation.fireParameterChangedEvent();

        return 0;
    }

    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return CORRELATION_OPERATOR;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return CORRELATION_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            ConjugateWishartStatisticsProvider ws = (ConjugateWishartStatisticsProvider) xo.getChild
                    (ConjugateWishartStatisticsProvider.class);

            MultivariateDistributionLikelihood prior;
            MatrixParameterInterface precMatrix;

            precMatrix = ws.getPrecisionParameter();
            prior = (MultivariateDistributionLikelihood) xo.getChild(MultivariateDistributionLikelihood.class);

            if (!(prior.getDistribution() instanceof WishartStatistics)) {
                throw new XMLParseException("Only a Wishart distribution is conjugate for Gibbs sampling");
            }

            // Make sure precMatrix is square and dim(precMatrix) = dim(parameter)
            if (precMatrix.getColumnDimension() != precMatrix.getRowDimension()) {
                throw new XMLParseException("The variance matrix is not square or of wrong dimension");
            }

            return new CorrelationMatrixGibbsOperator(
                    ws, precMatrix, (WishartStatistics) prior.getDistribution(),
                    null,
                    weight
            );
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
