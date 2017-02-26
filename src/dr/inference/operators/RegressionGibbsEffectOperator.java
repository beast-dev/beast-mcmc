/*
 * RegressionGibbsEffectOperator.java
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

package dr.inference.operators;

import dr.inference.distribution.LinearRegression;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Parameter;
import dr.inferencexml.distribution.GeneralizedLinearModelParser;
import dr.math.distributions.MultivariateDistribution;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.xml.*;


/**
 * @author Marc Suchard
 */
public class RegressionGibbsEffectOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String GIBBS_OPERATOR = "regressionGibbsEffectOperator";

    private LinearRegression linearModel;
    private Parameter effect;
    private Parameter indicators;
    private boolean hasNoIndicators = true;
    private MultivariateDistribution effectPrior;
    private int dim;
    private int effectNumber;
    private int N;
    private int numEffects;
    private double[][] X;

    private double[] mean = null;
    private double[][] variance = null;
    private double[][] precision = null;

    public RegressionGibbsEffectOperator(LinearRegression linearModel, Parameter effect, Parameter indicators,
                                   MultivariateDistributionLikelihood effectPrior) {
        super();
        this.linearModel = linearModel;
        this.effect = effect;
        this.indicators = indicators;
        if (indicators != null) {
            hasNoIndicators = false;
            if (indicators.getDimension() != effect.getDimension())
                throw new RuntimeException("Indicator and effect dimensions must match");
        }
        effectNumber = linearModel.getEffectNumber(effect);
        this.effectPrior = effectPrior.getDistribution();
        dim = effect.getDimension();
        N = linearModel.getDependentVariable().getDimension();
        numEffects = linearModel.getNumberOfFixedEffects();
        X = linearModel.getX(effectNumber);
    }

    public int getStepCount() {
        return 1;
    }

    public void computeForwardDensity(double[] outMean, double[][] outVariance, double[][] outPrecision) {

         double[] W = linearModel.getTransformedDependentParameter();
         double[] P = linearModel.getScale();  // outcome precision, fresh copy

         for (int k = 0; k < numEffects; k++) {
             if (k != effectNumber) {
                 double[] thisXBeta = linearModel.getXBeta(k);
                 for (int i = 0; i < N; i++)
                     W[i] -= thisXBeta[i];
             }
         }

         double[] priorBetaMean = effectPrior.getMean();
         double[][] priorBetaScale = effectPrior.getScaleMatrix();

         double[][] XtP = new double[dim][N];
         for (int j = 0; j < dim; j++) {
             if (hasNoIndicators || indicators.getParameterValue(j) == 1) {
                  for (int i = 0; i < N; i++)
                     XtP[j][i] = X[i][j] * P[i];
             } // else already filled with zeros
         }

         double[][] XtPX = new double[dim][dim];
         for (int i = 0; i < dim; i++) {
             if (hasNoIndicators || indicators.getParameterValue(i) == 1) {
                 for (int j = i; j < dim; j++) {// symmetric
                     if (hasNoIndicators || indicators.getParameterValue(j) == 1) {
                         for (int k = 0; k < N; k++)
                             XtPX[i][j] += XtP[i][k] * X[k][j];
                         XtPX[j][i] = XtPX[i][j]; // symmetric
                     }
                 }
             }
         }

         double[][] XtPX_plus_P0 = new double[dim][dim];
         for (int i = 0; i < dim; i++) {
             for (int j = i; j < dim; j++) // symmetric
                 XtPX_plus_P0[j][i] = XtPX_plus_P0[i][j] = XtPX[i][j] + priorBetaScale[i][j];
         }

         double[] XtPW = new double[dim];
         for (int i = 0; i < dim; i++) {
             for (int j = 0; j < N; j++)
                 XtPW[i] += XtP[i][j] * W[j];
         }

         double[] P0Mean0 = new double[dim];
         for (int i = 0; i < dim; i++) {
             for (int j = 0; j < dim; j++)
                 P0Mean0[i] += priorBetaScale[i][j] * priorBetaMean[j];
         }

         double[] unscaledMean = new double[dim];
         for (int i = 0; i < dim; i++)
             unscaledMean[i] = P0Mean0[i] + XtPW[i];

         double[][] variance = new SymmetricMatrix(XtPX_plus_P0).inverse().toComponents();

         for (int i = 0; i < dim; i++) {
             outMean[i] = 0.0;
             for (int j = 0; j < dim; j++) {
                 outMean[i] += variance[i][j] * unscaledMean[j];
                 outVariance[i][j] = variance[i][j];
                 outPrecision[i][j] = XtPX_plus_P0[i][j];
             }
         }
    }

    public double[] getLastMean() { return mean; }

    public double[][] getLastVariance() { return variance; }

    public double[][] getLastPrecision() { return precision; }

    public double doOperation() {

        if (mean == null)
            mean = new double[dim];

        if (variance == null)
            variance = new double[dim][dim];

        if (precision == null)
            precision = new double[dim][dim];


        computeForwardDensity(mean,variance,precision);

        double[] draw = MultivariateNormalDistribution.nextMultivariateNormalVariance(
                mean, variance);

        for (int i = 0; i < dim; i++)
            effect.setParameterValue(i, draw[i]);

        return 0;
    }

    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return GIBBS_OPERATOR;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return GIBBS_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            LinearRegression linearModel = (LinearRegression) xo.getChild(LinearRegression.class);
            Parameter effect = (Parameter) xo.getChild(Parameter.class);
            MultivariateDistributionLikelihood prior = (MultivariateDistributionLikelihood) xo.getChild(MultivariateDistributionLikelihood.class);

            if (prior.getDistribution().getType().compareTo(MultivariateNormalDistribution.TYPE) != 0)
                throw new XMLParseException("Only a multivariate normal prior is conjugate");

            XMLObject cxo = xo.getChild(GeneralizedLinearModelParser.INDICATOR);
            Parameter indicators = null;
            if (cxo != null) {
                indicators = (Parameter) cxo.getChild(Parameter.class);
            }

            RegressionGibbsEffectOperator operator = new RegressionGibbsEffectOperator(linearModel, effect, indicators, prior);
            operator.setWeight(weight);
            return operator;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a multivariate Gibbs operator on an internal node trait.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(Parameter.class),
                new ElementRule(MultivariateDistributionLikelihood.class),
                new ElementRule(LinearRegression.class),
                new ElementRule(GeneralizedLinearModelParser.INDICATOR,
                        new XMLSyntaxRule[] {
                                new ElementRule(Parameter.class)
                        },true)
        };

    };

}
