/*
 * RegressionMetropolizedIndicatorOperator.java
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
import dr.inferencexml.model.MaskedParameterParser;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.xml.*;

/**
 * @author Marc Suchard
 */
public class RegressionMetropolizedIndicatorOperator extends SimpleMCMCOperator {

    public static final String MH_OPERATOR = "regressionMetropolizedIndicatorOperator";

    public RegressionMetropolizedIndicatorOperator(LinearRegression linearModel, Parameter effect,
                                            Parameter indicators, MultivariateDistributionLikelihood effectPrior,
                                            Parameter mask) {

        effectOperator = new RegressionGibbsEffectOperator(linearModel, effect, indicators, effectPrior);
        this.effect = effect;
        this.indicators = indicators;
        this.mask = mask;
    }

    /**
     * @return a short descriptive message of the performance of this operator.
     */
    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return MH_OPERATOR;
    }

    public double doOperation() {

        double logHastingsRatio = 0.0;

        if (mask != null) {
            int sum = 0;
            for(int i=0; i<mask.getDimension(); i++)
                sum += mask.getParameterValue(i);
            if (sum == 0)
                throw new RuntimeException("Mask parameter has all zeros");
        }

        if (mean == null) {
            final int dim = effect.getDimension();
            mean = new double[dim];
            variance = new double[dim][dim];
            precision = new double[dim][dim];
        }

        // Compute back transition probability
        effectOperator.computeForwardDensity(mean, variance,precision);
        logHastingsRatio += MultivariateNormalDistribution.logPdf(effect.getParameterValues(),
                                mean, precision,
                                Math.log(MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(precision)),
                                1.0);

        // Do update
        int index;
        do {
            index = MathUtils.nextInt(indicators.getDimension());
        } while( (mask != null) && (mask.getParameterValue(index) == 0));

        indicators.setParameterValue(index, 1 - indicators.getParameterValue(index));        

        effectOperator.doOperation(); // Gibbs sample new coefficients given updated indicator

        // Get forward transition probability
        mean = effectOperator.getLastMean();
        variance = effectOperator.getLastVariance();
        precision = effectOperator.getLastPrecision();
        logHastingsRatio -= MultivariateNormalDistribution.logPdf(effect.getParameterValues(),
                                mean, precision,
                                Math.log(MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(precision)),
                                1.0);

        return logHastingsRatio;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MH_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            LinearRegression linearModel = (LinearRegression) xo.getChild(LinearRegression.class);
            Parameter effect = (Parameter) xo.getChild(Parameter.class);
            MultivariateDistributionLikelihood prior = (MultivariateDistributionLikelihood)
                    xo.getChild(MultivariateDistributionLikelihood.class);

            if (prior.getDistribution().getType().compareTo(MultivariateNormalDistribution.TYPE) != 0)
                throw new XMLParseException("Only a multivariate normal prior is conjugate");

            XMLObject cxo = (XMLObject) xo.getChild(GeneralizedLinearModelParser.INDICATOR);
            Parameter indicators = (Parameter) cxo.getChild(Parameter.class);

            cxo = (XMLObject) xo.getChild(MaskedParameterParser.MASKING);
            Parameter mask = null;
            if (cxo != null) {
                mask = (Parameter) cxo.getChild(Parameter.class);
                if (mask.getDimension() != indicators.getDimension())
                    throw new XMLParseException("Indicator and mask parameter must have the same dimension");
            }

            RegressionMetropolizedIndicatorOperator operator = new RegressionMetropolizedIndicatorOperator(linearModel, effect,
                    indicators, prior, mask);
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
                        new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class)
                        }),
                new ElementRule(MaskedParameterParser.MASKING,
                        new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class)
                        }, true)
        };

    };

    private Parameter mask;
    private Parameter indicators;
    private Parameter effect;
    private RegressionGibbsEffectOperator effectOperator;
    private double[] mean = null;
    private double[][] variance = null;
    private double[][] precision = null;

}
