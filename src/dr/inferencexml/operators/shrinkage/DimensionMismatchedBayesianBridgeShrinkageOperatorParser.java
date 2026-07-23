/*
 * DimensionMismatchedBayesianBridgeShrinkageOperatorParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inferencexml.operators.shrinkage;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.shrinkage.BayesianBridgeShrinkageOperator;
import dr.inference.operators.shrinkage.DimensionMismatchedBayesianBridgeShrinkageOperator;
import dr.math.distributions.GammaDistribution;
import dr.xml.*;

public class DimensionMismatchedBayesianBridgeShrinkageOperatorParser extends AbstractXMLObjectParser {

    public final static String BAYESIAN_BRIDGE_PARSER = "dimensionMismatchedBayesianBridgeGibbsOperator";
    public final static String MASK = "mask";
    public final static String WEIGHT = "weight";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(WEIGHT);

        BayesianBridgeStatisticsProvider bayesianBridge =
                (BayesianBridgeStatisticsProvider) xo.getChild(BayesianBridgeStatisticsProvider.class);


        GammaDistribution globalScalePrior = null;

        // This prior is actually on phi = globalScale^-exponent
        DistributionLikelihood prior = (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);
        if (prior != null) {
            if (prior.getDistribution() instanceof GammaDistribution) {
                globalScalePrior = (GammaDistribution) prior.getDistribution();
            } else {
                throw new XMLParseException("Gibbs sampler only implemented for a gamma prior on globalScale^(-exponent).");
            }
        }

        Parameter mask = (Parameter) xo.getElementFirstChild(MASK);

        BayesianBridgeShrinkageOperator operator = new BayesianBridgeShrinkageOperator(bayesianBridge, globalScalePrior, mask, Double.MIN_VALUE);

        return new DimensionMismatchedBayesianBridgeShrinkageOperator(bayesianBridge, globalScalePrior, mask, operator, weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WEIGHT),
            new ElementRule(BayesianBridgeStatisticsProvider.class),
            new ElementRule(DistributionLikelihood.class, true),
            new ElementRule(MASK, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),

            }, false),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return DimensionMismatchedBayesianBridgeShrinkageOperator.class;
    }

    @Override
    public String getParserName() {
        return BAYESIAN_BRIDGE_PARSER;
    }
}
