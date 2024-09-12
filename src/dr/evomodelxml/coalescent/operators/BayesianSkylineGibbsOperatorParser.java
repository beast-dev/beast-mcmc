/*
 * BayesianSkylineGibbsOperatorParser.java
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

package dr.evomodelxml.coalescent.operators;

import dr.evomodel.coalescent.BayesianSkylineLikelihood;
import dr.evomodel.coalescent.operators.BayesianSkylineGibbsOperator;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class BayesianSkylineGibbsOperatorParser extends AbstractXMLObjectParser {
    public static final String BAYESIAN_SKYLINE_GIBBS_OPERATOR = "generalizedSkylineGibbsOperator";
    public static final String POPULATION_SIZES = "populationSizes";
    public static final String GROUP_SIZES = "groupSizes";
    public static final String LOWER = "lower";
    public static final String UPPER = "upper";

    public static final String JEFFREYS = "Jeffreys";
    public static final String EXPONENTIALMARKOV = "exponentialMarkov";
    public static final String SHAPE = "shape";
    public static final String REVERSE = "reverse";
    public static final String ITERATIONS = "iterations";
    public static final String TYPE = "type";

    public static final String STEPWISE = "stepwise";
    public static final String LINEAR = "linear";
    public static final String EXPONENTIAL = "exponential";

    public String getParserName() {
        return BAYESIAN_SKYLINE_GIBBS_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        final double lowerBound = xo.getAttribute(LOWER, 0.0);
        final double upperBound = xo.getAttribute(UPPER, Double.MAX_VALUE);
        final boolean jeffreysPrior = xo.getAttribute(JEFFREYS, true);

        boolean exponentialMarkovPrior = xo.getAttribute(EXPONENTIALMARKOV, false);
        double shape = xo.getAttribute(SHAPE, 1.0);
        boolean reverse = xo.getAttribute(REVERSE, false);
        int iterations = xo.getAttribute(ITERATIONS, 1);

        BayesianSkylineLikelihood bayesianSkylineLikelihood = (BayesianSkylineLikelihood) xo
                .getChild(BayesianSkylineLikelihood.class);

        // This is the parameter on which this operator acts
        Parameter paramPops = (Parameter) xo.getChild(Parameter.class);

        Parameter paramGroups = bayesianSkylineLikelihood.getGroupSizeParameter();

        final int type = bayesianSkylineLikelihood.getType();

        if (type != BayesianSkylineLikelihood.STEPWISE_TYPE) {
            throw new XMLParseException(
                    "Need stepwise control points (set 'linear=\"false\"' in skyline Gibbs operator)");
        }

        return new BayesianSkylineGibbsOperator(bayesianSkylineLikelihood,
                paramPops, paramGroups, type, weight, lowerBound,
                upperBound, jeffreysPrior, exponentialMarkovPrior,
                shape, reverse, iterations);

    }

    // ************************************************************************
    // AbstractXMLObjectParser implementation
    // ************************************************************************

    public String getParserDescription() {
        return "This element returns a Gibbs operator for the joint distribution of population sizes.";
    }

    public Class getReturnType() {
        return BayesianSkylineGibbsOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newBooleanRule(LINEAR, true),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(LOWER),
            AttributeRule.newDoubleRule(UPPER),
            AttributeRule.newBooleanRule(JEFFREYS, true),
            AttributeRule.newBooleanRule(REVERSE, true),
            AttributeRule.newBooleanRule(EXPONENTIALMARKOV, true),
            AttributeRule.newDoubleRule(SHAPE),
            new ElementRule(BayesianSkylineLikelihood.class),
            new ElementRule(Parameter.class)
    };

}
