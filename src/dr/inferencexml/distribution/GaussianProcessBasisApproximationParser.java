/*
 * GaussianMarkovRandomFieldParser.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.distribution;

import dr.inference.distribution.RandomField;
import dr.inference.model.Parameter;
import dr.math.distributions.GaussianProcessBasisApproximation;
import dr.xml.*;


public class GaussianProcessBasisApproximationParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "gaussianProcessBasisApproximation";
    private static final String DIM = "dim";
    private static final String KNOTS = "knots";
    private static final String DEGREE = "degree";
    private static final String TIMES = "times";
    private static final String BOUNDARY = "boundary";
    private static final String MEAN = "mean";
    private static final String MARGINALVARIANCE = "marginalVariance";
    private static final String LENGTHSCALE = "lengthScale";
    private static final String COEFFICIENT = "coefficient";
    private static final String PRECISION = "precision";

    public String getParserName() { return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int knots = xo.getIntegerAttribute(KNOTS);

        int dim = xo.getIntegerAttribute(DIM);

        double degree = xo.getDoubleAttribute(DEGREE);

        double[] times = xo.getDoubleArrayAttribute(TIMES);

        double boundary = xo.getDoubleAttribute(BOUNDARY);

        Parameter mean = (Parameter) xo.getElementFirstChild(MEAN);

        Parameter marginalVariance = (Parameter) xo.getElementFirstChild(MARGINALVARIANCE);

        if (marginalVariance.getParameterValue(0) <= 0.0) {
            throw new XMLParseException("Marginal variance must be > 0.0");
        }

        Parameter lengthScale = (Parameter) xo.getElementFirstChild(LENGTHSCALE);

        if (lengthScale.getParameterValue(0) <= 0.0) {
            throw new XMLParseException("Length scale must be > 0.0");
        }


        Parameter coefficient = (Parameter) xo.getElementFirstChild(COEFFICIENT);


        Parameter precision = (Parameter) xo.getElementFirstChild(PRECISION);

        if (precision.getParameterValue(0) <= 0.0) {
            throw new XMLParseException("Noise must be > 0.0");
        }


//        RandomField.WeightProvider weights = parseWeightProvider(xo, dim);


        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        return new GaussianProcessBasisApproximation(id, dim, knots, degree, times, boundary, mean, marginalVariance,
                lengthScale, coefficient, precision);
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(DIM),
            AttributeRule.newIntegerRule(KNOTS),
            AttributeRule.newDoubleRule(DEGREE),
            AttributeRule.newDoubleArrayRule(TIMES),
            AttributeRule.newDoubleRule(BOUNDARY),
            new ElementRule(MEAN,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(MARGINALVARIANCE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(LENGTHSCALE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(COEFFICIENT,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(PRECISION,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
    };

    public String getParserDescription() { // TODO update
        return "Describes a normal distribution with a given mean and precision " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() { return RandomField.class; }
}