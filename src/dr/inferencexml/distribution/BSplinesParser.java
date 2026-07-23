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
import dr.math.distributions.BSplines;
import dr.xml.*;


public class BSplinesParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "BSplines";
    private static final String KNOTS = "knots";
    private static final String DEGREE = "degree";
    private static final String TIMES = "times";
    private static final String LOWERBOUNDARY = "lowerBoundary";
    private static final String UPPERBOUNDARY = "upperBoundary";
    private static final String COEFFICIENT = "coefficient";
    private static final String PRECISION = "precision";

    public String getParserName() { return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int degree = xo.getIntegerAttribute(DEGREE);

        double[] knots = xo.getDoubleArrayAttribute(KNOTS);

        double[] times = xo.getDoubleArrayAttribute(TIMES);

        double lowerBoundary = xo.getDoubleAttribute(LOWERBOUNDARY);

        double upperBoundary = xo.getDoubleAttribute(UPPERBOUNDARY);

        Parameter coefficient = (Parameter) xo.getElementFirstChild(COEFFICIENT);

        Parameter precision = (Parameter) xo.getElementFirstChild(PRECISION);

        if (precision.getParameterValue(0) <= 0.0) {
            throw new XMLParseException("Noise must be > 0.0");
        }


//        RandomField.WeightProvider weights = parseWeightProvider(xo, dim);


        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        return new BSplines(id, knots, times, degree, lowerBoundary, upperBoundary, coefficient, precision);
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(DEGREE),
            AttributeRule.newDoubleArrayRule(KNOTS),
            AttributeRule.newDoubleArrayRule(TIMES),
            AttributeRule.newDoubleRule(LOWERBOUNDARY),
            AttributeRule.newDoubleRule(UPPERBOUNDARY),
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