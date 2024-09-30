/*
 * GaussianMarkovRandomFieldParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.RandomField;
import dr.inference.model.Parameter;
import dr.math.distributions.GaussianMarkovRandomField;
import dr.xml.*;

import static dr.inferencexml.distribution.RandomFieldParser.WEIGHTS_RULE;
import static dr.inferencexml.distribution.RandomFieldParser.parseWeightProvider;

public class GaussianMarkovRandomFieldParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "gaussianMarkovRandomField";
    private static final String DIMENSION = "dim";
    private static final String PRECISION = "precision";
    private static final String MEAN = "mean";
    private static final String LAMBDA = "lambda";
    private static final String MATCH_PSEUDO_DETERMINANT = "matchPseudoDeterminant";

    public String getParserName() { return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int dim = xo.getIntegerAttribute(DIMENSION);

        Parameter incrementPrecision = (Parameter) xo.getElementFirstChild(PRECISION);

        if (incrementPrecision.getParameterValue(0) <= 0.0) {
            throw new XMLParseException("Scale must be > 0.0");
        }

        Parameter start = xo.hasChildNamed(MEAN) ?
                (Parameter) xo.getElementFirstChild(MEAN) : null;

        Parameter lambda = xo.hasChildNamed(LAMBDA) ?
                (Parameter) xo.getElementFirstChild(LAMBDA) : null;

//        RandomField.WeightProvider weights = parseWeightProvider(xo, dim);

        RandomField.WeightProvider weights = (RandomField.WeightProvider) xo.getChild(RandomField.WeightProvider.class);

        boolean matchPseudoDeterminant = xo.getAttribute(MATCH_PSEUDO_DETERMINANT, false);

        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        return new GaussianMarkovRandomField(id, dim, incrementPrecision, start, lambda, weights, matchPseudoDeterminant);
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(DIMENSION),
            new ElementRule(RandomField.WeightProvider.class, true),
            new ElementRule(PRECISION,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(MEAN,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(LAMBDA,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            AttributeRule.newBooleanRule(MATCH_PSEUDO_DETERMINANT, true),
            WEIGHTS_RULE,
    };

    public String getParserDescription() { // TODO update
        return "Describes a normal distribution with a given mean and precision " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() { return RandomField.class; }
}
