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
import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;
import dr.math.distributions.gp.GaussianProcessDistribution;
import dr.math.distributions.gp.Kernel;
import dr.xml.*;

import static dr.inferencexml.distribution.RandomFieldParser.WEIGHTS_RULE;
import static dr.inferencexml.distribution.RandomFieldParser.parseWeightProvider;

public class GaussianProcessFieldParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "gaussianProcessField";
    private static final String DIMENSION = "dim";
    private static final String MEAN = "mean";
    private static final String BASES = "bases";

    public String getParserName() { return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int dim = xo.getIntegerAttribute(DIMENSION);

        Parameter mean = xo.hasChildNamed(MEAN) ?
                (Parameter) xo.getElementFirstChild(MEAN) : null;

        RandomField.WeightProvider weights = parseWeightProvider(xo, dim);

        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        return new GaussianProcessDistribution(id, dim, mean,
                new Kernel.DotProduct(null, null),
                weights);
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(DIMENSION),
            new ElementRule(MEAN,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(BASES, new XMLSyntaxRule[] {
                    new ElementRule(DesignMatrix.class) },
                    0, Integer.MAX_VALUE),
            WEIGHTS_RULE,
    };

    public String getParserDescription() { // TODO update
        return "Describes a normal distribution with a given mean and precision " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() { return RandomField.class; }
}
