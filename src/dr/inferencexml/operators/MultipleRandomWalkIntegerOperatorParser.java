/*
 * MultipleRandomWalkIntegerOperatorParser.java
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

package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.MultipleRandomWalkIntegerOperator;
import dr.xml.*;

/**
 */
public class MultipleRandomWalkIntegerOperatorParser extends AbstractXMLObjectParser {
    public static final String MULTIPLE_RANDOM_WALK_INT_OP = "multipleRandomWalkIntegerOperator";
    public static final String SAMPLE_SIZE = "sampleSize";

    public String getParserName() {
        return MULTIPLE_RANDOM_WALK_INT_OP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        double w = xo.getDoubleAttribute(RandomWalkIntegerOperatorParser.WINDOW_SIZE);
        if (w != Math.floor(w)) {
            throw new XMLParseException("The window size of a randomWalkIntegerOperator should be an integer");
        }

        double s = xo.getDoubleAttribute(SAMPLE_SIZE);
        if (s != Math.floor(s)) {
            throw new XMLParseException("The window size of a randomWalkIntegerOperator should be an integer");
        }

        int windowSize = (int)w;
        int sampleSize = (int)s;
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        return new MultipleRandomWalkIntegerOperator(parameter, windowSize, sampleSize, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a random walk operator on a given parameter.";
    }

    public Class getReturnType() {
        return MultipleRandomWalkIntegerOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(RandomWalkIntegerOperatorParser.WINDOW_SIZE),
            AttributeRule.newDoubleRule(SAMPLE_SIZE),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(Parameter.class)
    };

}
