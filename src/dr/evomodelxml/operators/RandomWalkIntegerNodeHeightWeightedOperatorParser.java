/*
 * RandomWalkIntegerNodeHeightWeightedOperatorParser.java
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

package dr.evomodelxml.operators;

import dr.xml.*;
import dr.inference.operators.MCMCOperator;
import dr.evomodel.operators.RandomWalkIntegerNodeHeightWeightedOperator;
import dr.inference.model.Parameter;

/**
 * @author Chieh-Hsi Wu
 *
 * The parser for random walk integer node height weighted operator.
 */
public class RandomWalkIntegerNodeHeightWeightedOperatorParser extends AbstractXMLObjectParser {

    public static final String RANDOM_WALK_INT_NODE_HEIGHT_WGT_OP = "randomWalkIntegerNodeHeightWeightedOperator";

    public static final String WINDOW_SIZE = "windowSize";
    public static final String INTERNAL_NODE_HEIGHTS = "internalNodeHeights";

    public String getParserName() {
        return RANDOM_WALK_INT_NODE_HEIGHT_WGT_OP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        double d = xo.getDoubleAttribute(WINDOW_SIZE);
        if (d != Math.floor(d)) {
            throw new XMLParseException("The window size of a " + RANDOM_WALK_INT_NODE_HEIGHT_WGT_OP + " should be an integer");
        }

        int windowSize = (int)d;
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        Parameter internalNodeHeights = (Parameter)xo.getElementFirstChild(INTERNAL_NODE_HEIGHTS);
        
        return new RandomWalkIntegerNodeHeightWeightedOperator(parameter, windowSize, weight, internalNodeHeights);
    }

    public String getParserDescription() {
        return "This element returns a random walk node height weighted operator on a given parameter.";
    }

    public Class getReturnType() {
        return RandomWalkIntegerNodeHeightWeightedOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WINDOW_SIZE),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(Parameter.class),
            new ElementRule(INTERNAL_NODE_HEIGHTS, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };
}
