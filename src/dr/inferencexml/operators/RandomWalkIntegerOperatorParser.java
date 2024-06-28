/*
 * RandomWalkIntegerOperatorParser.java
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
import dr.inference.model.Variable;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.RandomWalkIntegerOperator;
import dr.xml.*;

/**
 */
public class RandomWalkIntegerOperatorParser extends AbstractXMLObjectParser {

    public static final String RANDOM_WALK_INTEGER_OPERATOR = "randomWalkIntegerOperator";

    public static final String WINDOW_SIZE = "windowSize";
    public static final String UPDATE_INDEX = "updateIndex";

    public String getParserName() {
        return RANDOM_WALK_INTEGER_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        double d = xo.getDoubleAttribute(WINDOW_SIZE);
        if (d != Math.floor(d)) {
            throw new XMLParseException("The window size of a " + RANDOM_WALK_INTEGER_OPERATOR + " should be an integer");
        }

        int windowSize = (int)d;
        Variable parameter = (Variable) xo.getChild(Variable.class);

        return new RandomWalkIntegerOperator(parameter, windowSize, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a random walk operator on a given parameter.";
    }

    public Class getReturnType() {
        return RandomWalkIntegerOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WINDOW_SIZE),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(Variable.class)
    };

}
