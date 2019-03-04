/*
 * DeltaExchangeOperatorParser.java
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

package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.DeltaExchangeOperator;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class DeltaExchangeOperatorParser extends AbstractXMLObjectParser {

    public static final String DELTA_EXCHANGE = "deltaExchange";
    public static final String DELTA = "delta";
    public static final String INTEGER_OPERATOR = "integer";
    public static final String PARAMETER_WEIGHTS = "parameterWeights";

    public String getParserName() {
        return DELTA_EXCHANGE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AdaptationMode mode = AdaptationMode.parseMode(xo);

        final boolean isIntegerOperator = xo.getAttribute(INTEGER_OPERATOR, false);

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        double delta = xo.getDoubleAttribute(DELTA);

        if (delta <= 0.0) {
            throw new XMLParseException("delta must be greater than 0.0");
        }

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);


        int[] parameterWeights;
        if (xo.hasAttribute(PARAMETER_WEIGHTS)) {
            parameterWeights = xo.getIntegerArrayAttribute(PARAMETER_WEIGHTS);
            System.out.print("Parameter weights for delta exchange are: ");
            for (int parameterWeight : parameterWeights) {
                System.out.print(parameterWeight + "\t");
            }
            System.out.println();

        } else {
            parameterWeights = new int[parameter.getDimension()];
            for (int i = 0; i < parameterWeights.length; i++) {
                parameterWeights[i] = 1;
            }
        }

        if (parameterWeights.length != parameter.getDimension()) {
            throw new XMLParseException("parameter weights have the same length as parameter");
        }


        return new DeltaExchangeOperator(parameter, parameterWeights, delta, weight, isIntegerOperator, mode);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a delta exchange operator on a given parameter.";
    }

    public Class getReturnType() {
        return MCMCOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(DELTA),
            AttributeRule.newIntegerArrayRule(PARAMETER_WEIGHTS, true),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            AttributeRule.newBooleanRule(INTEGER_OPERATOR, true),
            new ElementRule(Parameter.class)
    };
}
