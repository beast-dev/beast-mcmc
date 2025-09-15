/*
 * DeltaExchangeOperatorParser.java
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
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.FreeRateDeltaExchangeOperator;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class FreeRateDeltaExchangeOperatorParser extends AbstractXMLObjectParser {

    public static final String FREE_RATE_DELTA_EXCHANGE = "freeRateDeltaExchange";
    public static final String DELTA = "delta";
    public static final String WEIGHTS_PARAMETER = "weights";
    public static final String RATES_PARAMETER = "rates";

    public String getParserName() {
        return FREE_RATE_DELTA_EXCHANGE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AdaptationMode mode = AdaptationMode.parseMode(xo);

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        double delta = xo.getDoubleAttribute(DELTA);

        if (delta <= 0.0) {
            throw new XMLParseException("delta must be greater than 0.0");
        }

        Parameter weightsParameter = (Parameter) xo.getElementFirstChild(WEIGHTS_PARAMETER);
        Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RATES_PARAMETER);

        return new FreeRateDeltaExchangeOperator(weightsParameter, ratesParameter, delta, weight, mode);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a delta exchange operator on the weights of a FreeRate model. One rate is also adjusted to maintain the constraints.";
    }

    public Class getReturnType() {
        return MCMCOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(DELTA),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(WEIGHTS_PARAMETER, Parameter.class),
            new ElementRule(RATES_PARAMETER, Parameter.class)
    };
}
