/*
 * UpDownOperatorParser.java
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
import dr.inference.operators.*;
import dr.xml.*;

/**
 */
public class FreeRateUpDownOperatorParser extends AbstractXMLObjectParser {

    public static final String FREERATE_UP_DOWN_OPERATOR = "freeRateUpDownOperator";

    public static final String SCALE_FACTOR = ScaleOperatorParser.SCALE_FACTOR;

    public static final String WEIGHTS_PARAMETER = "weights";
    public static final String RATES_PARAMETER = "rates";

    public String getParserName() {
        return FREERATE_UP_DOWN_OPERATOR;
    }


    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AdaptationMode mode = AdaptationMode.parseMode(xo);

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

        if( (scaleFactor <= 0.0) | (scaleFactor >= 1.0) ) {
            throw new IllegalArgumentException("scale must be between 0 and 1");
        }

        Parameter weightsParameter = (Parameter) xo.getElementFirstChild(WEIGHTS_PARAMETER);
        Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RATES_PARAMETER);

        return new FreeRateUpDownOperator(weightsParameter, ratesParameter, scaleFactor, weight, mode);
    }

    public String getParserDescription() {
        return "This element represents an operator that scales two elements of the rates parameter of a FreeRate model" +
                " such that the constraints are followed. Weights are not touched.";
    }

    public Class getReturnType() {
        return FreeRateUpDownOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(SCALE_FACTOR),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(WEIGHTS_PARAMETER, Parameter.class),
            new ElementRule(RATES_PARAMETER, Parameter.class)
    };

}
