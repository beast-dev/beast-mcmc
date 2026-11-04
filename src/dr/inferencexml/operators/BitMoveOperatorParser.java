/*
 * BitMoveOperatorParser.java
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
import dr.inference.operators.BitMoveOperator;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 *
 */
public class BitMoveOperatorParser extends AbstractXMLObjectParser {

    public static final String BIT_MOVE_OPERATOR = "bitMoveOperator";
    public static final String NUM_BITS_TO_MOVE = "numBitsToMove";

    public String getParserName() {
        return BIT_MOVE_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        int numBitsToMove = xo.getIntegerAttribute(NUM_BITS_TO_MOVE);

        Parameter bitsParameter = (Parameter) xo.getElementFirstChild("bits");
        Parameter valuesParameter = null;


        if (xo.hasChildNamed("values")) {
            valuesParameter = (Parameter) xo.getElementFirstChild("values");
        }


        return new BitMoveOperator(bitsParameter, valuesParameter, numBitsToMove, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a bit-move operator on a given parameter.";
    }

    public Class getReturnType() {
        return BitMoveOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newIntegerRule(NUM_BITS_TO_MOVE),
            new ElementRule("bits", Parameter.class),
            new ElementRule("values", Parameter.class, "values parameter", true)
    };

}
