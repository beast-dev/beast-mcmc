/*
 * SwapOperatorParser.java
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
import dr.inference.operators.SwapOperator;
import dr.xml.*;

/**
 */
public class SwapOperatorParser extends AbstractXMLObjectParser {

    public final static String SWAP_OPERATOR = "swapOperator";

    public String getParserName() {
        return SWAP_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        double weight = xo.getDoubleAttribute("weight");
        int size = xo.getIntegerAttribute("size");

        boolean autoOptimize = xo.getBooleanAttribute("autoOptimize");
        if (autoOptimize) throw new XMLParseException("swapOperator can't be optimized!");

        System.out.println("Creating swap operator for parameter " + parameter.getParameterName() + " (weight=" + weight + ")");

        SwapOperator so = new SwapOperator(parameter, size);
        so.setWeight(weight);

        return so;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an operator that swaps values in a multi-dimensional parameter.";
    }

    public Class getReturnType() {
        return SwapOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{

            AttributeRule.newDoubleRule("weight"),
            AttributeRule.newIntegerRule("size"),
            AttributeRule.newBooleanRule("autoOptimize"),
            new ElementRule(Parameter.class)
    };

}
