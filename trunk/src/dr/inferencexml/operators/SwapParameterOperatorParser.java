/*
 * SwapParameterOperatorParser.java
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
import dr.inference.operators.SwapParameterOperator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 */
public class SwapParameterOperatorParser extends AbstractXMLObjectParser {

    public final static String SWAP_OPERATOR = "swapParameterOperator";

    public String getParserName() {
        return SWAP_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<Parameter> parameterList = new ArrayList<Parameter>();

        for (int i = 0; i < xo.getChildCount(); ++i) {
            parameterList.add((Parameter) xo.getChild(i));
        }

        double weight = xo.getDoubleAttribute("weight");

        return new SwapParameterOperator(parameterList, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an operator that swaps values in a multi-dimensional parameter.";
    }

    public Class getReturnType() {
        return SwapParameterOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule("weight"),
            new ElementRule(Parameter.class, 2, Integer.MAX_VALUE)
    };

}
