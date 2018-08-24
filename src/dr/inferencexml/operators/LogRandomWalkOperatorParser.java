/*
 * LogRandomWalkOperatorParser.java
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
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.LogRandomWalkOperator;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
@Deprecated
public class LogRandomWalkOperatorParser extends AbstractXMLObjectParser {

    public static final String LOGRANDOMWALK_OPERATOR = "logRandomWalkOperator";
    // Use same attributes as scale operator to help users
    public static final String SCALE_ALL = ScaleOperatorParser.SCALE_ALL;
    public static final String SCALE_ALL_IND = ScaleOperatorParser.SCALE_ALL_IND;
    public static final String WINDOW_SIZE = "window";

    public String getParserName() {
        return LOGRANDOMWALK_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CoercionMode mode = CoercionMode.parseMode(xo);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double size = xo.getDoubleAttribute(WINDOW_SIZE);
        final boolean scaleAll = xo.getAttribute(SCALE_ALL, false);
        final boolean scaleAllInd = xo.getAttribute(SCALE_ALL_IND, false);

        if( size <= 0.0 ) {
            throw new XMLParseException("size must be positive");
        }

        final Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        return new LogRandomWalkOperator(parameter, size, mode, weight, scaleAll, scaleAllInd);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a scale operator on a given parameter.";
    }

    public Class getReturnType() {
        return LogRandomWalkOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(WINDOW_SIZE),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),

            AttributeRule.newBooleanRule(SCALE_ALL, true),
            AttributeRule.newBooleanRule(SCALE_ALL_IND, true),

            new ElementRule(Parameter.class),
    };

}
