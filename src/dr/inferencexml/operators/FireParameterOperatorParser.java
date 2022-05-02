/*
 * FireParameterOperatorParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.inference.operators.DirtyLikelihoodOperator;
import dr.inference.operators.FireParameterOperator;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 *
 */
public class FireParameterOperatorParser extends AbstractXMLObjectParser {

    private static final String VALUE = "value";
    private static final String COPY_FROM = "copyFrom";

    public static final String FIRE_PARAMETER_OPERATOR = "fireParameterChanged";

    public String getParserName() {
        return FIRE_PARAMETER_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        double[] values = null;
        if (xo.hasAttribute(VALUE)) {
            values = xo.getDoubleArrayAttribute(VALUE);

        } else if (xo.hasChildNamed(COPY_FROM)) {
            Parameter copyFromParameter = (Parameter) xo.getChild(COPY_FROM).getChild(0);
            values = copyFromParameter.getParameterValues();
        }

        if (values != null) {
            System.out.println("\nWarning: when the operator " + FIRE_PARAMETER_OPERATOR + " is given a \"" + VALUE +
                    "\" attribute or \"" + COPY_FROM + "\" element, the resulting MCMC run will NOT result in a " +
                    "valid draw from the posterior. " +
                    "Only set the \"" + VALUE + "\" attribute for debugging purposes.\n");
        }

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        return new FireParameterOperator(parameter, values, weight);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a operator that forces the entire model likelihood recomputation";
    }

    public Class getReturnType() {
        return DirtyLikelihoodOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(Parameter.class),
            new XORRule(
                    AttributeRule.newDoubleArrayRule(VALUE),
                    new ElementRule(COPY_FROM, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                    true)
    };
}
