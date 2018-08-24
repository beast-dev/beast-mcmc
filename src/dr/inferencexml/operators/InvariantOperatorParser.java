/*
 * InvariantOperatorParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.InvariantOperator;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 *
 */
public class InvariantOperatorParser extends AbstractXMLObjectParser {

    public static final String OPERATOR_NAME = "invariantOperator";
    private static final String CHECK_LIKELIHOOD = "checkLikelihood";
    private static final String TRANSLATE = "translate";
    private static final String ROTATE = "rotate";
    private static final String DIMENSION = "dimension";

    public String getParserName() {
        return OPERATOR_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);

        boolean translate = xo.getAttribute(TRANSLATE, true);
        boolean rotate = xo.getAttribute(ROTATE, true);
        boolean checkLikelihood = xo.getAttribute(CHECK_LIKELIHOOD, true);

        int dim = xo.getIntegerAttribute(DIMENSION);

        return new InvariantOperator.Rotation(parameter, dim, weight, likelihood,
                translate, rotate, checkLikelihood);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a operator that forces the entire model likelihood recomputation";
    }

    public Class getReturnType() {
        return InvariantOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(Parameter.class),
            AttributeRule.newBooleanRule(TRANSLATE, true),
            AttributeRule.newBooleanRule(ROTATE, true),
            AttributeRule.newBooleanRule(CHECK_LIKELIHOOD, true),
            AttributeRule.newIntegerRule(DIMENSION),
            new ElementRule(Likelihood.class, true),
    };
}
