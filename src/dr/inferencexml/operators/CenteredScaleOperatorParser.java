/*
 * CenteredScaleOperatorParser.java
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
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.CenteredScaleOperator;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class CenteredScaleOperatorParser extends AbstractXMLObjectParser {

    public static final String CENTERED_SCALE = "centeredScale";
    public static final String SCALE_FACTOR = "scaleFactor";

    public String getParserName() {
        return CENTERED_SCALE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        double scale = xo.getDoubleAttribute(SCALE_FACTOR);
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        CenteredScaleOperator op = new CenteredScaleOperator(parameter);
        op.setWeight(weight);
        op.scaleFactor = scale;
        op.mode = AdaptationMode.parseMode(xo);
        return op;
    }

    public String getParserDescription() {
        return "A centered-scale operator. This operator scales the the values of a multi-dimensional parameter so as to perserve the mean. It does this by expanding or conrtacting the parameter values around the mean.";
    }

    public Class getReturnType() {
        return CenteredScaleOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(SCALE_FACTOR),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(Parameter.class)
    };

}
