/*
 * SelectorOperatorParser.java
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
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SelectorOperator;
import dr.xml.*;

/**
 *
 */
public class SelectorOperatorParser extends AbstractXMLObjectParser {
    public static String SELECTOR_OPERATOR = "selectorOperator";

    public String getParserName() {
        return SELECTOR_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final SelectorOperator op = new SelectorOperator(parameter);
        op.setWeight(weight);
        return op;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return ".";
    }

    public Class getReturnType() {
        return SelectorOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(Parameter.class),
        };
    }
}
