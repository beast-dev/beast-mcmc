/*
 * JointOperatorParser.java
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

import dr.inference.operators.JointOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

/**
 *
 */
public class JointOperatorParser extends AbstractXMLObjectParser {

    public static final String JOINT_OPERATOR = "jointOperator";
    public static final String WEIGHT = "weight";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";

    public String getParserName() {
        return JOINT_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(WEIGHT);

        final double targetProb = xo.getAttribute(TARGET_ACCEPTANCE, 0.2);

        if (targetProb <= 0.0 || targetProb >= 1.0)
            throw new RuntimeException("Target acceptance probability must be between 0.0 and 1.0");

        JointOperator operator = new JointOperator(weight, targetProb);

        for (int i = 0; i < xo.getChildCount(); i++) {
            operator.addOperator((SimpleMCMCOperator) xo.getChild(i));
        }

        return operator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an arbitrary list of operators; only the first is optimizable";
    }

    public Class getReturnType() {
        return JointOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SimpleMCMCOperator.class, 1, Integer.MAX_VALUE),
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true)
    };
}
