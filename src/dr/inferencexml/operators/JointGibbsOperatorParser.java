/*
 * JointGibbsOperatorParser.java
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

import dr.inference.operators.GibbsOperator;
import dr.inference.operators.JointGibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class JointGibbsOperatorParser extends AbstractXMLObjectParser {

    public static final String JOINT_GIBBS_OPERATOR = "jointGibbsOperator";
    public static final String WEIGHT = "weight";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";

    public String getParserName() {
        return JOINT_GIBBS_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(WEIGHT);

        JointGibbsOperator operator = new JointGibbsOperator(weight);

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof  GibbsOperator)
                operator.addOperator((SimpleMCMCOperator) xo.getChild(i));
            else
                throw new RuntimeException("Operator list must consist only of GibbsOperators");
        }

        return operator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a list of Gibbs Operators";
    }

    public Class getReturnType() {
        return JointGibbsOperator.class;
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
