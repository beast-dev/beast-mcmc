/*
 * JointParameterParser.java
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

package dr.inferencexml.model;

import dr.inference.model.CompoundParameter;
import dr.inference.model.JointParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Andrew Rambaut
 */
public class JointParameterParser extends AbstractXMLObjectParser {

    public static final String COMPOUND_PARAMETER = "jointParameter";

    public String getParserName() {
        return COMPOUND_PARAMETER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        JointParameter jointParameter = new JointParameter((String) null);

        for (int i = 0; i < xo.getChildCount(); i++) {
            jointParameter.addParameter((Parameter) xo.getChild(i));
        }

        return jointParameter;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A parameter that synchronises its component parameters.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules;{
        rules = new XMLSyntaxRule[]{
                new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
        };
    }

    public Class getReturnType() {
        return JointParameter.class;
    }
}
