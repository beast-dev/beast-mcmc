/*
 * CompoundParameterParser.java
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

package dr.inferencexml.model;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class CompoundParameterParser extends AbstractXMLObjectParser {

    public static final String COMPOUND_PARAMETER = "compoundParameter";

    public String getParserName() {
        return COMPOUND_PARAMETER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CompoundParameter compoundParameter = new CompoundParameter((String) null);

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        if (treeModel != null) {
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            for (int i = 0; i < treeModel.getNodeCount() - 1; i++) {
                compoundParameter.addParameter(new Parameter.Default(parameter.getParameterValue(0)));
            }
        } else {
            for (int i = 0; i < xo.getChildCount(); i++) {
                compoundParameter.addParameter((Parameter) xo.getChild(i));
            }
        }

        return compoundParameter;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A multidimensional parameter constructed from its component parameters.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules;{
        rules = new XMLSyntaxRule[]{
                new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
                new ElementRule(TreeModel.class, true)
        };
    }

    public Class getReturnType() {
        return CompoundParameter.class;
    }
}
