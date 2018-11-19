/*
 * RateScaleOperatorParser.java
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

package dr.evomodelxml.operators;

import dr.evomodel.operators.RateScaleOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class RateScaleOperatorParser extends AbstractXMLObjectParser {

    public static final String SCALE_OPERATOR = "rateScaleOperator";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String NO_ROOT = "noRoot";

    public String getParserName() {
        return SCALE_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AdaptationMode mode = AdaptationMode.parseMode(xo);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

        final boolean noRoot = xo.getBooleanAttribute(NO_ROOT);

        if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
            throw new XMLParseException("scaleFactor must be between 0.0 and 1.0");
        }

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        RateScaleOperator operator = new RateScaleOperator(treeModel, scaleFactor, noRoot, mode);
        operator.setWeight(weight);
        return operator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a rateScale operator on a given parameter.";
    }

    public Class getReturnType() {
        return RateScaleOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(SCALE_FACTOR),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            AttributeRule.newBooleanRule(NO_ROOT, true),
            new ElementRule(TreeModel.class),
    };

}
