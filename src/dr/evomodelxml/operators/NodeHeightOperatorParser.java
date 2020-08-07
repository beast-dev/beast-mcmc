/*
 * SubtreeLeapOperatorParser.java
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

import dr.evomodel.operators.NodeHeightOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * @author Andrew Rambaut
 */
public class NodeHeightOperatorParser extends AbstractXMLObjectParser {

    public static final String NODE_HEIGHT_OPERATOR = "nodeHeightOperator";

    public static final String SIZE = "size";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";
    public static final String OPERATOR_TYPE = "type";

    public String getParserName() {
        return NODE_HEIGHT_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AdaptationMode mode = AdaptationMode.parseMode(xo);

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        double tuningParameter = 0.75;

        if (xo.hasAttribute(SIZE)) {
            tuningParameter = xo.getDoubleAttribute(SIZE);
            if (tuningParameter <= 0.0) {
                throw new XMLParseException("The NodeHeightOperator size attribute must be positive and non-zero.");
            }
        }
        if (xo.hasAttribute(SCALE_FACTOR)) {
            tuningParameter = xo.getDoubleAttribute(SCALE_FACTOR);
            if (tuningParameter <= 0.0 || tuningParameter >= 1.0) {
                throw new XMLParseException("The NodeHeightOperator scaleFactor attribute must be between 0 and 1.");
            }
        }

        final double targetAcceptance = xo.getAttribute(TARGET_ACCEPTANCE, 0.234);

        NodeHeightOperator.OperatorType operatorType = NodeHeightOperator.OperatorType.UNIFORM;
        if (xo.hasAttribute(OPERATOR_TYPE)) {
            try {
                operatorType = NodeHeightOperator.OperatorType.valueOf(xo.getStringAttribute(OPERATOR_TYPE).trim().toUpperCase());
            } catch (IllegalArgumentException iae) {
                throw new XMLParseException("Unrecognised operator type attribute: " + xo.getStringAttribute(OPERATOR_TYPE));
            }
        }


        if (targetAcceptance <= 0.0 || targetAcceptance >= 1.0) {
            throw new XMLParseException("Target acceptance probability has to lie in (0, 1)");
        }

        return new NodeHeightOperator(treeModel, weight, tuningParameter, operatorType, mode, targetAcceptance);
    }

    public String getParserDescription() {
        return "An operator that moves node heights about.";
    }

    public Class getReturnType() {
        return NodeHeightOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(SIZE, true),
            AttributeRule.newDoubleRule(SCALE_FACTOR, true),
            AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true),
            AttributeRule.newStringRule(OPERATOR_TYPE, true),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(TreeModel.class)
    };

}