/*
 * SubtreeJumpOperatorParser.java
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

import dr.evomodel.operators.SubtreeJumpOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class SubtreeJumpOperatorParser extends AbstractXMLObjectParser {

    public static final String SUBTREE_JUMP = "subtreeJump";

    public static final String SIZE = "size";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";
    public static final String UNIFORM = "uniform";

    
    public String getParserName() {
        return SUBTREE_JUMP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CoercionMode mode = CoercionMode.parseMode(xo);

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        // size attribute is mandatory
        final double size = xo.getAttribute(SIZE, Double.NaN);
        final double targetAcceptance = xo.getAttribute(TARGET_ACCEPTANCE, 0.234);
        final boolean uniform = xo.getAttribute(UNIFORM, false);

        if (size <= 0.0) {
            throw new XMLParseException("The SubTreeLeap size attribute must be positive and non-zero.");
        }

        if (Double.isInfinite(size)) {
            // uniform so no auto optimize
            mode = CoercionMode.COERCION_OFF;
        }

        if (targetAcceptance <= 0.0 || targetAcceptance >= 1.0) {
            throw new XMLParseException("Target acceptance probability has to lie in (0, 1)");
        }

        SubtreeJumpOperator operator = new SubtreeJumpOperator(treeModel, weight, size, targetAcceptance, uniform, mode);
        operator.setTargetAcceptanceProbability(targetAcceptance);

        return operator;
    }

    public String getParserDescription() {
        return "An operator that jumps a subtree to another edge at the same height.";
    }

    public Class getReturnType() {
        return SubtreeJumpOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(SIZE, true),
            AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true),
            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(TreeModel.class)
    };

}
