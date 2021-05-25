/*
 * TreeBitRandomWalkOperatorParser.java
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

import dr.evomodel.operators.TreeBitRandomWalkOperator;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class TreeBitRandomWalkOperatorParser extends AbstractXMLObjectParser {

    public static final String BIT_RANDOM_WALK_OPERATOR = "treeBitRandomWalk";
    public static final String INDICTATOR_TRAIT = "indicatorTrait";
    public static final String TRAIT2 = "trait2";
    public static final String SWAP_TRAIT2 = "swapTrait2";

    public String getParserName() {
        return BIT_RANDOM_WALK_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        DefaultTreeModel treeModel = (DefaultTreeModel) xo.getChild(DefaultTreeModel.class);


        String trait1 = null;
        String trait2 = null;
        if (xo.hasAttribute(INDICTATOR_TRAIT)) trait1 = xo.getStringAttribute(INDICTATOR_TRAIT);
        if (xo.hasAttribute(TRAIT2)) trait2 = xo.getStringAttribute(TRAIT2);
        int k = xo.getAttribute("k", 1);
        boolean swapTrait2 = xo.getAttribute(SWAP_TRAIT2, true);

        return new TreeBitRandomWalkOperator(treeModel, trait1, trait2, weight, k, swapTrait2);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a bit-random walk operator on a random " +
                "indicator/variable pair in the tree.";
    }

    public Class getReturnType() {
        return MCMCOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(DefaultTreeModel.class),
            AttributeRule.newStringRule(INDICTATOR_TRAIT, true),
            AttributeRule.newStringRule(TRAIT2, true),
            AttributeRule.newBooleanRule(SWAP_TRAIT2, true),
            AttributeRule.newIntegerRule("k", true)
    };
}
