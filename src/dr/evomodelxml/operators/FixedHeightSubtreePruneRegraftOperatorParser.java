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

import dr.evomodel.operators.FixedHeightSubtreePruneRegraftOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class FixedHeightSubtreePruneRegraftOperatorParser extends AbstractXMLObjectParser {

    public static final String FIXED_HEIGHT_SUBTREE_PRUNE_REGRAFT = "fixedHeightSubtreePruneRegraft";


    public String getParserName() {
        return FIXED_HEIGHT_SUBTREE_PRUNE_REGRAFT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        return new FixedHeightSubtreePruneRegraftOperator(treeModel, weight);
    }

    public String getParserDescription() {
        return "An operator that jumps a subtree to another edge at the same height.";
    }

    public Class getReturnType() {
        return FixedHeightSubtreePruneRegraftOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(TreeModel.class)
    };

}
