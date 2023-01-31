/*
 * RateExchangeOperatorParser.java
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

import dr.evomodel.operators.RateExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class RateExchangeOperatorParser extends AbstractXMLObjectParser {

    public static final String RATE_EXCHANGE = "rateExchange";
    public static final String SWAP_TRAITS = "swapTraits";
    public static final String SWAP_RATES = "swapRates";
    public static final String SWAP_AT_ROOT = "swapAtRoot";
    public static final String MOVE_HEIGHT = "moveHeight";

    public String getParserName() {
        return RATE_EXCHANGE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DefaultTreeModel treeModel = (DefaultTreeModel) xo.getChild(DefaultTreeModel.class);
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        boolean swapRates = xo.getBooleanAttribute(SWAP_RATES);
        boolean swapTraits = xo.getBooleanAttribute(SWAP_TRAITS);
        boolean swapAtRoot = xo.getBooleanAttribute(SWAP_AT_ROOT);
        boolean moveHeight = xo.getBooleanAttribute(MOVE_HEIGHT);
        return new RateExchangeOperator(treeModel, weight, swapRates, swapTraits, swapAtRoot, moveHeight);
    }

    public String getParserDescription() {
        return "An operator that exchanges rates and traits on a tree.";
    }

    public Class getReturnType() {
        return RateExchangeOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(SWAP_RATES),
            AttributeRule.newBooleanRule(SWAP_TRAITS),
            AttributeRule.newBooleanRule(SWAP_AT_ROOT),
            AttributeRule.newBooleanRule(MOVE_HEIGHT),
            new ElementRule(DefaultTreeModel.class)
    };
}
