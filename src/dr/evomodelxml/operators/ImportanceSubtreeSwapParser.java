/*
 * ImportanceSubtreeSwapParser.java
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

import dr.evomodel.operators.ImportanceSubtreeSwap;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class ImportanceSubtreeSwapParser extends AbstractXMLObjectParser {

    public static final String IMPORTANCE_SUBTREE_SWAP = "ImportanceSubtreeSwap";

    public String getParserName() {
        return IMPORTANCE_SUBTREE_SWAP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final int samples = xo.getIntegerAttribute("samples");

        return new ImportanceSubtreeSwap(treeModel, weight, samples);
    }

    // ************************************************************************
    // AbstractXMLObjectParser implementation
    // ************************************************************************

    public String getParserDescription() {
        return "This element represents a importance guided subtree swap operator. "
                + "This operator swaps a random subtree with a second subtree guided by an importance distribution.";
    }

    public Class getReturnType() {
        return ImportanceSubtreeSwap.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules;{
        rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newIntegerRule("samples"),
                new ElementRule(TreeModel.class)
        };
    }

}
