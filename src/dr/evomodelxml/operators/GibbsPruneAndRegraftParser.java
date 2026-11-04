/*
 * GibbsPruneAndRegraftParser.java
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

package dr.evomodelxml.operators;

import dr.evomodel.operators.GibbsPruneAndRegraft;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class GibbsPruneAndRegraftParser extends AbstractXMLObjectParser {

    public static final String GIBBS_PRUNE_AND_REGRAFT = "GibbsPruneAndRegraft";

    public String getParserName() {
        return GIBBS_PRUNE_AND_REGRAFT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        boolean pruned = true;
        if (xo.hasAttribute("pruned")) {
            pruned = xo.getBooleanAttribute("pruned");
        }

        return new GibbsPruneAndRegraft(treeModel, pruned, weight);
    }

    // ************************************************************************
    // AbstractXMLObjectParser implementation
    // ************************************************************************

    public String getParserDescription() {
        return "This element represents a Gibbs sampler implemented through a prune and regraft operator. "
                + "This operator prunes a random subtree and regrafts it below a node chosen by an importance distribution which is the proportion of the likelihoods of the proposals.";
    }

    public Class getReturnType() {
        return GibbsPruneAndRegraft.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules;

    {
        rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newBooleanRule("pruned"),
                new ElementRule(TreeModel.class)};
    }

}
