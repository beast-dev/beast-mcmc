/*
 * SubtreeRateStatisticParser.java
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

package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.SubtreeRateStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Statistic;
import dr.xml.*;

import static dr.evomodelxml.tree.MonophylyStatisticParser.parseTaxonListOrTaxa;

/**
 */
public class SubtreeRateStatisticParser extends AbstractXMLObjectParser {

    public static final String RATE_STATISTIC = "SubtreeRateStatistic";
    public static final String MODE = "mode";
    public static final String MEAN = "mean";
    public static final String VARIANCE = "variance";
    public static final String COEFFICIENT_OF_VARIATION = "coefficientOfVariation";
    private static final String MRCA = "mrca";

    public String getParserName() {
        return RATE_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String name = xo.getAttribute(Statistic.NAME, xo.getId());

        final Tree tree = (Tree) xo.getChild(Tree.class);

        final BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        TaxonList mrcaTaxa = null;
        if (xo.hasChildNamed(MRCA)) {
            mrcaTaxa = parseTaxonListOrTaxa(xo.getChild(MRCA));
        }

        final boolean complement = xo.getAttribute("complement", false);

        final boolean includeStem = xo.getAttribute("includeStem", false);

        final String mode = xo.getStringAttribute(MODE);

        try {
            return new SubtreeRateStatistic(name, tree, branchRateModel, mode, mrcaTaxa, complement, includeStem);
        } catch (TreeUtils.MissingTaxonException e) {
            throw new XMLParseException("Unable to find taxon-set.");
        }

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that returns the average, variance, or coefficient of variations of the branch rates in a subtree (or its complement).";
    }

    public Class getReturnType() {
        return SubtreeRateStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new StringAttributeRule("name", "A name for this statistic primarily for the purposes of logging", true),
            new ElementRule(TreeModel.class, "The tree."),
            new ElementRule(BranchRateModel.class, "The model giving the branch rates to track on the tree."),
            new StringAttributeRule("mode", "This attribute determines how the rates are summarized, can be one of (mean, variance, coefficientOfVariance)", new String[]{MEAN, VARIANCE, COEFFICIENT_OF_VARIATION}, false),
            new ElementRule(MRCA,
                    new XMLSyntaxRule[]{new ElementRule(Taxa.class)}, "Taxa that define the subtree.", false),
            AttributeRule.newBooleanRule("complement", true, "If true, computes statistics for all branches not in the subtree."),
            AttributeRule.newBooleanRule("includeStem", true, "Should stem branch be included in the subtree defined by the MRCA of the taxa?")
    };

}
