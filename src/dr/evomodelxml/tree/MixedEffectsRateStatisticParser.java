/*
 * MixedEffectsRateStatisticParser.java
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

package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.tree.MixedEffectsRateStatistic;
import dr.evomodel.tree.RateStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class MixedEffectsRateStatisticParser extends AbstractXMLObjectParser {

    public static final String RATE_STATISTIC = "mixedEffectsRateStatistic";
    public static final String MODE = "mode";
    public static final String MEAN = "meanOfResiduals";
    public static final String VARIANCE = "varianceOfResiduals";
    public static final String PROP_EXPLAINED = "proportionOfVarianceExplained";

    public String getParserName() {
        return RATE_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String name = xo.getAttribute(Statistic.NAME, xo.getId());
        final Tree tree = (Tree) xo.getChild(Tree.class);
        final ArbitraryBranchRates branchRateModel = (ArbitraryBranchRates) xo.getChild(ArbitraryBranchRates.class);

        final boolean internal = xo.getAttribute("internal", true);
        final boolean external = xo.getAttribute("external", true);
        final boolean logScale = xo.getAttribute("logScale", true);

        if (!(internal || external)) {
            throw new XMLParseException("At least one of internal and external must be true!");
        }

        final String mode = xo.getStringAttribute(MODE);

        return new MixedEffectsRateStatistic(name, tree, branchRateModel, external, internal, logScale, mode);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that returns information about mixed-effects clock models.";
    }

    public Class getReturnType() {
        return MixedEffectsRateStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(ArbitraryBranchRates.class),
            AttributeRule.newBooleanRule("internal",true, "If true (default), the statistic considers internal branches. If false, internal branches are excluded from consideration."),
            AttributeRule.newBooleanRule("external",true, "If true (default), the statistic considers external branches. If false, external branches are excluded from consideration."),
            AttributeRule.newBooleanRule("logScale",true,"If true (default), rates and residuals are considered on the log scale."),
            new StringAttributeRule("mode", "This attribute determines how the rate residuals are summarized, can be the mean, variance, or proportion of explained variance.", new String[]{MEAN, VARIANCE, PROP_EXPLAINED}, false),
            new StringAttributeRule("name", "A name for this statistic primarily for the purposes of logging", true),
    };

}
