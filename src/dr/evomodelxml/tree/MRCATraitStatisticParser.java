/*
 * MRCATraitStatisticParser.java
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

import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.MRCATraitStatistic;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 */
public class MRCATraitStatisticParser extends AbstractXMLObjectParser {

    public static final String MRCA_TRAIT_STATISTIC = "mrcaTraitStatistic";
    public static final String MRCA = "mrca";
    public static final String NAME = "name";
    public static final String TRAIT = "trait";

    public String getParserName() {
        return MRCA_TRAIT_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(NAME, xo.getId());
        String trait = xo.getStringAttribute(TRAIT);

        DefaultTreeModel tree = (DefaultTreeModel) xo.getChild(DefaultTreeModel.class);
        TaxonList taxa = (TaxonList) xo.getElementFirstChild(MRCA);

        try {
            return new MRCATraitStatistic(name, trait, tree, taxa);
        } catch (TreeUtils.MissingTaxonException mte) {
            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that has as its value the height of the most recent common ancestor of a set of taxa in a given tree";
    }

    public Class getReturnType() {
        return MRCATraitStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(DefaultTreeModel.class),
            new StringAttributeRule("name", "A name for this statistic primarily for the purposes of logging", true),
            new StringAttributeRule("trait", "The name of the trait (can be rate)"),
            AttributeRule.newBooleanRule("rate", true),
            new ElementRule(MRCA,
                    new XMLSyntaxRule[]{new ElementRule(Taxa.class)})
    };

}
