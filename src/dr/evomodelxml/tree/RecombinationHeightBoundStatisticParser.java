/*
 * RecombinantHeightBoundStatisticParser.java
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
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.RecombinationHeightBoundStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 * @author Marc A Suchard
 * @author Philippe Lemey
 */
public class RecombinationHeightBoundStatisticParser extends AbstractXMLObjectParser {

    public static final String RECOMBINATION_HEIGHT_STATISTIC = "recombinationHeightStatistic";
    public static final String RECOMBINANT1 = "recombinant1";
    public static final String RECOMBINANT2 = "recombinant2";
    public static final String ABSOLUTE = TMRCAStatisticParser.ABSOLUTE;

    public String getParserName() {
        return RECOMBINATION_HEIGHT_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());

        Tree tree = (Tree) xo.getChild(Tree.class);

        TaxonList recombinant1 = (TaxonList) xo.getElementFirstChild(RECOMBINANT1); // TODO might break with list of taxa
        TaxonList recombinant2 = (TaxonList) xo.getElementFirstChild(RECOMBINANT2); // TODO might break with list of taxa

        boolean isAbsolute = xo.getAttribute(ABSOLUTE, false);

        try {
            return new RecombinationHeightBoundStatistic(name, tree, recombinant1, recombinant2, isAbsolute);
        } catch (TreeUtils.MissingTaxonException mte) {
            throw new XMLParseException(
                    "Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
        }
    }

    public String getParserDescription() {
        return "A statistic that has as its value the height of the most recent common ancestor " +
                "of a set of taxa in a given tree. ";
    }

    public Class getReturnType() {
        return RecombinationHeightBoundStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Tree.class),
            new StringAttributeRule("name",
                    "A name for this statistic primarily for the purposes of logging", true),
            AttributeRule.newBooleanRule(ABSOLUTE, true),
            new ElementRule(RECOMBINANT1,
                    new XMLSyntaxRule[]{
                            new ElementRule(Taxa.class, 1, Integer.MAX_VALUE)
            }),
            new ElementRule(RECOMBINANT2,
                    new XMLSyntaxRule[]{
                            new ElementRule(Taxa.class, 1, Integer.MAX_VALUE)
            }),
    };
}
