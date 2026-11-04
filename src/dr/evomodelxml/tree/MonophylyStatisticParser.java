/*
 * MonophylyStatisticParser.java
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
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.MonophylyStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class MonophylyStatisticParser extends AbstractXMLObjectParser {

    public static final String MONOPHYLY_STATISTIC = "monophylyStatistic";
    public static final String MRCA = "mrca";
    public static final String IGNORE = "ignore";
    public static final String INVERSE = "inverse";

    public String getParserName() {
        return MONOPHYLY_STATISTIC;
    }


    public static TaxonList parseTaxonListOrTaxa(XMLObject cxo) {
        TaxonList taxa = (TaxonList) cxo.getChild(TaxonList.class);
        if (taxa == null) {
            Taxa taxa1 = new Taxa();
            for (int i = 0; i < cxo.getChildCount(); i++) {
                Object ccxo = cxo.getChild(i);
                if (ccxo instanceof Taxon) {
                    taxa1.addTaxon((Taxon) ccxo);
                }
            }
            taxa = taxa1;
        }
        return taxa;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());

        Boolean inverse = xo.getAttribute(INVERSE, false);

        Tree tree = (Tree) xo.getChild(Tree.class);

        TaxonList taxa = parseTaxonListOrTaxa(xo.getChild(MRCA));

        TaxonList ignore = null;
        if (xo.hasChildNamed(IGNORE)) {
            ignore = parseTaxonListOrTaxa(xo.getChild(IGNORE));
        }

        try {
            return new MonophylyStatistic(name, tree, taxa, ignore, inverse);
        } catch (TreeUtils.MissingTaxonException mte) {
            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that returns true if a given set of taxa are monophyletic for a given tree";
    }

    public Class getReturnType() {
        return MonophylyStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new StringAttributeRule(Statistic.NAME, "A name for this statistic for the purpose of logging", true),
            AttributeRule.newBooleanRule(INVERSE, true, "inverse, returns 0/false when monophyletic and 1/true when not monophyletic"),
            // Any tree will do, no need to insist on a Tree Model
            new ElementRule(Tree.class),
            new ElementRule(MRCA, new XMLSyntaxRule[]{
                    new XORRule(
                            new ElementRule(Taxon.class, 1, Integer.MAX_VALUE),
                            new ElementRule(Taxa.class)
                    )
            }),
            new ElementRule(IGNORE, new XMLSyntaxRule[]{
                    new XORRule(
                            new ElementRule(Taxon.class, 1, Integer.MAX_VALUE),
                            new ElementRule(Taxa.class)
                    )
            }, "An optional list of taxa to ignore from the test of monophyly", true)
    };

}
