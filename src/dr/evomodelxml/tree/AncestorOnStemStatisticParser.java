/*
 * MonophylyStatisticParser.java
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
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.continuous.AncestralTaxonInTree;
import dr.evomodel.tree.AncestorOnStemStatistic;
import dr.evomodel.tree.AncestralTraitTreeModel;
import dr.evomodel.tree.MonophylyStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

import static dr.evomodelxml.tree.MonophylyStatisticParser.parseTaxonListOrTaxa;

/**
 */
public class AncestorOnStemStatisticParser extends AbstractXMLObjectParser {

    public static final String ANCESTOR_STATISTIC = "ancestorOnStemStatistic";
    private static final String ANCESTOR = "ancestor";
    private static final String MRCA = "mrca";

    public String getParserName() {
        return ANCESTOR_STATISTIC;
    }


    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());

        AncestralTraitTreeModel tree = (AncestralTraitTreeModel)
                xo.getChild(AncestralTraitTreeModel.class);

        TaxonList mrca = parseTaxonListOrTaxa(xo.getChild(MRCA));

        Taxon ancestor = (Taxon) xo.getChild(ANCESTOR).getChild(Taxon.class);

        try {
            return new AncestorOnStemStatistic(name, tree, ancestor, mrca);
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
            new ElementRule(AncestralTraitTreeModel.class),
            new ElementRule(MRCA, new XMLSyntaxRule[]{
                    new XORRule(
                            new ElementRule(Taxon.class, 1, Integer.MAX_VALUE),
                            new ElementRule(Taxa.class)
                    )
            }),
    };

}
