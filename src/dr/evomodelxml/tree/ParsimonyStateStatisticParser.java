/*
 * ParsimonyStateStatisticParser.java
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
import dr.evomodel.tree.ParsimonyStateStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class ParsimonyStateStatisticParser extends AbstractXMLObjectParser {

    public static final String PARSIMONY_STATE_STATISTIC = "parsimonyStateStatistic";
    public static final String STATE = "state";
    public static final String MRCA = "mrca";

    public String getParserName() {
        return PARSIMONY_STATE_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());

        Tree tree = (Tree) xo.getChild(Tree.class);
        XMLObject cxo = xo.getChild(STATE);
        TaxonList stateTaxa = (TaxonList) cxo.getChild(TaxonList.class);

        cxo = xo.getChild(MRCA);
        TaxonList mrcaTaxa = (TaxonList) cxo.getChild(TaxonList.class);

        try {
            return new ParsimonyStateStatistic(name, tree, stateTaxa, mrcaTaxa);
        } catch (TreeUtils.MissingTaxonException mte) {
            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that has as its value the parsimony state reconstruction of a " +
                "binary state defined by a set of taxa at a given MRCA of a tree";
    }

    public Class getReturnType() {
        return ParsimonyStateStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(Statistic.NAME, "A name for this statistic for the purposes of logging", true),
            new ElementRule(TreeModel.class),
            new ElementRule(STATE,
                    new XMLSyntaxRule[]{new ElementRule(Taxa.class)}),
            new ElementRule(MRCA,
                    new XMLSyntaxRule[]{new ElementRule(Taxa.class)})
    };

}
