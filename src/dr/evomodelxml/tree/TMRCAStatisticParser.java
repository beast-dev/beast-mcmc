/*
 * TMRCAStatisticParser.java
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
import dr.evomodel.tree.TMRCAStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 *
 * To get the age of the root in absolute time:
 * <tmrcaStatistic id="age" absolute="true">
 *     <tree idref="tree"/>
 * </tmrcaStatistic>
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class TMRCAStatisticParser extends AbstractXMLObjectParser {

    public static final String TMRCA_STATISTIC = "tmrcaStatistic";
    public static final String MRCA = "mrca";
    public static final String ABSOLUTE = "absolute";

    // The tmrcaStatistic will represent that age of the parent node of the MRCA, rather than the MRCA itself
    public static final String PARENT = "forParent";
    public static final String STEM = "includeStem";


    public String getParserName() {
        return TMRCA_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());
        Tree tree = (Tree) xo.getChild(Tree.class);
        TaxonList taxa = null;

        if (xo.hasChildNamed(MRCA)) {
            taxa = (TaxonList) xo.getElementFirstChild(MRCA);
        }
        boolean isAbsolute = xo.getAttribute(ABSOLUTE, false);
        boolean includeStem = false;
        if (xo.hasAttribute(PARENT) && xo.hasAttribute(STEM)) {
             throw new XMLParseException("Please use either " + PARENT + " or " + STEM + "!");
        } else if (xo.hasAttribute(PARENT)) {
             includeStem = xo.getBooleanAttribute(PARENT);
        } else if (xo.hasAttribute(STEM)) {
             includeStem = xo.getBooleanAttribute(STEM);
        }

        try {
            return new TMRCAStatistic(name, tree, taxa, isAbsolute, includeStem);
        } catch (TreeUtils.MissingTaxonException mte) {
            throw new XMLParseException(
                    "Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that has as its value the height of the most recent common ancestor " +
                "of a set of taxa in a given tree. ";
    }

    public Class getReturnType() {
        return TMRCAStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Tree.class),
            new StringAttributeRule("name",
                    "A name for this statistic primarily for the purposes of logging", true),
            AttributeRule.newBooleanRule(ABSOLUTE, true),
            new ElementRule(MRCA,
                    new XMLSyntaxRule[]{new ElementRule(Taxa.class)}, true),
            new OrRule(
                    new XMLSyntaxRule[]{
                            AttributeRule.newBooleanRule(PARENT, true),
                            AttributeRule.newBooleanRule(STEM, true)
                    }
            )
    };

}
