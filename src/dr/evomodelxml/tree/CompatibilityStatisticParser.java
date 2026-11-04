/*
 * CompatibilityStatisticParser.java
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
import dr.evomodel.tree.CompatibilityStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class CompatibilityStatisticParser extends AbstractXMLObjectParser {

    public static final String COMPATIBILITY_STATISTIC = "compatibilityStatistic";
    public static final String COMPATIBLE_WITH = "compatibleWith";

    public String getParserName() {
        return COMPATIBILITY_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());

        Tree tree1 = (Tree) xo.getChild(Tree.class);

        XMLObject cxo = xo.getChild(COMPATIBLE_WITH);
        Tree tree2 = (Tree) cxo.getChild(Tree.class);

        try {
            return new CompatibilityStatistic(name, tree1, tree2);
        } catch (TreeUtils.MissingTaxonException mte) {
            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was in the source tree but not the constraints tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that returns true if a pair of trees are compatible";
    }

    public Class getReturnType() {
        return CompatibilityStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(Statistic.NAME, "A name for this statistic for the purpose of logging", true),
            new ElementRule(Tree.class),
            new ElementRule(COMPATIBLE_WITH, new XMLSyntaxRule[]{
                    new ElementRule(Tree.class)
            }),
            new ElementRule(Tree.class)
    };

}
