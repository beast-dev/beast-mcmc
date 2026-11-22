/*
 * TreeIntervalsParser.java
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

package dr.evomodelxml.coalescent.intervals;

import dr.evolution.coalescent.TreeIntervalList;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.intervals.FastTreeIntervals;
import dr.evomodel.coalescent.intervals.TreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class TreeIntervalsParser extends AbstractXMLObjectParser{

    public static final String TREE_INTERVALS = "treeIntervals";
    public static final String INCLUDE = "include";
    public static final String EXCLUDE = "exclude";

    public static final boolean USE_FAST_INTERVALS = true;

    public String getParserName() {
        return TREE_INTERVALS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree) xo.getChild(Tree.class);

        boolean useFastTreeIntervals = tree instanceof TreeModel;

        TaxonList includeSubtree = null;

        if (xo.hasChildNamed(INCLUDE)) {
            includeSubtree = (TaxonList) xo.getElementFirstChild(INCLUDE);
            // FastTreeIntervals doesn't support include/exclude
            useFastTreeIntervals = false;
        }

        List<TaxonList> excludeSubtrees = new ArrayList<>();

        if (xo.hasChildNamed(EXCLUDE)) {
            XMLObject cxo = xo.getChild(EXCLUDE);
            for (int i = 0; i < cxo.getChildCount(); i++) {
                excludeSubtrees.add((TaxonList) cxo.getChild(i));
            }
            // FastTreeIntervals doesn't support include/exclude
            useFastTreeIntervals = false;
        }

        try {
            if (useFastTreeIntervals) {
                return new FastTreeIntervals((TreeModel)tree);
            } else {
                TreeIntervals intervals = new TreeIntervals(tree, includeSubtree, excludeSubtrees, false);
                if (!intervals.isMonophyly()) {
                    throw new XMLParseException("The included or excluded clades in TreeLineages with id, " + xo.getId() + ", are not monophyletic..");
                }
                return intervals;
            }

        } catch (TreeUtils.MissingTaxonException mte) {
            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Parser for TreeIntervals.";
    }

    public Class getReturnType() {
        return TreeIntervalList.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(INCLUDE, new XMLSyntaxRule[]{
                    new ElementRule(Taxa.class)
            }, "An optional subset of taxa on which to calculate the likelihood (should be monophyletic)", true),

            new ElementRule(EXCLUDE, new XMLSyntaxRule[]{
                    new ElementRule(Taxa.class, 1, Integer.MAX_VALUE)
            }, "One or more subsets of taxa which should be excluded from calculate the likelihood (should be monophyletic)", true)
    };

}