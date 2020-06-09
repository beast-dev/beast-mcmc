/*
 * CoalescentLikelihoodParser.java
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

package dr.evomodelxml.coalescent;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.*;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * MultiTreeIntervalsParser
 *
 * @author Andrew Rambaut
 */
public class MultiTreeIntervalsParser extends AbstractXMLObjectParser {

    public static final String MULTI_TREE_INTERVALS = "multiTreeIntervals";
    public static final String TREES = "trees";
    public static final String SINGLETONS = "singletons";
    public static final String INCLUDE_STEMS = "includeStems";
    public static final String CUTOFF = "cutoff";

    public String getParserName() {
        return MULTI_TREE_INTERVALS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(TREES);

        List<Tree> trees = new ArrayList<Tree>(cxo.getAllChildren(Tree.class));

        Taxa singletonTaxa = null;
        if(xo.hasChildNamed(SINGLETONS)){
            singletonTaxa = (Taxa)xo.getElementFirstChild(SINGLETONS);
        }


        boolean includeStems = xo.getBooleanAttribute(INCLUDE_STEMS);

        double cutoffTime = 0.0;
        if (includeStems) {
            if (!xo.hasAttribute(CUTOFF)) {
                throw new XMLParseException("MultiTreeIntervals needs a cutoff time if it is to include stems");
            }
            cutoffTime = xo.getDoubleAttribute(CUTOFF);
        }

        return new MultiTreeIntervals(trees, singletonTaxa, includeStems, cutoffTime);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element stores a set of coalescent intervals over multiple trees.";
    }

    public Class getReturnType() {
        return MultiTreeIntervals.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(INCLUDE_STEMS, false),
            AttributeRule.newDoubleRule(CUTOFF, true),

            new ElementRule(TREES, new XMLSyntaxRule[] {
                    new ElementRule(Tree.class, 1, Integer.MAX_VALUE)
            }, "Tree(s) to compute intervals for for", false),

            new ElementRule(SINGLETONS, new XMLSyntaxRule[]{
                    new ElementRule(Taxa.class)
            }, "An optional set of taxa which represent singleton trees", true)
    };
}
