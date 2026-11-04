/*
 * IGCoalescentLikelihoodParser.java
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

import dr.evolution.coalescent.IntervalList;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.IGCoalescentLikelihood;
import dr.evomodel.coalescent.TreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Parser for the coalescent likelihood with a constant population size integrated out
 * under an inverse-Gamma(alpha, beta) prior.
 */
public class IGCoalescentLikelihoodParser extends AbstractXMLObjectParser {

    public static final String IGCOALESCENT_LIKELIHOOD = "coalescentLikelihoodIG";
    public static final String POPULATION_TREE = "populationTree";

    public static final String ALPHA = "alpha";
    public static final String BETA = "beta";

    public static final String INCLUDE = "include";
    public static final String EXCLUDE = "exclude";

    public String getParserName() {
        return IGCOALESCENT_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double alpha = xo.getDoubleAttribute(ALPHA);
        final double beta = xo.getDoubleAttribute(BETA);

        final XMLObject treeElement = xo.getChild(POPULATION_TREE);
        final TreeModel treeModel = (TreeModel) treeElement.getChild(TreeModel.class);

        TaxonList includeSubtree = null;
        if (xo.hasChildNamed(INCLUDE)) {
            includeSubtree = (TaxonList) xo.getElementFirstChild(INCLUDE);
        }

        List<TaxonList> excludeSubtrees = new ArrayList<TaxonList>();
        if (xo.hasChildNamed(EXCLUDE)) {
            XMLObject cxo = xo.getChild(EXCLUDE);
            for (int i = 0; i < cxo.getChildCount(); i++) {
                excludeSubtrees.add((TaxonList) cxo.getChild(i));
            }
        }

        final IntervalList intervalList;
        try {
            intervalList = new TreeIntervals(treeModel, includeSubtree, excludeSubtrees, false);
        } catch (TreeUtils.MissingTaxonException mte) {
            throw new XMLParseException("treeModel missing a taxon from taxon list in " + getParserName() + " element");
        }

        return new IGCoalescentLikelihood(intervalList, alpha, beta);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the coalescent likelihood of a tree with the (constant) " +
                "population size integrated out under an inverse-Gamma prior.";
    }

    public Class getReturnType() {
        return IGCoalescentLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(ALPHA),
            AttributeRule.newDoubleRule(BETA),

            new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
                    new ElementRule(TreeModel.class)
            }, "The tree to compute the likelihood for"),

            new ElementRule(INCLUDE, new XMLSyntaxRule[]{
                    new ElementRule(Taxa.class)
            }, "An optional subset of taxa on which to calculate the likelihood (should be monophyletic)", true),

            new ElementRule(EXCLUDE, new XMLSyntaxRule[]{
                    new ElementRule(Taxa.class, 1, Integer.MAX_VALUE)
            }, "One or more subsets of taxa which should be excluded from calculate the likelihood (should be monophyletic)", true)
    };
}
