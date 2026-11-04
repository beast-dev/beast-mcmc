/*
 * SpeciationLikelihoodParser.java
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

package dr.evomodelxml.birthdeath;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.birthdeath.BirthDeathLikelihood;
import dr.evomodel.birthdeath.BirthDeathModel;
import dr.evomodel.birthdeath.EfficientBirthDeathLikelihood;
import dr.xml.*;

import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

/**
 * Parser for the likelihood of a tree given a birth death model.
 */
public class BirthDeathLikelihoodParser extends AbstractXMLObjectParser {

    public static final String BIRTH_DEATH_LIKELIHOOD = "birthDeathLikelihood";
    public static final String MODEL = "model";
    public static final String TREE = "tree";
    public static final String INCLUDE = "include";
    public static final String EXCLUDE = "exclude";

    public static final String POINT = "point";

    private static final String USE_NEW_LOOP = "useNewLoop";

    public String getParserName() {
        return BIRTH_DEATH_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(MODEL);
        final BirthDeathModel specModel = (BirthDeathModel) cxo.getChild(BirthDeathModel.class);

        cxo = xo.getChild(TREE);
        final Tree tree = (Tree) cxo.getChild(Tree.class);

        Set<Taxon> excludeTaxa = null;

        if (xo.hasChildNamed(INCLUDE)) {
            excludeTaxa = new HashSet<Taxon>();
            for (int i = 0; i < tree.getTaxonCount(); i++) {
                excludeTaxa.add(tree.getTaxon(i));
            }

            cxo = xo.getChild(INCLUDE);
            for (int i = 0; i < cxo.getChildCount(); i++) {
                TaxonList taxonList = (TaxonList) cxo.getChild(i);
                for (int j = 0; j < taxonList.getTaxonCount(); j++) {
                    excludeTaxa.remove(taxonList.getTaxon(j));
                }
            }
        }

        if (xo.hasChildNamed(EXCLUDE)) {
            excludeTaxa = new HashSet<Taxon>();
            cxo = xo.getChild(EXCLUDE);
            for (int i = 0; i < cxo.getChildCount(); i++) {
                TaxonList taxonList = (TaxonList) cxo.getChild(i);
                for (int j = 0; j < taxonList.getTaxonCount(); j++) {
                    excludeTaxa.add(taxonList.getTaxon(j));
                }
            }
        }
        if (excludeTaxa != null) {
            Logger.getLogger("dr.evomodel").info("Speciation model excluding " + excludeTaxa.size() + " taxa from prior - " +
                    (tree.getTaxonCount() - excludeTaxa.size()) + " taxa remaining.");
        }

        boolean useNewLoop = xo.getAttribute(USE_NEW_LOOP, false);

        if (useNewLoop) {
            return new EfficientBirthDeathLikelihood(tree, specModel, excludeTaxa, null);
        } else {
            return new BirthDeathLikelihood(tree, specModel, excludeTaxa, null);
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of the tree given the speciation.";
    }

    public Class getReturnType() {
        return BirthDeathLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }


    private final XMLSyntaxRule[] rules = {
            new ElementRule(MODEL, new XMLSyntaxRule[]{
                    new ElementRule(BirthDeathModel.class)
            }),
            new ElementRule(TREE, new XMLSyntaxRule[]{
                    new ElementRule(Tree.class)
            }),

            new ElementRule(INCLUDE, new XMLSyntaxRule[]{
                    new ElementRule(Taxa.class, 1, Integer.MAX_VALUE)
            }, "One or more subsets of taxa which should be included from calculate the likelihood (the remaining taxa are excluded)", true),

            new ElementRule(EXCLUDE, new XMLSyntaxRule[]{
                    new ElementRule(Taxa.class, 1, Integer.MAX_VALUE)
            }, "One or more subsets of taxa which should be excluded from calculate the likelihood (which is calculated on the remaining subtree)", true),

            AttributeRule.newBooleanRule(USE_NEW_LOOP, true),
    };

}
