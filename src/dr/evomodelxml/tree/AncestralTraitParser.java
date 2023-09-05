/*
 * AncestralTraitParser.java
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

import dr.evolution.tree.*;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.AncestralTrait;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 */
public class AncestralTraitParser extends AbstractXMLObjectParser {

    public static final String ANCESTRAL_TRAIT = "ancestralTrait";
    public static final String ANCESTRAL_STATE = "ancestralState";
    public static final String NAME = "name";
    public static final String MRCA = "mrca";
    public static final String TRAIT_NAME = "traitName";
    public static final String STATES = "states";

    public String getParserName() {
        return ANCESTRAL_TRAIT;
    }

    public String[] getParserNames() {
        // provide a synonym to maintain backwards compatibility
        return new String[]{
                getParserName(), ANCESTRAL_STATE
        };
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, STATES);
        String name = xo.getAttribute(NAME, traitName);
        Tree tree = (Tree) xo.getChild(Tree.class);
        TreeTraitProvider treeTraitProvider = (TreeTraitProvider) xo.getChild(TreeTraitProvider.class);

        TaxonList taxa = null;
        if (xo.hasChildNamed(MRCA)) {
            taxa = (TaxonList) xo.getElementFirstChild(MRCA);
        }

        TreeTrait trait = treeTraitProvider.getTreeTrait(traitName);
        if (trait == null) {
            throw new XMLParseException("A trait called, " + traitName + ", was not available from the TreeTraitProvider supplied to " + getParserName() + (xo.hasId() ? ", with ID " + xo.getId() : ""));
        }
        try {
            return new AncestralTrait(name, trait, tree, taxa);
        } catch (TreeUtils.MissingTaxonException mte) {
            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that has as its value the height of the most recent common ancestor of a set of taxa in a given tree";
    }

    public Class getReturnType() {
        return AncestralTrait.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(NAME, "A name for this statistic primarily for the purposes of logging", true),
            new StringAttributeRule(TRAIT_NAME, "The name of the trait to log", true),
            new ElementRule(Tree.class),
            new ElementRule(TreeTraitProvider.class),
            new ElementRule(MRCA, new XMLSyntaxRule[]{new ElementRule(Taxa.class)},  "The MRCA to reconstruct the trait at (default root node)", true)
    };
}
