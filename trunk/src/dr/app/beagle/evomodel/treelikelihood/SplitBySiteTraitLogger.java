/*
 * SplitBySiteTraitLogger.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.evomodel.treelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AncestralStateTreeLikelihood;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Takes a TreeTraitProvider that returns multiple-site traits and provides a new Trait for logging by site
 *
 * @author Marc A. Suchard
 */
public class SplitBySiteTraitLogger extends TreeTraitProvider.Helper implements Citable {

    public static final String TRAIT_LOGGER = "splitTraitBySite";
    public static final String TRAIT_NAME = "traitName";
    public static final String SCALE = "scaleByBranchLength";

    public SplitBySiteTraitLogger(final AncestralStateBeagleTreeLikelihood treeLikelihood, String traitName, boolean scale) throws XMLParseException {
        TreeTrait trait = treeLikelihood.getTreeTrait(traitName);
        if (trait == null) {
            throw new XMLParseException("TreeTraitProvider does not provide trait named '" + traitName + ".");
        }

        TreeModel tree = treeLikelihood.getTreeModel();
        int length;
        Object obj = trait.getTrait(tree, tree.getNode(0));

        boolean isDoubleArray = false;
        boolean isIntegerArray = false;

        if (obj instanceof double[]) {
            length = ((double[]) obj).length;
            isDoubleArray = true;
        } else if (obj instanceof int[]) {
            length = ((int[])obj).length;
            isIntegerArray = true;
        } else {
            throw new XMLParseException("Unknown trait type to split");
        }

        TreeTrait[] partitionedTraits = new TreeTrait[length];
        if (isDoubleArray) {
            for (int i = 0; i < length; i++) {
                if (scale) {
                    partitionedTraits[i] = new TreeTrait.PickEntryDAndScale(trait, i);
                } else {
                    partitionedTraits[i] = new TreeTrait.PickEntryD(trait, i);
                }
            }
        } else if (isIntegerArray) {
            for (int i = 0; i < length; ++i) {
                if (traitName.compareTo(AncestralStateTreeLikelihood.STATES_KEY) == 0) {
                    partitionedTraits[i] = new TreeTrait.PickEntryI(trait, i) {
                        public String getTraitString(Tree tree, NodeRef node) {
                            int[] state = new int[1];
                            state[0] = getTrait(tree, node);
                            return treeLikelihood.formattedState(state);
                        }
                    };
                } else {
                    partitionedTraits[i] = new TreeTrait.PickEntryI(trait, i);
                }
            }
        }
        addTraits(partitionedTraits);        

        Logger.getLogger("dr.app.beagle").info("\tConstructing a split logger with " + length + " partitions;  please cite:\n"
                + Citable.Utils.getCitationString(this));

    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TRAIT_LOGGER;
        }

        /**
         * @return an object based on the XML element it was passed.
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String traitName = xo.getStringAttribute(TRAIT_NAME);
            AncestralStateBeagleTreeLikelihood tree = (AncestralStateBeagleTreeLikelihood) xo.getChild(AncestralStateBeagleTreeLikelihood.class);
            boolean scale = xo.getAttribute(SCALE, false);

            return new SplitBySiteTraitLogger(tree, traitName, scale);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(AncestralStateBeagleTreeLikelihood.class, "The tree which is to be logged"),
                AttributeRule.newStringRule(TRAIT_NAME),
                AttributeRule.newBooleanRule(SCALE, true),
        };

        public String getParserDescription() {
            return null;
        }

        public String getExample() {
            return null;
        }

        public Class getReturnType() {
            return TreeTraitProvider.class;
        }
    };

    /**
     * @return a list of citations associated with this object
     */
    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(
                CommonCitations.SUCHARD_2012
        );
        return citations;
    }
}
