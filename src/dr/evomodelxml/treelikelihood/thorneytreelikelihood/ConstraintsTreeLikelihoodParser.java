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

package dr.evomodelxml.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.treelikelihood.thorneytreelikelihood.ConstraintsTreeLikelihood;
import dr.xml.*;

/**
 */
public class ConstraintsTreeLikelihoodParser extends AbstractXMLObjectParser {

    public static final String CONSTRAINTS_TREE_LIKELIHOOD = "constraintsTreeLikelihood";
    public static final String CONSTRAINTS = "constraints";
    public static final String TREE = "tree";

    public String getParserName() {
        return CONSTRAINTS_TREE_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree constraintsTree = (Tree) xo.getElementFirstChild(CONSTRAINTS);
        Tree targetTree = (Tree) xo.getChild(Tree.class);

        try {
            return new ConstraintsTreeLikelihood(CONSTRAINTS_TREE_LIKELIHOOD, targetTree, constraintsTree);
        } catch (TreeUtils.MissingTaxonException e) {
            throw new XMLParseException("Target tree and constraints tree do not contain the same taxa");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This binary likelihood that assesses whether a resolved tree is compatible with a constraints tree.";
    }

    public Class getReturnType() {
        return ConstraintsTreeLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(CONSTRAINTS, new XMLSyntaxRule[]{
                    new ElementRule(Tree.class)
            }, "Less resolved tree to provide constraints"),

            new ElementRule(Tree.class,"Tree to assess for compatibility with constraints")

    };
}
