/*
 * ConstrainedTreeModelParser.java
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

package dr.evomodelxml.bigfasttree.thorney;

import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

import java.util.logging.Logger;

public class ConstrainedTreeModelParser extends AbstractXMLObjectParser {

    public static final String CONSTRAINED_TREE_MODEL = "constrainedTreeModel";
    public static final String CONSTRAINTS_TREE = "constraintsTree";

    public String getParserName() {
        return CONSTRAINED_TREE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree) xo.getChild(Tree.class);
        Tree constraintsTree = (Tree) xo.getElementFirstChild(CONSTRAINTS_TREE);
        Logger.getLogger("dr.evomodel").info("\nCreating the constrained tree model based on big fast tree model, '" + xo.getId() + "'");
        ConstrainedTreeModel treeModel =  new ConstrainedTreeModel(xo.getId(),tree, constraintsTree);
        Logger.getLogger("dr.evomodel").info("  taxon count = " + treeModel.getExternalNodeCount());
        Logger.getLogger("dr.evomodel").info("  tree height = " + treeModel.getNodeHeight(treeModel.getRoot()));
        return treeModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the a tree model that with defined clades that may not be broken.";
    }

    public Class getReturnType() {
        return TreeModel.class;
    }

    public static final XMLSyntaxRule[] rules = {
            new ElementRule(Tree.class),
            new ElementRule(CONSTRAINTS_TREE, new XMLSyntaxRule[]{
                    new ElementRule(Tree.class)
            }),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
