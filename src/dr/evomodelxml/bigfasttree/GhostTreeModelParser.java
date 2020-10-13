/*
 * TreeModelParser.java
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

package dr.evomodelxml.bigfasttree;

import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.bigfasttree.GhostTreeModel;

import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 */
public class GhostTreeModelParser extends AbstractXMLObjectParser {

    public static final String GHOST_TAXA = "ghostTaxa";
    public String getParserName() {
        return GhostTreeModel.GHOST_TREE_MODEL;
    }

    /**
     * @return a tree object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree) xo.getChild(Tree.class);

        TaxonList ghostTaxa = (TaxonList) xo.getElementFirstChild(GHOST_TAXA);

        GhostTreeModel treeModel = new GhostTreeModel(xo.getId(), tree, ghostTaxa);

        Logger.getLogger("dr.evomodel").info("\nCreating the tree model, '" + xo.getId() + "'" +
                "\n\nwith " + ghostTaxa.getTaxonCount() + " ghost lineages.");

        Logger.getLogger("dr.evomodel").info("  taxon count = " + treeModel.getExternalNodeCount());
        Logger.getLogger("dr.evomodel").info("  tree height = " + treeModel.getNodeHeight(treeModel.getRoot()));

        return treeModel;
    }


    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a model of the tree with ghost lineages (branches without sequence data).";
    }

    public Class getReturnType() {
        return GhostTreeModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Tree.class),
            new ElementRule(GHOST_TAXA, TaxonList.class, "A list of taxa which are the ghost lineages", false),
    };
}
