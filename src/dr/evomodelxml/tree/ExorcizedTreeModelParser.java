/*
 * CorporealTreeModelParser.java
 *
 * Copyright © 2002-2026, the BEAST Development Team.
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
 */

package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.ExorcizedTreeModel;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 */
public class ExorcizedTreeModelParser extends AbstractXMLObjectParser {

    public static final String EXORCIZED_TREE_MODEL = "exorcizedTreeModel";
    public static final String CORPOREAL_TAXA = "corporealTaxa";
    public static final String HAUNTED_TREE = "hauntedTree";

    public String getParserName() {
        return ExorcizedTreeModelParser.EXORCIZED_TREE_MODEL;
    }

    /**
     * @return a tree object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree) xo.getChild(Tree.class);

        TaxonList corporealTaxa = (TaxonList) xo.getElementFirstChild(CORPOREAL_TAXA);
        Tree ghostlyTree = (Tree) xo.getElementFirstChild(HAUNTED_TREE);

        ExorcizedTreeModel treeModel = new ExorcizedTreeModel(xo.getId(), tree, corporealTaxa);

        Logger.getLogger("dr.evomodel").info("\nCreating the corporeal tree model, '" + xo.getId() + "'" +
                "\n\nwith " + corporealTaxa.getTaxonCount() + " taxa.");

        return treeModel;
    }


    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This a tree with ghost lineages (branches without sequence data) and extracts the non-ghost subtree.";
    }

    public Class getReturnType() {
        return ExorcizedTreeModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(HAUNTED_TREE, Tree.class, "The tree containing ghost lineages", false),
            new ElementRule(CORPOREAL_TAXA, TaxonList.class, "A list of taxa which are the non-ghost lineages", false),
    };
}
