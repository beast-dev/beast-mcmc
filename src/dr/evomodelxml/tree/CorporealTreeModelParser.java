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

package dr.evomodelxml.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.GhostTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 */
public class CorporealTreeModelParser extends AbstractXMLObjectParser {

    public static final String CORPOREAL_TREE_MODEL = "corporealTreeModel";

    public String getParserName() {
        return CORPOREAL_TREE_MODEL;
    }

    /**
     * @return a tree object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        GhostTreeModel ghostTreeModel = (GhostTreeModel) xo.getChild(GhostTreeModel.class);

        return ghostTreeModel.getCorporealTreeModel();
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element extracts the corporeal TreeModel from a GhostTreeModel.";
    }

    public Class getReturnType() {
        return TreeModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(TreeModel.class),
    };
}
