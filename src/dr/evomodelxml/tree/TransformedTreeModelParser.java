/*
 * TransformedTreeModelParser.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evomodel.tree.TransformedTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */
public class TransformedTreeModelParser extends AbstractXMLObjectParser {

    public static final String STAR_TREE_MODEL = "transformedTreeModel";

    public String getParserName() {
        return STAR_TREE_MODEL;
    }

    /**
     * @return a tree object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        Parameter scale = (Parameter) xo.getChild(Parameter.class);

        String id = tree.getId();
        if (!xo.hasId()) {
//            System.err.println("No check!");
            id = "transformed." + id;
        } else {
//            System.err.println("Why am I here?");
            id = xo.getId();
        }
        Logger.getLogger("dr.evomodel").info("Creating a transformed tree model, '" + id + "'");

        return new TransformedTreeModel(id, tree, scale);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a transformed model of the tree.";
    }

    public Class getReturnType() {
        return TransformedTreeModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules =
            new XMLSyntaxRule[]{
                    new ElementRule(TreeModel.class),
                    new ElementRule(Parameter.class),
            };
}
