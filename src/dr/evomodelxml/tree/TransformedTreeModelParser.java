/*
 * TransformedTreeModelParser.java
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

import dr.evolution.tree.MutableTreeModel;
import dr.evomodel.tree.*;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */
public class TransformedTreeModelParser extends AbstractXMLObjectParser {

    public static final String TRANSFORMED_TREE_MODEL = "transformedTreeModel";
    public static final String VERSION = "version";

    public String getParserName() {
        return TRANSFORMED_TREE_MODEL;
    }

    /**
     * @return a tree object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        MutableTreeModel tree = (MutableTreeModel) xo.getChild(MutableTreeModel.class);
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

        TreeTransform transform;

        String version = xo.getAttribute(VERSION, "generic");

        if (version.compareTo("new") == 0) {
            transform = new ProgressiveScalarTreeTransform(scale);
        } else if (version.compareTo("branch") == 0) {
            transform = new ProgressiveScalarTreeTransform(tree, scale);
        } else if (version.compareTo("ou") == 0) {
            transform = new OuScalarTreeTransform(scale);
        } else {
            transform = new SingleScalarTreeTransform(scale);
        }

        return new TransformedTreeModel(id, tree, transform);
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
                    new ElementRule(MutableTreeModel.class),
                    new ElementRule(Parameter.class),
                    AttributeRule.newStringRule(VERSION, true),
            };
}
