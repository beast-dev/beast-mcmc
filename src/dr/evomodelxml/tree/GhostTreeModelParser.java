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
import dr.evolution.util.Taxa;
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
public class GhostTreeModelParser extends AbstractXMLObjectParser {

    public static final String CORPOREAL_TAXA = "corporealTaxa";
    public static final String ROOT_HEIGHT = "rootHeight";
    public static final String LEAF_HEIGHT = "leafHeight";
    public static final String LEAF_TRAIT = "leafTrait";

    public static final String NODE_HEIGHTS = "nodeHeights";
    public static final String NODE_RATES = "nodeRates";

    public static final String ROOT_NODE = "rootNode";
    public static final String INTERNAL_NODES = "internalNodes";
    public static final String LEAF_NODES = "leafNodes";

    public static final String TAXON = "taxon";
    public static final String NAME = "name";

    public String getParserName() {
        return TreeModel.TREE_MODEL;
    }

    /**
     * @return a tree object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree) xo.getChild(Tree.class);

        TaxonList corporealTaxa = (TaxonList) xo.getElementFirstChild(CORPOREAL_TAXA);

        TreeModel subTreeModel = new TreeModel(TreeModel.TREE_MODEL);

        GhostTreeModel treeModel = new GhostTreeModel(xo.getId(), tree, subTreeModel);

        Logger.getLogger("dr.evomodel").info("\nCreating the tree model, '" + xo.getId() + "'");

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof XMLObject) {

                XMLObject cxo = (XMLObject) xo.getChild(i);

                switch (cxo.getName()) {
                    case ROOT_HEIGHT:

                        ParameterParser.replaceParameter(cxo, treeModel.getRootHeightParameter());

                        break;
                    case LEAF_HEIGHT: {

                        String taxonName;
                        if (cxo.hasAttribute(TAXON)) {
                            taxonName = cxo.getStringAttribute(TAXON);
                        } else {
                            throw new XMLParseException("taxa element missing from leafHeight element in treeModel element");
                        }

                        int index = treeModel.getTaxonIndex(taxonName);
                        if (index == -1) {
                            throw new XMLParseException("taxon " + taxonName + " not found for leafHeight element in treeModel element");
                        }
                        NodeRef node = treeModel.getExternalNode(index);

                        Parameter newParameter = treeModel.getLeafHeightParameter(node);

                        ParameterParser.replaceParameter(cxo, newParameter);

                        Taxon taxon = treeModel.getTaxon(index);

                        setUncertaintyBounds(newParameter, taxon);

                        break;
                    }
                    case NODE_HEIGHTS:

                        boolean rootNode = cxo.getAttribute(ROOT_NODE, false);
                        boolean internalNodes = cxo.getAttribute(INTERNAL_NODES, false);
                        boolean leafNodes = cxo.getAttribute(LEAF_NODES, false);

                        if (!rootNode && !internalNodes && !leafNodes) {
                            throw new XMLParseException("one or more of root, internal or leaf nodes must be selected for the nodeHeights element");
                        }

                        ParameterParser.replaceParameter(cxo, treeModel.createNodeHeightsParameter(rootNode, internalNodes, leafNodes));

                        break;
                    case LEAF_TRAIT: {

                        String name = cxo.getAttribute(NAME, "trait");

                        String taxonName;
                        if (cxo.hasAttribute(TAXON)) {
                            taxonName = cxo.getStringAttribute(TAXON);
                        } else {
                            throw new XMLParseException("taxa element missing from leafTrait element in treeModel element");
                        }

                        int index = treeModel.getTaxonIndex(taxonName);
                        if (index == -1) {
                            throw new XMLParseException("taxon '" + taxonName + "' not found for leafTrait element in treeModel element");
                        }
                        NodeRef node = treeModel.getExternalNode(index);

                        Parameter parameter = treeModel.getNodeTraitParameter(node, name);

                        if (parameter == null)
                            throw new XMLParseException("trait '" + name + "' not found for leafTrait (taxon, " + taxonName + ") element in treeModel element");

                        ParameterParser.replaceParameter(cxo, parameter);

                        break;
                    }
                    default:
                        throw new XMLParseException("illegal child element in " + getParserName() + ": " + cxo.getName());
                }

            } else if (xo.getChild(i) instanceof Tree) {
                // do nothing - already handled
            } else {
                throw new XMLParseException("illegal child element in  " + getParserName() + ": " + xo.getChildName(i) + " " + xo.getChild(i));
            }
        }

//        Logger.getLogger("dr.evomodel").info("  initial tree topology = " + TreeUtils.uniqueNewick(treeModel, treeModel.getRoot()));
        Logger.getLogger("dr.evomodel").info("  taxon count = " + treeModel.getExternalNodeCount());
        Logger.getLogger("dr.evomodel").info("  tree height = " + treeModel.getNodeHeight(treeModel.getRoot()));

        return treeModel;
    }

    private void setUncertaintyBounds(Parameter newParameter, Taxon taxon) {
        Date date = taxon.getDate();
        if (date != null) {
            double uncertainty = date.getUncertainty();
            if (uncertainty > 0.0) {
                // taxon date not specified to exact value so add appropriate bounds
                double upper = Taxon.getHeightFromDate(date);
                double lower = Taxon.getHeightFromDate(date);
                if (date.isBackwards()) {
                    upper += uncertainty;
                } else {
                    lower -= uncertainty;
                }

                // set the bounds for the given precision
                newParameter.addBounds(new Parameter.DefaultBounds(upper, lower, 1));

                // set the initial value to be mid-point
                newParameter.setParameterValue(0, (upper + lower) / 2);
            }
        }
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
            new ElementRule(CORPOREAL_TAXA, TaxonList.class, "A list of taxa which are the non-ghost tips", false),
            new ElementRule(ROOT_HEIGHT, Parameter.class, "A parameter definition with id only (cannot be a reference!)", false),
            new ElementRule(NODE_HEIGHTS,
                    new XMLSyntaxRule[]{
                            AttributeRule.newBooleanRule(ROOT_NODE, true, "If true the root height is included in the parameter"),
                            AttributeRule.newBooleanRule(INTERNAL_NODES, true, "If true the internal node heights (minus the root) are included in the parameter"),
                            new ElementRule(Parameter.class, "A parameter definition with id only (cannot be a reference!)")
                    }, 1, Integer.MAX_VALUE),
            new ElementRule(LEAF_HEIGHT,
                    new XMLSyntaxRule[]{
                            AttributeRule.newStringRule(TAXON, false, "The name of the taxon for the leaf"),
                            new ElementRule(Parameter.class, "A parameter definition with id only (cannot be a reference!)")
                    }, 0, Integer.MAX_VALUE)
    };
}
