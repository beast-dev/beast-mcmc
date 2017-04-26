/*
 * StarTreeModelParser.java
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
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.StarTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */
public class StarTreeModelParser extends AbstractXMLObjectParser {

    public static final String STAR_TREE_MODEL = "starTreeModel";

    public static final String SHARE_ROOT = "sharedRootHeight";

    public static final String ROOT_HEIGHT = "rootHeight";
    public static final String LEAF_HEIGHT = "leafHeight";
    public static final String LEAF_TRAIT = "leafTrait";

    //    public static final String NODE_HEIGHTS = "nodeHeights";
    public static final String NODE_RATES = "nodeRates";
    public static final String NODE_TRAITS = "nodeTraits";
    public static final String MULTIVARIATE_TRAIT = "traitDimension";
    public static final String INITIAL_VALUE = "initialValue";

    public static final String ROOT_NODE = "rootNode";
    public static final String INTERNAL_NODES = "internalNodes";
    public static final String LEAF_NODES = "leafNodes";
    public static final String FIRE_TREE_EVENTS = "fireTreeEvents";

    public static final String TAXON = "taxon";
    public static final String NAME = "name";

    public StarTreeModelParser() {
        rules = new XMLSyntaxRule[]{
                new ElementRule(Tree.class),
                new XORRule(
                        new ElementRule(ROOT_HEIGHT, Parameter.class, "A parameter definition with id only (cannot be a reference!)", false),
                        new ElementRule(SHARE_ROOT, TreeModel.class)
                ),
                new ElementRule(LEAF_HEIGHT,
                        new XMLSyntaxRule[]{
                                AttributeRule.newStringRule(TAXON, false, "The name of the taxon for the leaf"),
                                new ElementRule(Parameter.class, "A parameter definition with id only (cannot be a reference!)")
                        }, 0, Integer.MAX_VALUE),
                new ElementRule(LEAF_TRAIT,
                        new XMLSyntaxRule[]{
                                AttributeRule.newStringRule(TAXON, false, "The name of the taxon for the leaf"),
                                AttributeRule.newStringRule(NAME, false, "The name of the trait attribute in the taxa"),
                                new ElementRule(Parameter.class, "A parameter definition with id only (cannot be a reference!)")
                        }, 0, Integer.MAX_VALUE)
        };
    }

    public String getParserName() {
        return STAR_TREE_MODEL;
    }

    /**
     * @return a tree object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree) xo.getChild(Tree.class);
        StarTreeModel treeModel = new StarTreeModel(xo.getId(), tree);

        Logger.getLogger("dr.evomodel").info("Creating a star tree model, '" + xo.getId() + "'");

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof XMLObject) {

                XMLObject cxo = (XMLObject) xo.getChild(i);

                if (cxo.getName().equals(ROOT_HEIGHT)) {

                    if (cxo.getRawChild(0) instanceof Reference) {
                        // Co-opt existing parameter
//                        treeModel.setRootHeightParameter((Parameter) cxo.getChild(Parameter.class));
                        throw new XMLParseException("Can not provide idref to a new root height parameter");
                    } else {
                        ParameterParser.replaceParameter(cxo, treeModel.getRootHeightParameter());
                    }

                } else if (cxo.getName().equals(SHARE_ROOT)) {
                    TreeModel sharedRoot = (TreeModel) cxo.getChild(TreeModel.class);
                    treeModel.setSharedRootHeightParameter(sharedRoot);
                } else if (cxo.getName().equals(LEAF_HEIGHT)) {

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
                    ParameterParser.replaceParameter(cxo, treeModel.getLeafHeightParameter(node));

//                } else if (cxo.getName().equals(NODE_HEIGHTS)) {
//
//                    boolean rootNode = cxo.getAttribute(ROOT_NODE, false);
//                    boolean internalNodes = cxo.getAttribute(INTERNAL_NODES, false);
//                    boolean leafNodes = cxo.getAttribute(LEAF_NODES, false);
//
//                    if (!rootNode && !internalNodes && !leafNodes) {
//                        throw new XMLParseException("one or more of root, internal or leaf nodes must be selected for the nodeHeights element");
//                    }
//
//                    ParameterParser.replaceParameter(cxo, treeModel.createNodeHeightsParameter(rootNode, internalNodes, leafNodes));
//
//                } else if (cxo.getName().equals(NODE_RATES)) {
//
//                    boolean rootNode = cxo.getAttribute(ROOT_NODE, false);
//                    boolean internalNodes = cxo.getAttribute(INTERNAL_NODES, false);
//                    boolean leafNodes = cxo.getAttribute(LEAF_NODES, false);
//                    double[] initialValues = null;
//
//                    if (cxo.hasAttribute(INITIAL_VALUE)) {
//                        initialValues = cxo.getDoubleArrayAttribute(INITIAL_VALUE);
//                    }
//
//                    if (!rootNode && !internalNodes && !leafNodes) {
//                        throw new XMLParseException("one or more of root, internal or leaf nodes must be selected for the nodeRates element");
//                    }
//
//                    ParameterParser.replaceParameter(cxo, treeModel.createNodeRatesParameter(initialValues, rootNode, internalNodes, leafNodes));
//
//                } else if (cxo.getName().equals(NODE_TRAITS)) {
//
//                    boolean rootNode = cxo.getAttribute(ROOT_NODE, false);
//                    boolean internalNodes = cxo.getAttribute(INTERNAL_NODES, false);
//                    boolean leafNodes = cxo.getAttribute(LEAF_NODES, false);
//                    boolean fireTreeEvents = cxo.getAttribute(FIRE_TREE_EVENTS, false);
//                    String name = cxo.getAttribute(NAME, "trait");
//                    int dim = cxo.getAttribute(MULTIVARIATE_TRAIT, 1);
//
//                    double[] initialValues = null;
//                    if (cxo.hasAttribute(INITIAL_VALUE)) {
//                        initialValues = cxo.getDoubleArrayAttribute(INITIAL_VALUE);
//                    }
//
//                    if (!rootNode && !internalNodes && !leafNodes) {
//                        throw new XMLParseException("one or more of root, internal or leaf nodes must be selected for the nodeTraits element");
//                    }
//
//                    ParameterParser.replaceParameter(cxo, treeModel.createNodeTraitsParameter(name, dim, initialValues, rootNode, internalNodes, leafNodes, fireTreeEvents));

                } else if (cxo.getName().equals(LEAF_TRAIT)) {

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


                } else {
                    throw new XMLParseException("illegal child element in " + getParserName() + ": " + cxo.getName());
                }

            } else if (xo.getChild(i) instanceof Tree) {
                // do nothing - already handled
            } else {
                throw new XMLParseException("illegal child element in  " + getParserName() + ": " + xo.getChildName(i) + " " + xo.getChild(i));
            }
        }

        Logger.getLogger("dr.evomodel").info("  initial tree topology = " + TreeUtils.uniqueNewick(treeModel, treeModel.getRoot()));
        Logger.getLogger("dr.evomodel").info("  tree height = " + treeModel.getNodeHeight(treeModel.getRoot()));
        return treeModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a model of the tree. The tree model includes and attributes of the nodes " +
                "including the age (or <i>height</i>) and the rate of evolution at each node in the tree.";
    }

    public String getExample() {
        return
                "<!-- the tree model as special sockets for attaching parameters to various aspects of the tree     -->\n" +
                        "<!-- The treeModel below shows the standard setup with a parameter associated with the root height -->\n" +
                        "<!-- a parameter associated with the internal node heights (minus the root height) and             -->\n" +
                        "<!-- a parameter associates with all the internal node heights                                     -->\n" +
                        "<!-- Notice that these parameters are overlapping                                                  -->\n" +
                        "<!-- The parameters are subsequently used in operators to propose changes to the tree node heights -->\n" +
                        "<treeModel id=\"treeModel1\">\n" +
                        "	<tree idref=\"startingTree\"/>\n" +
                        "	<rootHeight>\n" +
                        "		<parameter id=\"treeModel1.rootHeight\"/>\n" +
                        "	</rootHeight>\n" +
                        "	<nodeHeights internalNodes=\"true\" rootNode=\"false\">\n" +
                        "		<parameter id=\"treeModel1.internalNodeHeights\"/>\n" +
                        "	</nodeHeights>\n" +
                        "	<nodeHeights internalNodes=\"true\" rootNode=\"true\">\n" +
                        "		<parameter id=\"treeModel1.allInternalNodeHeights\"/>\n" +
                        "	</nodeHeights>\n" +
                        "</treeModel>";

    }

    public Class getReturnType() {
        return StarTreeModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules;
}
