/*
 * TreeTrace.java
 *
 * Copyright (C) 2002-2007 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 */
public class TreeModelParser extends AbstractXMLObjectParser {

    public static final String ROOT_HEIGHT = "rootHeight";
    public static final String LEAF_HEIGHT = "leafHeight";
    public static final String LEAF_TRAIT = "leafTrait";

    public static final String NODE_HEIGHTS = "nodeHeights";
    public static final String NODE_RATES = "nodeRates";
    public static final String NODE_TRAITS = "nodeTraits";
    public static final String MULTIVARIATE_TRAIT = "traitDimension";

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
        TreeModel treeModel = new TreeModel(tree);

        Logger.getLogger("dr.evomodel").info("Creating the tree model, '" + xo.getId() + "'");

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof XMLObject) {

                XMLObject cxo = (XMLObject) xo.getChild(i);

                if (cxo.getName().equals(ROOT_HEIGHT)) {

                    replaceParameter(cxo, treeModel.getRootHeightParameter());

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
                    replaceParameter(cxo, treeModel.getLeafHeightParameter(node));

                } else if (cxo.getName().equals(NODE_HEIGHTS)) {

                    boolean rootNode = false;
                    boolean internalNodes = false;
                    boolean leafNodes = false;

                    if (cxo.hasAttribute(ROOT_NODE)) {
                        rootNode = cxo.getBooleanAttribute(ROOT_NODE);
                    }

                    if (cxo.hasAttribute(INTERNAL_NODES)) {
                        internalNodes = cxo.getBooleanAttribute(INTERNAL_NODES);
                    }

                    if (cxo.hasAttribute(LEAF_NODES)) {
                        leafNodes = cxo.getBooleanAttribute(LEAF_NODES);
                    }

                    if (!rootNode && !internalNodes && !leafNodes) {
                        throw new XMLParseException("one or more of root, internal or leaf nodes must be selected for the nodeHeights element");
                    }

                    replaceParameter(cxo, treeModel.createNodeHeightsParameter(rootNode, internalNodes, leafNodes));

                } else if (cxo.getName().equals(NODE_RATES)) {

                    boolean rootNode = false;
                    boolean internalNodes = false;
                    boolean leafNodes = false;

                    if (cxo.hasAttribute(ROOT_NODE)) {
                        rootNode = cxo.getBooleanAttribute(ROOT_NODE);
                    }

                    if (cxo.hasAttribute(INTERNAL_NODES)) {
                        internalNodes = cxo.getBooleanAttribute(INTERNAL_NODES);
                    }

                    if (cxo.hasAttribute(LEAF_NODES)) {
                        leafNodes = cxo.getBooleanAttribute(LEAF_NODES);
                    }

                    //if (rootNode) {
                    //	throw new XMLParseException("root node does not have a rate parameter");
                    //}

                    if (!rootNode && !internalNodes && !leafNodes) {
                        throw new XMLParseException("one or more of root, internal or leaf nodes must be selected for the nodeRates element");
                    }

                    replaceParameter(cxo, treeModel.createNodeRatesParameter(rootNode, internalNodes, leafNodes));

                } else if (cxo.getName().equals(NODE_TRAITS)) {

                    boolean rootNode = false;
                    boolean internalNodes = false;
                    boolean leafNodes = false;
                    String name = "trait";

                    int dim = 1;

                    if (cxo.hasAttribute(MULTIVARIATE_TRAIT)) {
                        dim = cxo.getIntegerAttribute(MULTIVARIATE_TRAIT);
                    }

                    if (cxo.hasAttribute(ROOT_NODE)) {
                        rootNode = cxo.getBooleanAttribute(ROOT_NODE);
                    }

                    if (cxo.hasAttribute(INTERNAL_NODES)) {
                        internalNodes = cxo.getBooleanAttribute(INTERNAL_NODES);
                    }

                    if (cxo.hasAttribute(LEAF_NODES)) {
                        leafNodes = cxo.getBooleanAttribute(LEAF_NODES);
                    }

                    if (cxo.hasAttribute(NAME)) {
                        name = cxo.getStringAttribute(NAME);
                    }

                    if (!rootNode && !internalNodes && !leafNodes) {
                        throw new XMLParseException("one or more of root, internal or leaf nodes must be selected for the nodeTraits element");
                    }

                    replaceParameter(cxo, treeModel.createNodeTraitsParameter(name, dim, rootNode, internalNodes, leafNodes));

                } else if (cxo.getName().equals(LEAF_TRAIT)) {

                    String name = "trait";

                    String taxonName;
                    if (cxo.hasAttribute(TAXON)) {
                        taxonName = cxo.getStringAttribute(TAXON);
                    } else {
                        throw new XMLParseException("taxa element missing from leafTrait element in treeModel element");
                    }

                    int index = treeModel.getTaxonIndex(taxonName);
                    if (index == -1) {
                        throw new XMLParseException("taxon " + taxonName + " not found for leafTrait element in treeModel element");
                    }
                    NodeRef node = treeModel.getExternalNode(index);

                    Parameter parameter = treeModel.getNodeTraitParameter(node, name);
                    replaceParameter(cxo, parameter);


                } else {
                    throw new XMLParseException("illegal child element in " + getParserName() + ": " + cxo.getName());
                }

            } else if (xo.getChild(i) instanceof Tree) {
                // do nothing - already handled
            } else {
                throw new XMLParseException("illegal child element in  " + getParserName() + ": " + xo.getChildName(i) + " " + xo.getChild(i));
            }
        }

        treeModel.setupHeightBounds();

        Logger.getLogger("dr.evomodel").info("  initial tree topology = " + Tree.Utils.uniqueNewick(treeModel, treeModel.getRoot()));
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
        return TreeModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Tree.class),
            new ElementRule(ROOT_HEIGHT, Parameter.class, "A parameter definition with id only (cannot be a reference!)", false),
            new ElementRule(NODE_HEIGHTS,
                    new XMLSyntaxRule[]{
                            AttributeRule.newBooleanRule(ROOT_NODE, true, "If true the root height is included in the parameter"),
                            AttributeRule.newBooleanRule(INTERNAL_NODES, true, "If true the internal node heights (minus the root) are included in the parameter"),
                            new ElementRule(Parameter.class, "A parameter definition with id only (cannot be a reference!)")
                    }, 1, Integer.MAX_VALUE)
    };

    public Parameter getParameter(XMLObject xo) throws XMLParseException {

        int paramCount = 0;
        Parameter param = null;
        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof Parameter) {
                param = (Parameter) xo.getChild(i);
                paramCount += 1;
            }
        }

        if (paramCount == 0) {
            throw new XMLParseException("no parameter element in treeModel " + xo.getName() + " element");
        } else if (paramCount > 1) {
            throw new XMLParseException("More than one parameter element in treeModel " + xo.getName() + " element");
        }

        return param;
    }

    public void replaceParameter(XMLObject xo, Parameter newParam) throws XMLParseException {

        for (int i = 0; i < xo.getChildCount(); i++) {

            if (xo.getChild(i) instanceof Parameter) {

                XMLObject rxo;
                Object obj = xo.getRawChild(i);

                if (obj instanceof Reference) {
                    rxo = ((Reference) obj).getReferenceObject();
                } else if (obj instanceof XMLObject) {
                    rxo = (XMLObject) obj;
                } else {
                    throw new XMLParseException("object reference not available");
                }

                if (rxo.getChildCount() > 0) {
                    throw new XMLParseException("No child elements allowed in parameter element.");
                }

                if (rxo.hasAttribute(XMLParser.IDREF)) {
                    throw new XMLParseException("References to " + xo.getName() + " parameters are not allowed in treeModel.");
                }

                if (rxo.hasAttribute(ParameterParser.VALUE)) {
                    throw new XMLParseException("Parameters in " + xo.getName() + " have values set automatically.");
                }

                if (rxo.hasAttribute(ParameterParser.UPPER)) {
                    throw new XMLParseException("Parameters in " + xo.getName() + " have bounds set automatically.");
                }

                if (rxo.hasAttribute(ParameterParser.LOWER)) {
                    throw new XMLParseException("Parameters in " + xo.getName() + " have bounds set automatically.");
                }

                if (rxo.hasAttribute(XMLParser.ID)) {

                    newParam.setId(rxo.getStringAttribute(XMLParser.ID));
                }

                rxo.setNativeObject(newParam);

                return;
            }
        }
    }
}
