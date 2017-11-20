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
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
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
    public static final String INITIAL_VALUE = "initialValue";
    public static final String AS_MATRIX = "asMatrix";

    public static final String ROOT_NODE = "rootNode";
    public static final String INTERNAL_NODES = "internalNodes";
    public static final String LEAF_NODES = "leafNodes";

    public static final String LEAF_HEIGHTS = "leafHeights";

    public static final String FIRE_TREE_EVENTS = "fireTreeEvents";
    public static final String FIX_HEIGHTS = "fixHeights";
    public static final String FIX_TREE = "fixTree";

    public static final String TAXON = "taxon";
    public static final String NAME = "name";

    public static ElementRule nodeTraitsRule = new ElementRule(NODE_TRAITS,
            new XMLSyntaxRule[]{
                    AttributeRule.newStringRule(NAME, false, "The name of the trait attribute in the taxa"),
                    AttributeRule.newBooleanRule(ROOT_NODE, true, "If true the root trait is included in the parameter"),
                    AttributeRule.newBooleanRule(INTERNAL_NODES, true, "If true the internal node traits (minus the root) are included in the parameter"),
                    AttributeRule.newBooleanRule(LEAF_NODES, true, "If true the leaf node traits are included in the parameter"),
                    AttributeRule.newIntegerRule(MULTIVARIATE_TRAIT, true, "The number of dimensions (if multivariate)"),
                    AttributeRule.newDoubleRule(INITIAL_VALUE, true, "The initial value(s)"),
                    AttributeRule.newBooleanRule(FIRE_TREE_EVENTS, true, "Whether to fire tree events if the traits change"),
                    AttributeRule.newBooleanRule(AS_MATRIX, true, "Whether to return parameter as a matrix"),
                    new ElementRule(Parameter.class, "A parameter definition with id only (cannot be a reference!)")
            }, 0, Integer.MAX_VALUE);

    public TreeModelParser() {
        rules = new XMLSyntaxRule[]{
                new ElementRule(Tree.class),
                new ElementRule(ROOT_HEIGHT, Parameter.class, "A parameter definition with id only (cannot be a reference!)", false),
                AttributeRule.newBooleanRule(FIX_HEIGHTS, true),
                AttributeRule.newBooleanRule(FIX_TREE, true),
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
                        }, 0, Integer.MAX_VALUE),
                nodeTraitsRule,
                new ElementRule(NODE_RATES,
                        new XMLSyntaxRule[]{
                                AttributeRule.newBooleanRule(ROOT_NODE, true, "If true the root rate is included in the parameter"),
                                AttributeRule.newBooleanRule(INTERNAL_NODES, true, "If true the internal node rate (minus the root) are included in the parameter"),
                                AttributeRule.newBooleanRule(LEAF_NODES, true, "If true the leaf node rate are included in the parameter"),
                                AttributeRule.newDoubleRule(INITIAL_VALUE, true, "The initial value(s)"),
                                new ElementRule(Parameter.class, "A parameter definition with id only (cannot be a reference!)")
                        }, 0, Integer.MAX_VALUE),
                new ElementRule(LEAF_TRAIT,
                        new XMLSyntaxRule[]{
                                AttributeRule.newStringRule(TAXON, false, "The name of the taxon for the leaf"),
                                AttributeRule.newStringRule(NAME, false, "The name of the trait attribute in the taxa"),
                                new ElementRule(Parameter.class, "A parameter definition with id only (cannot be a reference!)")
                        }, 0, Integer.MAX_VALUE),
                new ElementRule(LEAF_HEIGHTS,
                        new XMLSyntaxRule[]{
                                new ElementRule(TaxonList.class, "A set of taxa for which leaf heights are required"),
                                new ElementRule(Parameter.class, "A compound parameter containing the leaf heights")
                        }, true)
        };
    }

    public String getParserName() {
        return TreeModel.TREE_MODEL;
    }

    /**
     * @return a tree object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree) xo.getChild(Tree.class);
        boolean fixHeights = xo.getAttribute(FIX_HEIGHTS, false);
        boolean fixTree = xo.getAttribute(FIX_TREE, false);

        TreeModel treeModel = new TreeModel(xo.getId(), tree, fixHeights, fixTree);

        Logger.getLogger("dr.evomodel").info("\nCreating the tree model, '" + xo.getId() + "'");

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof XMLObject) {

                XMLObject cxo = (XMLObject) xo.getChild(i);

                if (cxo.getName().equals(ROOT_HEIGHT)) {

                    ParameterParser.replaceParameter(cxo, treeModel.getRootHeightParameter());

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

                    Parameter newParameter = treeModel.getLeafHeightParameter(node);

                    ParameterParser.replaceParameter(cxo, newParameter);

                    Taxon taxon = treeModel.getTaxon(index);

                    setUncertaintyBounds(newParameter, taxon);

                } else if (cxo.getName().equals(LEAF_HEIGHTS)) {
                    // get a set of leaf height parameters out as a compound parameter...

                    TaxonList taxa = (TaxonList)cxo.getChild(TaxonList.class);
                    Parameter offsetParameter = (Parameter)cxo.getChild(Parameter.class);

                    CompoundParameter leafHeights = new CompoundParameter("leafHeights");
                    for (Taxon taxon : taxa) {
                        int index = treeModel.getTaxonIndex(taxon);
                        if (index == -1) {
                            throw new XMLParseException("taxon " + taxon.getId() + " not found for leafHeight element in treeModel element");
                        }
                        NodeRef node = treeModel.getExternalNode(index);

                        Parameter newParameter = treeModel.getLeafHeightParameter(node);

                        leafHeights.addParameter(newParameter);

                        setUncertaintyBounds(newParameter, taxon);
                    }

                    ParameterParser.replaceParameter(cxo, leafHeights);

                } else if (cxo.getName().equals(NODE_HEIGHTS)) {

                    boolean rootNode = cxo.getAttribute(ROOT_NODE, false);
                    boolean internalNodes = cxo.getAttribute(INTERNAL_NODES, false);
                    boolean leafNodes = cxo.getAttribute(LEAF_NODES, false);

                    if (!rootNode && !internalNodes && !leafNodes) {
                        throw new XMLParseException("one or more of root, internal or leaf nodes must be selected for the nodeHeights element");
                    }

                    ParameterParser.replaceParameter(cxo, treeModel.createNodeHeightsParameter(rootNode, internalNodes, leafNodes));

                } else if (cxo.getName().equals(NODE_RATES)) {

                    boolean rootNode = cxo.getAttribute(ROOT_NODE, false);
                    boolean internalNodes = cxo.getAttribute(INTERNAL_NODES, false);
                    boolean leafNodes = cxo.getAttribute(LEAF_NODES, false);
                    double[] initialValues = null;

                    if (cxo.hasAttribute(INITIAL_VALUE)) {
                        initialValues = cxo.getDoubleArrayAttribute(INITIAL_VALUE);
                    }

                    if (!rootNode && !internalNodes && !leafNodes) {
                        throw new XMLParseException("one or more of root, internal or leaf nodes must be selected for the nodeRates element");
                    }

                    ParameterParser.replaceParameter(cxo, treeModel.createNodeRatesParameter(initialValues, rootNode, internalNodes, leafNodes));

                } else if (cxo.getName().equals(NODE_TRAITS)) {

                    parseNodeTraits(cxo, treeModel);

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
        
//        Logger.getLogger("dr.evomodel").info("  initial tree topology = " + TreeUtils.uniqueNewick(treeModel, treeModel.getRoot()));
        Logger.getLogger("dr.evomodel").info("  taxon count = " + treeModel.getExternalNodeCount());
        Logger.getLogger("dr.evomodel").info("  tree height = " + treeModel.getNodeHeight(treeModel.getRoot()));
        return treeModel;
    }

    public static void parseNodeTraits(XMLObject cxo, TreeModel treeModel) throws XMLParseException {

        boolean rootNode = cxo.getAttribute(ROOT_NODE, false);
        boolean internalNodes = cxo.getAttribute(INTERNAL_NODES, false);
        boolean leafNodes = cxo.getAttribute(LEAF_NODES, false);
        boolean fireTreeEvents = cxo.getAttribute(FIRE_TREE_EVENTS, false);
        boolean asMatrix = cxo.getAttribute(AS_MATRIX, false);
        String name = cxo.getAttribute(NAME, "trait");
        int dim = cxo.getAttribute(MULTIVARIATE_TRAIT, 1);

        double[] initialValues = null;
        if (cxo.hasAttribute(INITIAL_VALUE)) {
            initialValues = cxo.getDoubleArrayAttribute(INITIAL_VALUE);
        }

        if (!rootNode && !internalNodes && !leafNodes) {
            throw new XMLParseException("one or more of root, internal or leaf nodes must be selected for the nodeTraits element");
        }

        Parameter newParameter = asMatrix ?
                treeModel.createNodeTraitsParameterAsMatrix(name, dim, initialValues, rootNode, internalNodes, leafNodes, fireTreeEvents) :
                treeModel.createNodeTraitsParameter(name, dim, initialValues, rootNode, internalNodes, leafNodes, fireTreeEvents);

        ParameterParser.replaceParameter(cxo, newParameter);
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

    private final XMLSyntaxRule[] rules;
}
