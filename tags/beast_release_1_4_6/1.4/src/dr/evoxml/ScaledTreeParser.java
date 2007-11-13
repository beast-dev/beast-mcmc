/*
 * ScaledTreeParser.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evoxml;

import dr.evolution.distance.DistanceMatrix;
import dr.evolution.tree.*;
import dr.evolution.util.TimeScale;
import dr.evolution.util.TaxonList;
import dr.xml.*;
import dr.evomodel.sitemodel.SiteModel;

import java.util.Set;
import java.util.Map;
import java.util.HashMap;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class ScaledTreeParser extends AbstractXMLObjectParser {

	//
	// Public stuff
	//

	public static final String SCALED_TREE = "scaledTree";
    public static final String ROOT_UPPER = "rootUpper";
    public static final String ROOT_LOWER = "rootLower";
    public static final String MRCA = "mrca";
    public static final String UPPER = "upper";
    public static final String LOWER = "lower";

	public String getParserName() { return SCALED_TREE; }

	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

	    FlexibleTree tree = new FlexibleTree((Tree)xo.getChild(Tree.class));

        // set the tips to their assigned dates...
        dr.evolution.util.Date mostRecent = null;
        for (int i = 0; i < tree.getTaxonCount(); i++) {

            dr.evolution.util.Date date = (dr.evolution.util.Date)tree.getTaxonAttribute(i, dr.evolution.util.Date.DATE);

            if (date == null) {
                date = (dr.evolution.util.Date)tree.getNodeAttribute(tree.getExternalNode(i), dr.evolution.util.Date.DATE);
            }

            if (date != null && ((mostRecent == null) || date.after(mostRecent))) {
                mostRecent = date;
            }
        }

        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            dr.evolution.util.Date date = (dr.evolution.util.Date)tree.getNodeAttribute(tree.getInternalNode(i), dr.evolution.util.Date.DATE);

            if (date != null && ((mostRecent == null) || date.after(mostRecent))) {
                mostRecent = date;
            }
        }

        if (mostRecent != null) {
            TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());

            for (int i = 0; i < tree.getTaxonCount(); i++) {
                dr.evolution.util.Date date = (dr.evolution.util.Date)tree.getTaxonAttribute(i, dr.evolution.util.Date.DATE);

                if (date == null) {
                    date = (dr.evolution.util.Date)tree.getNodeAttribute(tree.getExternalNode(i), dr.evolution.util.Date.DATE);
                }

                if (date != null) {
                    double height = timeScale.convertTime(date.getTimeValue(), date);
                    tree.setNodeHeight(tree.getExternalNode(i), height);
                }
            }

            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                dr.evolution.util.Date date = (dr.evolution.util.Date)tree.getNodeAttribute(tree.getInternalNode(i), dr.evolution.util.Date.DATE);

                if (date != null) {
                    double height = timeScale.convertTime(date.getTimeValue(), date);
                    tree.setNodeHeight(tree.getInternalNode(i), height);
                }
            }

            MutableTree.Utils.correctHeightsForTips(tree);
        }

        Map nodeConstraints = new HashMap();

        // First add the constraints on the root...
        NodeConstraint constraint = new NodeConstraint();

		if (xo.hasAttribute(ROOT_UPPER)) {
			constraint.upper = xo.getDoubleAttribute(ROOT_UPPER);
		}

        if (xo.hasAttribute(ROOT_LOWER)) {
            constraint.lower = xo.getDoubleAttribute(ROOT_LOWER);
        }

        if (constraint.upper < constraint.lower) {
            throw new XMLParseException("Error parsing " + getParserName() + ": constraints on the root are not compatible.");
        }

        nodeConstraints.put(tree.getRoot(), constraint);

        // Now add any node constraints...
        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof XMLObject) {

                XMLObject xoc = (XMLObject)xo.getChild(i);
                if (xoc.getName().equals(MRCA)) {

                    TaxonList taxonList = (TaxonList)xoc.getChild(TaxonList.class);

                    double upper = Double.POSITIVE_INFINITY;
                    double lower = Double.NEGATIVE_INFINITY;

                    if (xoc.hasAttribute(UPPER)) {
                        upper = xoc.getDoubleAttribute(UPPER);
                    }

                    if (xoc.hasAttribute(LOWER)) {
                        lower = xoc.getDoubleAttribute(LOWER);
                    }

                    try {
                        Set leaves = Tree.Utils.getLeavesForTaxa(tree, taxonList);
                        NodeRef node = Tree.Utils.getCommonAncestorNode(tree, leaves);

                        constraint = (NodeConstraint)nodeConstraints.get(node);
                        if (constraint == null) {
                            constraint = new NodeConstraint();
                        }

                        if (upper < constraint.upper) {
                            constraint.upper = upper;
                        }

                        if (lower > constraint.lower) {
                            constraint.lower = lower;
                        }

                        if (constraint.upper < constraint.lower) {
                            throw new XMLParseException("Error parsing " + getParserName() + ": constraints on an MRCA are not compatible.");
                        }
                        nodeConstraints.put(node, constraint);

                    } catch (Tree.MissingTaxonException e) {
                        throw new XMLParseException("Error parsing " + getParserName() + ": " + e.getMessage());
                    }
                }

            }
        }

        rescaleNode(tree, tree.getRoot(), nodeConstraints);

        return tree;
	}

    private void rescaleNode(Tree tree, NodeRef node, Map nodeConstraints) {
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);
            if (!tree.isExternal(child)) {
                rescaleNode(tree, child, nodeConstraints);
            }
        }

        NodeConstraint constraint = (NodeConstraint)nodeConstraints.get(node);
        if (constraint != null) {

        }
    }



    class NodeConstraint {
        double upper = Double.POSITIVE_INFINITY;
        double lower = Double.NEGATIVE_INFINITY;
    }

    //************************************************************************
	// AbstractXMLObjectParser implementation
	//************************************************************************

	public XMLSyntaxRule[] getSyntaxRules() { return rules; }

	private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			AttributeRule.newDoubleRule(ROOT_UPPER, true),
            AttributeRule.newDoubleRule(ROOT_LOWER, true),
			new ElementRule(Tree.class, "The source tree"),
            new ElementRule(MRCA, new XMLSyntaxRule[] {
                    AttributeRule.newDoubleRule(UPPER, true),
                    AttributeRule.newDoubleRule(LOWER, true),
                    new ElementRule(TaxonList.class, "A list of taxa that defines the MRCA")
            }, "A constraint on an MRCA", 0, Integer.MAX_VALUE)
    };

	public String getParserDescription() {
		return "This element returns a tree with the same topology as the source but with node heights scaled to fit a set of constraints.";
	}

	public Class getReturnType() { return UPGMATree.class; }
}