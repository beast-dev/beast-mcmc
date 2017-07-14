/*
 * NewickParser.java
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

package dr.evoxml;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evoxml.util.XMLUnits;
import dr.xml.*;

import java.io.IOException;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: NewickParser.java,v 1.7 2006/04/25 14:41:08 rambaut Exp $
 */
public class NewickParser extends AbstractXMLObjectParser {

    public static final String NEWICK = "newick";
    public static final String UNITS = "units";
    public static final String RESCALE_HEIGHT = "rescaleHeight";
    public static final String RESCALE_LENGTH = "rescaleLength";
    public static final String USING_DATES = SimpleTreeParser.USING_DATES;
    public static final String USING_HEIGHTS = "usingHeights";

    public NewickParser() {
        rules = new XMLSyntaxRule[]{
                AttributeRule.newBooleanRule(USING_DATES, true),
                AttributeRule.newBooleanRule(USING_HEIGHTS, true),
                AttributeRule.newDoubleRule(RESCALE_HEIGHT, true, "Attempt to rescale the tree to the given root height"),
                AttributeRule.newDoubleRule(RESCALE_LENGTH, true, "Attempt to rescale the tree to the given total length"),
                new StringAttributeRule(UNITS, "The branch length units of this tree", Units.UNIT_NAMES, true),
                new ElementRule(String.class, "The NEWICK format tree. Tip labels are taken to be Taxon IDs")
        };
    }

    public String getParserName() {
        return NEWICK;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

//        boolean usingDates = xo.getAttribute(USING_DATES, true);

        boolean usingDates = true;
        if (xo.hasAttribute(USING_DATES)) {
            usingDates = xo.getAttribute(USING_DATES, true);
        }

        boolean usingHeights = false;
        if (xo.hasAttribute(USING_HEIGHTS)) {
            usingHeights = xo.getAttribute(USING_HEIGHTS, true);
        }

//        System.out.println("UsingDates=" + usingDates + " usingHeights= " + usingHeights);

        if (usingDates && usingHeights) {
            throw new XMLParseException("Unable to use both dates and node heights. Specify value of usingDates attribute.");
        }
//		else if (!usingDates && !usingHeights) {
//			System.out.println("Tree is assumed to be ultrametric");
//		}

        StringBuffer buffer = new StringBuffer();

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof String) {
                buffer.append((String) xo.getChild(i));
            } else {
                throw new XMLParseException("illegal element in newick element");
            }
        }

        java.io.Reader reader = new java.io.StringReader(buffer.toString());
        NewickImporter importer = new NewickImporter(reader);

        FlexibleTree tree;

        try {
            tree = (FlexibleTree) importer.importTree(null);
        } catch (IOException ioe) {
            throw new XMLParseException("error parsing tree in newick element");
        } catch (NewickImporter.BranchMissingException bme) {
            throw new XMLParseException("branch missing in tree in newick element");
        } catch (Importer.ImportException ime) {
            throw new XMLParseException("error parsing tree in newick element - " + ime.getMessage());
        }

        if (tree == null) {
            throw new XMLParseException("Failed to read tree");
        }

        tree.setUnits(units);

        for (int i = 0; i < tree.getTaxonCount(); i++) {

            FlexibleNode node = (FlexibleNode) tree.getExternalNode(i);

            String id = node.getTaxon().getId();
            Taxon taxon = null;

            XMLObject obj = getStore().get(id);

            if (obj != null && obj.getNativeObject() instanceof Taxon) {

                taxon = (Taxon) obj.getNativeObject();
            }

            if (taxon != null) {

                node.setTaxon(taxon);

            } else {
                throw new XMLParseException("unknown taxon, " + id + ", in newick tree");
            }
        }

        if (usingDates) {

            // are all the tips just being translated by a fixed amount?
            // in which case we can just translate the internal nodes.
            double fixedDiff = 0.0;
            boolean translateNodes = true;

            for (int i = 0; i < tree.getTaxonCount(); i++) {

                NodeRef node = tree.getExternalNode(i);

                dr.evolution.util.Date date = (dr.evolution.util.Date) tree.getTaxonAttribute(i, dr.evolution.util.Date.DATE);

                if (date == null) {
                    date = (dr.evolution.util.Date) tree.getNodeAttribute(tree.getExternalNode(i), dr.evolution.util.Date.DATE);
                }

                double height = 0.0;
                double nodeHeight = tree.getNodeHeight(node);
                if (date != null) {
                    height = Taxon.getHeightFromDate(date);
                }

                double diff = height - nodeHeight;

                if (i == 0) {
                    fixedDiff = diff;
                } else if (Math.abs(diff - fixedDiff) > 1e-5) {
                    translateNodes = false;
                }

                if (Math.abs(diff) > 1e-8 && (i == 0 || !translateNodes) ) {

                    System.out.println("  Changing height of node " + tree.getTaxon(node.getNumber()) + " from " + nodeHeight + " to " + height);
                    tree.setNodeHeight(node, height);
                }
            }

            if (translateNodes) {
                System.out.println("  Changing height of all nodes by " + fixedDiff);
            }

            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                NodeRef node = tree.getInternalNode(i);

                dr.evolution.util.Date date = (dr.evolution.util.Date) tree.getNodeAttribute(node, dr.evolution.util.Date.DATE);

                if (date != null) {
                    double height = Taxon.getHeightFromDate(date);
                    tree.setNodeHeight(node, height);
                } else if (translateNodes) {
                    tree.setNodeHeight(node, tree.getNodeHeight(node) + fixedDiff);
                }

            }// END: i loop

            MutableTree.Utils.correctHeightsForTips(tree);

        } else if (!usingDates && !usingHeights) {

            System.out.println("Tree is assumed to be ultrametric");

            // not using dates or heights
            for (int i = 0; i < tree.getTaxonCount(); i++) {
                final NodeRef leaf = tree.getExternalNode(i);
                final double h = tree.getNodeHeight(leaf);

                if (h != 0.0) {
                    double zero = 0.0;
                    System.out.println("  Changing height of leaf node " + tree.getTaxon(leaf.getNumber()) + " from " + h + " to " + zero);
                    tree.setNodeHeight(leaf, zero);
                }

            }// END: i loop

        } else {

            System.out.println("Using node heights.");

        }// END: usingDates check

        if (xo.hasAttribute(RESCALE_HEIGHT)) {
            double rescaleHeight = xo.getDoubleAttribute(RESCALE_HEIGHT);
            double scale = rescaleHeight / tree.getNodeHeight(tree.getRoot());
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                NodeRef n = tree.getInternalNode(i);
                tree.setNodeHeight(n, tree.getNodeHeight(n) * scale);
            }
        }

        if (xo.hasAttribute(RESCALE_LENGTH)) {
            double rescaleLength = xo.getDoubleAttribute(RESCALE_LENGTH);
            double scale = rescaleLength / TreeUtils.getTreeLength(tree, tree.getRoot());
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                NodeRef n = tree.getInternalNode(i);
                tree.setNodeHeight(n, tree.getNodeHeight(n) * scale);
            }
        }
        //System.out.println("Constructed newick tree = " + Tree.Utils.uniqueNewick(tree, tree.getRoot()));
        return tree;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Constructs a tree from a NEWICK format tree description";
    }

    public String getExample() {
        return "<" + getParserName() + " " + UNITS + "=\"" + Units.Utils.getDefaultUnitName(Units.Type.YEARS) + "\">" + " ((A:1.0, B:1.0):1.0,(C:2.0, D:2.0):1.0); </" + getParserName() + ">";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules;

    public Class getReturnType() {
        return Tree.class;
    }
}
