/*
 * UPGMATreeParser.java
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

import dr.evolution.distance.DistanceMatrix;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.UPGMATree;
import dr.evolution.util.TimeScale;
import dr.evomodelxml.tree.TreeModelParser;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 * @version $Id: UPGMATreeParser.java,v 1.6 2006/07/28 11:27:32 rambaut Exp $
 */
public class UPGMATreeParser extends AbstractXMLObjectParser {

    //
    // Public stuff
    //

    public static final String UPGMA_TREE = "upgmaTree";
    //public static final String DISTANCES = "distances";
    public static final String ROOT_HEIGHT = TreeModelParser.ROOT_HEIGHT;
    public static final String RANDOMIZE = "nonzeroBranchLengths";

    public String getParserName() {
        return UPGMA_TREE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean usingDatesSpecified = false;
        boolean usingDates = true;
        double rootHeight = xo.getAttribute(ROOT_HEIGHT, -1.0);

        if (xo.hasAttribute(SimpleTreeParser.USING_DATES)) {
            usingDatesSpecified = true;
            usingDates = xo.getBooleanAttribute(SimpleTreeParser.USING_DATES);
        }

        DistanceMatrix distances = (DistanceMatrix) xo.getChild(DistanceMatrix.class);

        UPGMATree tree = new UPGMATree(distances);

        if (rootHeight > 0) {
            double scaleFactor = rootHeight / tree.getNodeHeight(tree.getRoot());

            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                NodeRef node = tree.getInternalNode(i);
                double height = tree.getNodeHeight(node);
                tree.setNodeHeight(node, height * scaleFactor);
            }
        }

        if (usingDates) {

            dr.evolution.util.Date mostRecent = null;
            for (int i = 0; i < tree.getTaxonCount(); i++) {

                dr.evolution.util.Date date = (dr.evolution.util.Date) tree.getTaxonAttribute(i, dr.evolution.util.Date.DATE);

                if (date == null) {
                    date = (dr.evolution.util.Date) tree.getNodeAttribute(tree.getExternalNode(i), dr.evolution.util.Date.DATE);
                }

                if (date != null && ((mostRecent == null) || date.after(mostRecent))) {
                    mostRecent = date;
                }
            }

            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                dr.evolution.util.Date date = (dr.evolution.util.Date) tree.getNodeAttribute(tree.getInternalNode(i), dr.evolution.util.Date.DATE);

                if (date != null && ((mostRecent == null) || date.after(mostRecent))) {
                    mostRecent = date;
                }
            }

            if (mostRecent == null) {
                if (usingDatesSpecified) {
                    throw new XMLParseException("no date elements in tree (and usingDates attribute set)");
                }
            } else {
                TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());

                for (int i = 0; i < tree.getTaxonCount(); i++) {
                    dr.evolution.util.Date date = (dr.evolution.util.Date) tree.getTaxonAttribute(i, dr.evolution.util.Date.DATE);

                    if (date == null) {
                        date = (dr.evolution.util.Date) tree.getNodeAttribute(tree.getExternalNode(i), dr.evolution.util.Date.DATE);
                    }

                    if (date != null) {
                        double height = timeScale.convertTime(date.getTimeValue(), date);
                        tree.setNodeHeight(tree.getExternalNode(i), height);
                    }
                }

                for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                    dr.evolution.util.Date date = (dr.evolution.util.Date) tree.getNodeAttribute(tree.getInternalNode(i), dr.evolution.util.Date.DATE);

                    if (date != null) {
                        double height = timeScale.convertTime(date.getTimeValue(), date);
                        tree.setNodeHeight(tree.getInternalNode(i), height);
                    }
                }

                MutableTree.Utils.correctHeightsForTips(tree);
            }

        }

        if (rootHeight > 0) {
            double scaleFactor = rootHeight / tree.getNodeHeight(tree.getRoot());

            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                NodeRef node = tree.getInternalNode(i);
                double height = tree.getNodeHeight(node);
                tree.setNodeHeight(node, height * scaleFactor);
            }
        }

        if (xo.getAttribute(RANDOMIZE, false)) {
            shakeTree(tree);
        }


        return tree;
    }

    private boolean shakeNode(UPGMATree tree, NodeRef node) {
        if (tree.isRoot(node) || tree.isExternal(node)) {
            return false;
        }

        boolean shake = false;
        if (tree.getBranchLength(node) <= tolerance) {
            shake = true;
        }
        double maxHeight = tree.getNodeHeight(tree.getParent(node));
        double minHeight = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);
            if (tree.getBranchLength(child) <= tolerance) {
                shake = true;
            }
            double thisHeight = tree.getNodeHeight(child);
            if (thisHeight > minHeight) {
                minHeight = thisHeight;
            }
        }
        if (shake) {
            double draw = minHeight + (maxHeight - minHeight) * MathUtils.nextDouble();
            tree.setNodeHeight(node, draw);
        }
        return shake;
    }

    private void shakeTree(UPGMATree tree) {
        boolean shake = true;
        int[] permutation = new int[tree.getNodeCount()];
        for (int i = 0; i < tree.getNodeCount(); i++) {
            permutation[i] = i;
        }
        while (shake) {
            Logger.getLogger("dr.evomodelxml").info("Adjusting heights in UPGMA tree");
            MathUtils.permute(permutation);
            shake = false;
            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(permutation[i]);
                if (shakeNode(tree, node)) {
                    shake = true;
                }
            }
        }
    }

    private static double tolerance = 0.0;

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(SimpleTreeParser.USING_DATES, true),
            AttributeRule.newDoubleRule(ROOT_HEIGHT, true),
            AttributeRule.newBooleanRule(RANDOMIZE, true),
            new ElementRule(DistanceMatrix.class)
    };

    public String getParserDescription() {
        return "This element returns a UPGMA tree generated from the given distances.";
    }

    public Class getReturnType() {
        return UPGMATree.class;
    }
}
