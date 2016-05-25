/*
 * TreeUtils.java
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

package dr.app.tempest;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TimeScale;
import dr.stats.Variate;

/**
 * Class for getting information from trees.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeUtils.java,v 1.7 2005/07/11 14:07:25 rambaut Exp $
 */
public class TreeUtils {

    /**
     * Gets the root to tip distances from a tree.
     */
    public static void getRootToTipDistances(Tree tree, Variate distances) {

        double rootHeight = tree.getNodeHeight(tree.getRoot());
        double height;

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {

            height = tree.getNodeHeight(tree.getExternalNode(i));
            distances.add(rootHeight - height);
        }
    }

    /**
     * Gets the tip dates from a tree.
     */
    public static void getTipDates(Tree tree, Variate dates) {

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {

            Taxon taxon = tree.getNodeTaxon(tree.getExternalNode(i));
            Object date = taxon.getAttribute("date");
            if (date != null) {
                if (date instanceof Date) {
                    dates.add(((Date) date).getTimeValue());
                } else {
                    try {
                        dates.add(Double.parseDouble(date.toString()));
                    } catch (NumberFormatException nfe) {
                        dates.add(0.0);
                    }
                }
            } else {
                dates.add(0.0);
            }
        }
    }

    /**
     * Gets the root to tip distances from a node.
     */
    public static void getRootToTipDistances(Tree tree, NodeRef node, double rootHeight,
                                             Variate distances) {

        if (tree.isExternal(node)) {
            double height = tree.getNodeHeight(node);
            distances.add(rootHeight - height);
        } else {

            getRootToTipDistances(tree, tree.getChild(node, 0), rootHeight, distances);
            getRootToTipDistances(tree, tree.getChild(node, 1), rootHeight, distances);
        }
    }

    /**
     * Gets the tip dates from a node.
     */
    public static void getTipDates(Tree tree, NodeRef node, Variate dates) {

        if (tree.isExternal(node)) {

            String date = (String) tree.getNodeTaxon(node).getAttribute("date");
            if (date != null) {
                dates.add(Double.parseDouble(date));
            } else {
                dates.add(0.0);
            }
        } else {

            getTipDates(tree, tree.getChild(node, 0), dates);
            getTipDates(tree, tree.getChild(node, 1), dates);
        }
    }

    /**
     * Returns -1 if there is no number suffix
     */
    public static double guessDate(String s) {

        int i = s.length();
        char c;
        do {
            i--;
            c = s.charAt(i);
        } while (i >= 0 && (Character.isDigit(c) || c == '.'));

        if (i == s.length()) {
            return 0.0;
        }

        return Double.parseDouble(s.substring(i + 1));
    }

    /**
     * Sets the tip heights from the tip dates
     */
    public static void setHeightsFromDates(MutableTree tree) {

        Date mostRecent = null;

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {

            Taxon taxon = tree.getNodeTaxon(tree.getExternalNode(i));

            Date date = (Date) taxon.getAttribute("date");

            if (date != null) {

                if ((mostRecent == null) || date.after(mostRecent)) {
                    mostRecent = date;
                }
            }
        }

        TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);
            Taxon taxon = tree.getNodeTaxon(node);
            Date date = (Date) taxon.getAttribute("date");

            if (date != null) {
                double height = timeScale.convertTime(date.getTimeValue(), date);
                tree.setNodeHeight(node, height);
            } else {
                tree.setNodeHeight(node, 0.0);
            }
        }

        adjustInternalHeights(tree, tree.getRoot());

        if (mostRecent != null) {
            tree.setUnits(mostRecent.getUnits());
        }
    }

    // **************************************************************
    // Private static methods
    // **************************************************************

    private static void adjustInternalHeights(MutableTree tree, NodeRef node) {

        if (!tree.isExternal(node)) {
            // pre-order recursion
            for (int i = 0; i < tree.getChildCount(node); i++) {
                adjustInternalHeights(tree, tree.getChild(node, i));
            }
        }

        NodeRef parent = tree.getParent(node);

        if (parent != null) {

            if (tree.getNodeHeight(parent) < tree.getNodeHeight(node)) {
                tree.setNodeHeight(parent, tree.getNodeHeight(node));
            }
        }
    }


}
