/*
 * UPGMATreeParser.java
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
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.UPGMATree;
import dr.evolution.util.TimeScale;
import dr.evomodelxml.TreeModelParser;
import dr.xml.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: UPGMATreeParser.java,v 1.6 2006/07/28 11:27:32 rambaut Exp $
 */
public class UPGMATreeParser extends AbstractXMLObjectParser {

    //
    // Public stuff
    //

    public static final String UPGMA_TREE = "upgmaTree";
    //public static final String DISTANCES = "distances";
    public static final String ROOT_HEIGHT = TreeModelParser.ROOT_HEIGHT;

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


        return tree;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(SimpleTreeParser.USING_DATES, true),
            AttributeRule.newDoubleRule(ROOT_HEIGHT, true),
            new ElementRule(DistanceMatrix.class)
    };

    public String getParserDescription() {
        return "This element returns a UPGMA tree generated from the given distances.";
    }

    public Class getReturnType() {
        return UPGMATree.class;
    }
}
