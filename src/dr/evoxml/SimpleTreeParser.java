/*
 * SimpleTreeParser.java
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

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.TimeScale;
import dr.evolution.util.Units;
import dr.evoxml.util.XMLUnits;
import dr.xml.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: SimpleTreeParser.java,v 1.3 2005/05/24 20:25:59 rambaut Exp $
 */
public class SimpleTreeParser extends AbstractXMLObjectParser {

    //
    // Public stuff
    //

    public static final String SIMPLE_TREE = "tree";
    public static final String USING_DATES = "usingDates";
    public static final String DATE = "date";

    public String getParserName() {
        return SIMPLE_TREE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);
        boolean usingDates = xo.getAttribute(USING_DATES, false);

        SimpleNode root = (SimpleNode) xo.getChild(SimpleNode.class);

        if (root == null) {
            throw new XMLParseException("node element missing from tree");
        }

        SimpleTree tree = new SimpleTree(root);

        if (usingDates) {

            dr.evolution.util.Date mostRecent = null;
            for (int i = 0; i < tree.getTaxonCount(); i++) {

                dr.evolution.util.Date date = (dr.evolution.util.Date) tree.getTaxonAttribute(i, DATE);

                if (date == null) {
                    date = (dr.evolution.util.Date) tree.getNodeAttribute(tree.getExternalNode(i), DATE);
                }

                if (date != null && ((mostRecent == null) || date.after(mostRecent))) {
                    mostRecent = date;
                }
            }

            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                dr.evolution.util.Date date = (dr.evolution.util.Date) tree.getNodeAttribute(tree.getInternalNode(i), DATE);

                if (date != null && ((mostRecent == null) || date.after(mostRecent))) {
                    mostRecent = date;
                }
            }

            if (mostRecent == null) {
                throw new XMLParseException("no date elements in tree (and usingDates attribute set)");
            }

            TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());

            for (int i = 0; i < tree.getTaxonCount(); i++) {
                dr.evolution.util.Date date = (dr.evolution.util.Date) tree.getTaxonAttribute(i, DATE);

                if (date == null) {
                    date = (dr.evolution.util.Date) tree.getNodeAttribute(tree.getExternalNode(i), DATE);
                }

                if (date != null) {
                    double height = timeScale.convertTime(date.getTimeValue(), date);
                    tree.setNodeHeight(tree.getExternalNode(i), height);
                }
            }

            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                dr.evolution.util.Date date = (dr.evolution.util.Date) tree.getNodeAttribute(tree.getInternalNode(i), DATE);

                if (date != null) {
                    double height = timeScale.convertTime(date.getTimeValue(), date);
                    tree.setNodeHeight(tree.getInternalNode(i), height);
                }
            }

            MutableTree.Utils.correctHeightsForTips(tree);
        }

        tree.setUnits(units);

        return tree;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a rooted binary tree and associated attributes.";
    }

    public Class getReturnType() {
        return Tree.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            XMLUnits.SYNTAX_RULES[0],
            new ElementRule(SimpleNode.class)
    };
}
