/*
 * TreeMetricStatisticParser.java
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

import dr.evolution.tree.Tree;
import dr.evolution.tree.treemetrics.*;
import dr.evomodel.tree.TreeMetricStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class TreeMetricStatisticParser extends AbstractXMLObjectParser {

    public static final String TREE_METRIC_STATISTIC = "treeMetricStatistic";
    public static final String TARGET = "target";
    public static final String REFERENCE = "reference";
    public static final String TYPE = "type";

    public String getParserName() {
        return TREE_METRIC_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeMetric.Type type = null;
        if (xo.hasAttribute(TYPE)) {
            final String s = xo.getStringAttribute(TYPE);
            for (TreeMetric.Type t : TreeMetric.Type.values()) {
                if (t.toString().toLowerCase().equals(s.toLowerCase())) {
                    type = t;
                    break;
                }
            }

            if (type == null) {
                throw new XMLParseException("Tree metric type, " + s + ", is not recognised");
            }
        }

        final String name = xo.getAttribute(Statistic.NAME, xo.hasId() ? xo.getId() :
                (type == null ? "topology" : type.toString()));
        final Tree target = (Tree) xo.getElementFirstChild(TARGET);
        final Tree reference = (Tree) xo.getElementFirstChild(REFERENCE);

        TreeMetric treeMetric = null;
        if (type != null) {
            switch (type) {
                case ROBINSON_FOULDS:
                    treeMetric = new RobinsonFouldsMetric();
                    break;
                case BRANCH_SCORE:
                    treeMetric = new BranchScoreMetric();
                    break;
                case ROOTED_BRANCH_SCORE:
                    treeMetric = new RootedBranchScoreMetric();
                    break;
                case CLADE_HEIGHT:
                    treeMetric = new CladeHeightMetric();
                    break;
                case KENDALL_COLIJN:
                    treeMetric = new KendallColijnPathDifferenceMetric(0.5);
                    break;
                case STEEL_PENNY:
                    treeMetric = new SteelPennyPathDifferenceMetric();
                    break;
                default:
                    throw new IllegalArgumentException("Unknown tree metric type");
            }
        }

        return new TreeMetricStatistic(name, reference, target, treeMetric);
    }

    // ************************************************************************
    // AbstractXMLObjectParser
    // implementation
    // ************************************************************************

    public String getParserDescription() {
        return "A statistic that returns the distance between two trees. "
                + " with method=\"topology\", return a 0 for identity and a 1 for difference. "
                + "With other methods return the distance metric associated with that method.";
    }

    public Class getReturnType() {
        return TreeMetricStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new StringAttributeRule(Statistic.NAME,
                    "A name for this statistic primarily for the purposes of logging",
                    true),
            new StringAttributeRule(TYPE, "tree metric name", true),
            new ElementRule(TARGET, new XMLSyntaxRule[]{new ElementRule(
                    Tree.class)}),
            new ElementRule(REFERENCE, new XMLSyntaxRule[]{new ElementRule(
                    Tree.class)}),
    };
}
