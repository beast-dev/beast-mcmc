/*
 * TemporalStress.java
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

package dr.app.pathogen;

import java.util.*;

import dr.evolution.tree.TreeUtils;
import dr.stats.DiscreteStatistics;
import dr.math.UnivariateFunction;
import dr.math.UnivariateMinimum;
import dr.evolution.tree.*;
import dr.evolution.util.Taxon;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TemporalStress {

    public static Set<Taxon> annotateStress(MutableTree tree, NodeRef node) {
        Set<Taxon> taxa = new HashSet<Taxon>();

        if (!tree.isExternal(node)) {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                taxa.addAll(annotateStress(tree, child));
            }

            if (taxa.size() > 2) {
                Tree subtree = null;

                double stress = findGlobalRoot(subtree);
                tree.setNodeAttribute(node, "stress", stress);

            }
        } else {
            taxa.add(tree.getNodeTaxon(node));
        }

        return taxa;
    }


    private static double findGlobalRoot(Tree source) {

        FlexibleTree bestTree = new FlexibleTree(source);
        double minF = findLocalRoot(bestTree);

        for (int i = 0; i < source.getNodeCount(); i++) {
            FlexibleTree tmpTree = new FlexibleTree(source);
            NodeRef node = tmpTree.getNode(i);
            if (!tmpTree.isRoot(node)) {
                double length = tmpTree.getBranchLength(node);
                tmpTree.changeRoot(node, length * 0.5, length * 0.5);

                double f = findLocalRoot(tmpTree);
                if (f < minF) {
                    minF = f;
                    bestTree = tmpTree;
                }
            }
        }
        return minF;
    }

    private static double findLocalRoot(final FlexibleTree tree) {

        NodeRef node1 = tree.getChild(tree.getRoot(), 0);
        NodeRef node2 = tree.getChild(tree.getRoot(), 1);

        final double length1 = tree.getBranchLength(node1);
        final double length2 = tree.getBranchLength(node2);

        final double sumLength = length1 + length2;

        final Set<NodeRef> tipSet1 = TreeUtils.getExternalNodes(tree, node1);
        final Set<NodeRef> tipSet2 = TreeUtils.getExternalNodes(tree, node2);

        final double[] y = new double[tree.getExternalNodeCount()];

        UnivariateFunction f = new UnivariateFunction() {
            public double evaluate(double argument) {
                double l1 = argument * sumLength;

                for (NodeRef tip : tipSet1) {
                    y[tip.getNumber()] = getRootToTipDistance(tree, tip) - length1 + l1;
                }

                double l2 = (1.0 - argument) * sumLength;

                for (NodeRef tip : tipSet2) {
                    y[tip.getNumber()] = getRootToTipDistance(tree, tip) - length2 + l2;
                }

                return DiscreteStatistics.variance(y);
            }

            public double getLowerBound() { return 0; }
            public double getUpperBound() { return 1.0; }
        };

        UnivariateMinimum minimum = new UnivariateMinimum();

        double x = minimum.findMinimum(f);

        double fminx = minimum.fminx;

        double l1 = x * sumLength;
        double l2 = (1.0 - x) * sumLength;

        tree.setBranchLength(node1, l1);
        tree.setBranchLength(node2, l2);

        return fminx;
    }

    private static double getRootToTipDistance(Tree tree, NodeRef node) {
        double distance = 0;
        while (node != null) {
            distance += tree.getBranchLength(node);
            node = tree.getParent(node);
        }
        return distance;
    }


}
