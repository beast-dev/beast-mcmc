/*
 * HIPSTRTreeBuilder.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.tools.newtreeannotator;

import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.util.Pair;

import java.util.HashMap;
import java.util.Map;

public class HIPSTRTreeBuilder {
    private Map<Clade, Double> credibilityCache = new HashMap<>();

    public Tree getHIPSTRTree(CladeSystem cladeSystem) {
        BiClade rootClade = (BiClade)cladeSystem.getRootClade();

        credibilityCache.clear();

        score = findHIPSTRTree(cladeSystem, rootClade);

        SimpleTree tree = new SimpleTree(buildHIPSTRTree(cladeSystem, rootClade));

        return tree;
    }

    private double findHIPSTRTree(CladeSystem cladeSystem, BiClade clade) {

        double logCredibility = Math.log(clade.getCredibility());

        if (clade.getSize() > 1) {
            double bestLogCredibility = Double.NEGATIVE_INFINITY;

            for (Pair<BiClade, BiClade> subClade : clade.getSubClades()) {

                BiClade left = (BiClade)cladeSystem.getClade(subClade.first.getKey());
                if (left == null) {
                    throw new IllegalArgumentException("no clade found");
                }

                double leftLogCredibility = credibilityCache.getOrDefault(left, Double.NaN);
                if (Double.isNaN(leftLogCredibility)) {
                    leftLogCredibility = findHIPSTRTree(cladeSystem, left);
                    credibilityCache.put(left, leftLogCredibility);
                }
                BiClade right = (BiClade)cladeSystem.getClade(subClade.second.getKey());
                if (right == null) {
                    throw new IllegalArgumentException("no clade found");
                }
                double rightLogCredibility = credibilityCache.getOrDefault(right, Double.NaN);
                if (Double.isNaN(rightLogCredibility)) {
                    rightLogCredibility = findHIPSTRTree(cladeSystem, right);
                    credibilityCache.put(right, rightLogCredibility);
                }

                if (leftLogCredibility + rightLogCredibility > bestLogCredibility) {
                    bestLogCredibility = leftLogCredibility + rightLogCredibility;
                    clade.bestLeft = left;
                    clade.bestRight = right;
                }
            }

            logCredibility += bestLogCredibility;
            clade.bestSubTreeCredibility = logCredibility;

        }

        return logCredibility;
    }

    private SimpleNode buildHIPSTRTree(CladeSystem cladeSystem, Clade clade) {
        SimpleNode newNode = new SimpleNode();
        if (clade.getSize() == 1) {
            newNode.setTaxon(clade.getTaxon());
        } else {
            newNode.addChild(buildHIPSTRTree(cladeSystem, clade.getBestLeft()));
            newNode.addChild(buildHIPSTRTree(cladeSystem, clade.getBestRight()));
        }
        return newNode;
    }

    public double getScore() {
        return score;
    }

    private double score = Double.NaN;
}
