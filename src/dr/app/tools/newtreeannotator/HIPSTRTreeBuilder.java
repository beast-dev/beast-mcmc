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

import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.util.Pair;

import java.util.HashMap;
import java.util.Map;

public class HIPSTRTreeBuilder {
    private Map<Clade, Double> credibilityCache = new HashMap<>();

    public MutableTree getHIPSTRTree(CladeSystem cladeSystem, TaxonList taxonList) {
        BiClade rootClade = (BiClade)cladeSystem.getRootClade();

        credibilityCache.clear();

        score = findHIPSTRTree(rootClade);

        // create a map so that tip numbers are in the same order as the taxon list
        Map<Taxon, Integer> taxonNumberMap = new HashMap<>();
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            taxonNumberMap.put(taxonList.getTaxon(i), i);
        }

        FlexibleTree tree = new FlexibleTree(buildHIPSTRTree(rootClade), taxonNumberMap);
        tree.setLengthsKnown(false);
        return tree;
    }

    private double findHIPSTRTree(BiClade clade) {

        double logCredibility = Math.log(clade.getCredibility());

        if (clade.getSize() > 1) {
            double bestLogCredibility = Double.NEGATIVE_INFINITY;

            for (Pair<BiClade, BiClade> subClade : clade.getSubClades()) {
                BiClade left = subClade.first;

                double leftLogCredibility = credibilityCache.getOrDefault(left, Double.NaN);
                if (Double.isNaN(leftLogCredibility)) {
                    leftLogCredibility = findHIPSTRTree(left);
                    credibilityCache.put(left, leftLogCredibility);
                }

                BiClade right = subClade.second;

                double rightLogCredibility = credibilityCache.getOrDefault(right, Double.NaN);
                if (Double.isNaN(rightLogCredibility)) {
                    rightLogCredibility = findHIPSTRTree(right);
                    credibilityCache.put(right, rightLogCredibility);
                }

                if (leftLogCredibility + rightLogCredibility > bestLogCredibility) {
                    bestLogCredibility = leftLogCredibility + rightLogCredibility;
                    if (left.getSize() < right.getSize()) {
                        // sort by clade size because why not...
                        clade.bestLeft = left;
                        clade.bestRight = right;
                    } else {
                        clade.bestLeft = right;
                        clade.bestRight = left;
                    }
                }
            }

            logCredibility += bestLogCredibility;
            clade.bestSubTreeCredibility = logCredibility;

        }

        return logCredibility;
    }

    private FlexibleNode buildHIPSTRTree(Clade clade) {
        FlexibleNode newNode = new FlexibleNode();
        if (clade.getSize() == 1) {
            newNode.setTaxon(clade.getTaxon());
            newNode.setNumber(clade.getIndex());
        } else {
            newNode.addChild(buildHIPSTRTree(clade.getBestLeft()));
            newNode.addChild(buildHIPSTRTree(clade.getBestRight()));
            newNode.setNumber(-1);
        }
        return newNode;
    }

    public double getScore() {
        return score;
    }

    private double score = Double.NaN;
}
