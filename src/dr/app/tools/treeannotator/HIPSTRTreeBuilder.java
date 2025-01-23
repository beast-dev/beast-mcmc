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

package dr.app.tools.treeannotator;

import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.util.Pair;

import java.util.HashMap;
import java.util.Map;

public class HIPSTRTreeBuilder {
    static final boolean breakTies = false;

    private final Map<Clade, Double> credibilityCache = new HashMap<>();

    public MutableTree getHIPSTRTree(CladeSystem cladeSystem, TaxonList taxonList, double penaltyThreshold) {
        BiClade rootClade = (BiClade)cladeSystem.getRootClade();

        credibilityCache.clear();

        cladeSystem.embiggenBiClades(1, 5);
//        cladeSystem.embiggenBiCladesThreaded(1, 5);

        score = findHIPSTRTree(rootClade, penaltyThreshold);

        // create a map so that tip numbers are in the same order as the taxon list
        Map<Taxon, Integer> taxonNumberMap = new HashMap<>();
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            taxonNumberMap.put(taxonList.getTaxon(i), i);
        }

        FlexibleTree tree = new FlexibleTree(buildHIPSTRTree(rootClade), taxonNumberMap);
        tree.setLengthsKnown(false);
        return tree;
    }

    private double findHIPSTRTree(BiClade clade, double penaltyThreshold) {

        assert clade.getSize() > 1;

        double cladeScore = Math.log(clade.getCredibility() + penaltyThreshold);
//        double cladeScore = Math.log(clade.getCredibility() + penaltyThreshold) + clade.getCredibility() >= 0.5 ? 1.0E10 : 0.0;

        Map<Pair<BiClade, BiClade>, Double> ties;
        if (breakTies) {
            ties = new HashMap<>();
        }

        // if it is not a tip
        if (clade.getSize() > 2) {
            // more than two tips in this clade
            double bestSubtreeScore = Double.NEGATIVE_INFINITY;

            for (Pair<BiClade, BiClade> subClade : clade.getSubClades()) {
                BiClade left = subClade.first;

                double leftScore = Math.log(1.0 + penaltyThreshold);
//                double leftScore = 0.0);
                if (left.getSize() > 1) {
                    leftScore = credibilityCache.getOrDefault(left, Double.NaN);
                    if (Double.isNaN(leftScore)) {
                        leftScore = findHIPSTRTree(left, penaltyThreshold);
                        credibilityCache.put(left, leftScore);
                    }
                }

                BiClade right = subClade.second;

                double rightScore = Math.log(1.0 + penaltyThreshold);
//                double rightScore = 0.0;
                if (right.getSize() > 1) {
                    rightScore = credibilityCache.getOrDefault(right, Double.NaN);
                    if (Double.isNaN(rightScore)) {
                        rightScore = findHIPSTRTree(right, penaltyThreshold);
                        credibilityCache.put(right, rightScore);
                    }
                }

                if (breakTies) {
                    if (leftScore + rightScore >= bestSubtreeScore) {
                        if (leftScore + rightScore > bestSubtreeScore) {
                            ties.clear();
                            bestSubtreeScore = leftScore + rightScore;
                            if (left.getSize() < right.getSize()) {
                                // sort by clade size because why not...
                                clade.bestLeft = left;
                                clade.bestRight = right;
                            } else {
                                clade.bestLeft = right;
                                clade.bestRight = left;
                            }
                            ties.put(new Pair<>(left, right), (left.getCredibility() + penaltyThreshold) * (right.getCredibility() + penaltyThreshold));
                        }
                    }
                } else {
                    if (leftScore + rightScore > bestSubtreeScore) {
                        bestSubtreeScore = leftScore + rightScore;
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

                if (breakTies && ties.size() > 1) {
                    double best = Double.NEGATIVE_INFINITY;
                    for (Pair<BiClade, BiClade> key : ties.keySet()) {
                        double value = ties.get(key);
                        if (value > best) {
                            clade.bestLeft = key.first;
                            clade.bestRight = key.second;
                            best = value;
                        }
                    }
                }
            }

            cladeScore = cladeScore + bestSubtreeScore;
        } else {
            // two tips so there will only be one pair and their sum log cred will be 0.0
            assert clade.getSubClades().size() == 1;
            Pair<BiClade, BiClade> subClade = clade.getSubClades().stream().findFirst().get();
            clade.bestLeft = subClade.first;
            clade.bestRight = subClade.second;
            cladeScore += 2 * Math.log(1.0 + penaltyThreshold);
        }

        clade.bestSubTreeScore = cladeScore;

        return cladeScore;
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
