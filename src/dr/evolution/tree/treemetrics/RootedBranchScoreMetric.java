/*
 * RootedBranchScoreMetric.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

/**
 *
 */
package dr.evolution.tree.treemetrics;

import dr.evolution.tree.Clade;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * @author Guy Baele
 *
 */
public class RootedBranchScoreMetric extends BranchScoreMetric implements Citable {

    public static Type TYPE = Type.ROOTED_BRANCH_SCORE;

    public RootedBranchScoreMetric() {
    }

    @Override
    public double getMetric(Tree tree1, Tree tree2) {

        TreeMetric.Utils.checkTreeTaxa(tree1, tree2);

        List<Clade> clades1 = Clade.getCladeList(tree1);
        List<Clade> clades2 = Clade.getCladeList(tree2);

        return Math.sqrt(Math.pow(getDistance(clades1, clades2), 2) + getExternalDistance(tree1, tree2));
    }

    private double getExternalDistance(Tree tree1, Tree tree2) {

        double distance = 0.0;

        for (int i = 0; i < tree1.getExternalNodeCount(); i++) {
            NodeRef node1 = tree1.getExternalNode(i);
            NodeRef node2 = tree2.getExternalNode(i);
            NodeRef parent1 = tree1.getParent(node1);
            NodeRef parent2 = tree2.getParent(node2);
            distance += Math.pow(
                    (tree1.getNodeHeight(parent1) - tree1.getNodeHeight(node1))-
                            (tree2.getNodeHeight(parent2) - tree2.getNodeHeight(node2)), 2);
        }

        return distance;

    }

    @Override
    public Type getType() {
        return TYPE;
    }

    @Override
    public String toString() {
        return getType().getShortName();
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_METRICS;
    }

    @Override
    public String getDescription() {
        return "Normalized Rooted Branch Score";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("J", "Heled"),
                    new Author("AJ", "Drummond")
            },
            "Bayesian inference of species trees from multilocus data",
            2010,
            "Mol. Biol. Evol.",
            27,
            570, 580,
            Citation.Status.PUBLISHED
    );

}
