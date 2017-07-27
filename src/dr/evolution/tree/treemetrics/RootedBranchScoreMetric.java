/*
 * RootedBranchScoreMetric.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Clade;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import jebl.evolution.graphs.Node;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.RootedTree;

import java.io.IOException;
import java.util.*;

/**
 * @author Guy Baele
 *
 */
public class RootedBranchScoreMetric extends BranchScoreMetric {

    public RootedBranchScoreMetric() {
    }

    @Override
    public double getMetric(Tree tree1, Tree tree2) {

        Utils.checkTreeTaxa(tree1, tree2);

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


    public static void main(String[] args) {

        try {

            NewickImporter importer = new NewickImporter("((('C':0.03365591238,'A':0.7225157402):0.306488578,'B':0.4572411443):0.4673149632,('D':0.7966438427,'E':0.8063645191):0.7478901469)");
            Tree treeOne = importer.importNextTree();
            System.out.println("tree 1: " + treeOne);

            importer = new NewickImporter("(('A':0.2333369483,'B':0.3468381313):0.5562255983,('C':0.8732210915,('D':0.9124725792,'E':0.1983703848):0.5252404297):0.2000638912)");
            Tree treeTwo = importer.importNextTree();
            System.out.println("tree 2: " + treeTwo + "\n");

            double metric = (new RootedBranchScoreMetric().getMetric(treeOne, treeTwo));

            System.out.println("rooted branch score metric = " + metric);

        } catch(Importer.ImportException ie) {
            System.err.println(ie);
        } catch(IOException ioe) {
            System.err.println(ioe);
        }

    }

}
