/*
 * TopologyTracer.java
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

package dr.evolution.tree;

import jebl.evolution.trees.RootedTree;

/**
 * @author Guy Baele
 * Path difference metric according to Steel & Penny (1993)
 */
public class SPPathDifferenceMetric {

    public SPPathDifferenceMetric() {

    }

    public double getMetric(Tree tree1, Tree tree2) {

        int dim = tree1.getExternalNodeCount()*(tree1.getExternalNodeCount()-1);

        double[] pathOne = new double[dim];
        double[] pathTwo = new double[dim];

        //check if taxon lists are in the same order!!
        if (tree1.getExternalNodeCount() != tree2.getExternalNodeCount()) {
            throw new RuntimeException("Different number of taxa in both trees.");
        } else {
            for (int i = 0; i < tree1.getExternalNodeCount(); i++) {
                if (tree1.getNodeTaxon(tree1.getExternalNode(i)) != tree2.getNodeTaxon(tree2.getExternalNode(i))) {
                    throw new RuntimeException("Mismatch between taxa in both trees: " + tree1.getNodeTaxon(tree1.getExternalNode(i)) + " vs. " + tree2.getNodeTaxon(tree2.getExternalNode(i)));
                }
            }
        }

        int index = 0;
        for (int i = 0; i < tree1.getExternalNodeCount(); i++) {
            for (int j = i; j < tree1.getExternalNodeCount(); j++) {
                //get two leaf nodes
                NodeRef nodeOne = tree1.getExternalNode(i);
                NodeRef nodeTwo = tree1.getExternalNode(j);

                //get common ancestor of 2 leaf nodes
                NodeRef MRCA = Tree.Utils.getCommonAncestor(tree1, nodeOne, nodeTwo);

                double pathLength = 0.0;
                while (nodeOne != MRCA) {
                    pathLength += tree1.getNodeHeight(tree1.getParent(nodeOne)) - tree1.getNodeHeight(nodeOne);
                    nodeOne = tree1.getParent(nodeOne);
                }
                while (nodeTwo != MRCA) {
                    pathLength += tree1.getNodeHeight(tree1.getParent(nodeTwo)) - tree1.getNodeHeight(nodeTwo);
                    nodeTwo = tree1.getParent(nodeTwo);
                }

                pathOne[index] = pathLength;
                index++;
            }
        }

        index = 0;
        for (int i = 0; i < tree2.getExternalNodeCount(); i++) {
            for (int j = i; j < tree2.getExternalNodeCount(); j++) {
                //get two leaf nodes
                NodeRef nodeOne = tree2.getExternalNode(i);
                NodeRef nodeTwo = tree2.getExternalNode(j);

                //get common ancestor of 2 leaf nodes
                NodeRef MRCA = Tree.Utils.getCommonAncestor(tree2, nodeOne, nodeTwo);

                double pathLength = 0.0;
                while (nodeOne != MRCA) {
                    pathLength += tree2.getNodeHeight(tree2.getParent(nodeOne)) - tree2.getNodeHeight(nodeOne);
                    nodeOne = tree2.getParent(nodeOne);
                }
                while (nodeTwo != MRCA) {
                    pathLength += tree2.getNodeHeight(tree2.getParent(nodeTwo)) - tree2.getNodeHeight(nodeTwo);
                    nodeTwo = tree2.getParent(nodeTwo);
                }

                pathTwo[index] = pathLength;
                index++;
            }
        }

        double metric = 0.0;
        for (int i = 0; i < dim; i++) {
            metric += Math.pow(pathOne[i] - pathTwo[i],2);
        }
        metric = Math.sqrt(metric);

        return metric;
    }

}
