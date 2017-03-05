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

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import jebl.evolution.trees.RootedTree;

import java.io.IOException;
import java.util.ArrayList;

/**
 * @author Guy Baele
 * Path difference metric according to Steel & Penny (1993)
 */
public class SPPathDifferenceMetric {

    public SPPathDifferenceMetric() {

    }

    public double getMetric(Tree tree1, Tree tree2) {

        int dim = (tree1.getExternalNodeCount()-2)*(tree1.getExternalNodeCount()-1);

        double[] pathOne = new double[dim];
        double[] pathTwo = new double[dim];

        //check if taxon lists are in the same order!!
        if (tree1.getExternalNodeCount() != tree2.getExternalNodeCount()) {
            throw new RuntimeException("Different number of taxa in both trees.");
        } else {
            for (int i = 0; i < tree1.getExternalNodeCount(); i++) {
                if (!tree1.getNodeTaxon(tree1.getExternalNode(i)).getId().equals(tree2.getNodeTaxon(tree2.getExternalNode(i)).getId())) {
                    throw new RuntimeException("Mismatch between taxa in both trees: " + tree1.getNodeTaxon(tree1.getExternalNode(i)).getId() + " vs. " + tree2.getNodeTaxon(tree2.getExternalNode(i)).getId());
                }
            }
        }

        int index = 0;
        for (int i = 0; i < tree1.getExternalNodeCount(); i++) {
            for (int j = i+1; j < tree1.getExternalNodeCount(); j++) {
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

        /*for (int i = 0; i < pathOne.length; i++) {
            System.out.print(pathOne[i] + " ");
        }
        System.out.println();*/

        index = 0;
        for (int i = 0; i < tree2.getExternalNodeCount(); i++) {
            for (int j = i+1; j < tree2.getExternalNodeCount(); j++) {
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

        /*for (int i = 0; i < pathTwo.length; i++) {
            System.out.print(pathTwo[i] + " ");
        }
        System.out.println();*/

        double metric = 0.0;
        for (int i = 0; i < dim; i++) {
            metric += Math.pow(pathOne[i] - pathTwo[i],2);
        }
        metric = Math.sqrt(metric);

        return metric;
    }

    public static void main(String[] args) {

        try {
            NewickImporter importer = new NewickImporter("(('A':1.2,'B':0.8):0.5,('C':0.8,'D':1.0):1.1)");
            Tree treeOne = importer.importNextTree();
            System.out.println("tree 1: " + treeOne);

            importer = new NewickImporter("((('A':0.8,'B':1.4):0.3,'C':0.7):0.9,'D':1.0)");
            Tree treeTwo = importer.importNextTree();
            System.out.println("tree 2: " + treeTwo + "\n");

            double metric = (new SPPathDifferenceMetric().getMetric(treeOne, treeTwo));

            System.out.println("path difference = " + metric);

        } catch(Importer.ImportException ie) {
            System.err.println(ie);
        } catch(IOException ioe) {
            System.err.println(ioe);
        }

    }

}
