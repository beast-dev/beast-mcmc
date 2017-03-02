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

import java.util.ArrayList;

/**
 * @author Guy Baele
 * Path difference metric according to Kendall & Colijn (2015)
 */
public class KCPathDifferenceMetric {

    public KCPathDifferenceMetric() {

    }

    public ArrayList<Double> getMetric(Tree tree1, Tree tree2, ArrayList<Double> lambda) {

        int dim = tree1.getExternalNodeCount()*tree1.getExternalNodeCount();

        double[] smallMOne = new double[dim];
        double[] largeMOne = new double[dim];
        double[] smallMTwo = new double[dim];
        double[] largeMTwo = new double[dim];

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
            for (int j = i+1; j < tree1.getExternalNodeCount(); j++) {
                //get two leaf nodes
                NodeRef nodeOne = tree1.getExternalNode(i);
                NodeRef nodeTwo = tree1.getExternalNode(j);

                //get common ancestor of 2 leaf nodes
                NodeRef MRCA = Tree.Utils.getCommonAncestor(tree1, nodeOne, nodeTwo);

                int edges = 0;
                double branchLengths = 0.0;
                while (MRCA != tree1.getRoot()) {
                    edges++;
                    branchLengths += tree1.getNodeHeight(tree1.getRoot()) - tree1.getNodeHeight(MRCA);
                    MRCA = tree1.getParent(MRCA);
                }
                smallMOne[index] = edges;
                largeMOne[index] = branchLengths;
                index++;
            }
        }

        //fill out arrays further
        index = 0;
        for (int i = tree1.getExternalNodeCount()*(tree1.getExternalNodeCount()-1); i < dim; i++) {
            smallMOne[i] = 1.0;
            largeMOne[i] = tree1.getNodeHeight(tree1.getExternalNode(index));
            index++;
        }

        index = 0;
        for (int i = 0; i < tree2.getExternalNodeCount(); i++) {
            for (int j = i; j < tree2.getExternalNodeCount(); j++) {
                //get two leaf nodes
                NodeRef nodeOne = tree2.getExternalNode(i);
                NodeRef nodeTwo = tree2.getExternalNode(j);

                //get common ancestor of 2 leaf nodes
                NodeRef MRCA = Tree.Utils.getCommonAncestor(tree2, nodeOne, nodeTwo);

                int edges = 0;
                double branchLengths = 0.0;
                while (MRCA != tree2.getRoot()) {
                    edges++;
                    branchLengths += tree2.getNodeHeight(tree2.getRoot()) - tree2.getNodeHeight(MRCA);
                    MRCA = tree2.getParent(MRCA);
                }
                smallMTwo[index] = edges;
                largeMTwo[index] = branchLengths;
                index++;
            }
        }

        //fill out arrays further
        index = 0;
        for (int i = tree2.getExternalNodeCount()*(tree2.getExternalNodeCount()-1); i < dim; i++) {
            smallMTwo[i] = 1.0;
            largeMTwo[i] = tree2.getNodeHeight(tree2.getExternalNode(index));
            index++;
        }

        double[] vArrayOne = new double[dim];
        double[] vArrayTwo = new double[dim];

        ArrayList<Double> results = new ArrayList<Double>();

        for (Double l : lambda) {
            double distance = 0.0;
            //calculate Euclidean distance for this lambda value
            for (int i = 0; i < dim; i++) {
                vArrayOne[i] = (1.0 - l)*smallMOne[i] + l*largeMOne[i];
                vArrayTwo[i] = (1.0 - l)*smallMTwo[i] + l*largeMTwo[i];
                distance += Math.pow(vArrayOne[i] - vArrayTwo[i],2);
            }
            distance = Math.sqrt(distance);
            results.add(distance);
        }

        return results;
    }

}
