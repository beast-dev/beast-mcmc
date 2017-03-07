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

import java.io.IOException;
import java.util.ArrayList;

/**
 * @author Guy Baele
 * Path difference metric according to Kendall & Colijn (2015)
 */
public class KCPathDifferenceMetric {

    private Tree focalTree;
    private int dim, externalNodeCount;
    private double[] focalSmallM, focalLargeM;

    public KCPathDifferenceMetric() {

    }

    public KCPathDifferenceMetric(Tree focalTree) {
        this.focalTree = focalTree;
        this.externalNodeCount = focalTree.getExternalNodeCount();
        this.dim = (externalNodeCount-2)*(externalNodeCount-1)+externalNodeCount;
        this.focalSmallM = new double[dim];
        this.focalLargeM = new double[dim];

        int index = 0;
        for (int i = 0; i < externalNodeCount; i++) {
            for (int j = i+1; j < externalNodeCount; j++) {
                //get two leaf nodes
                NodeRef nodeOne = focalTree.getExternalNode(i);
                NodeRef nodeTwo = focalTree.getExternalNode(j);

                //get common ancestor of 2 leaf nodes
                NodeRef MRCA = Tree.Utils.getCommonAncestor(focalTree, nodeOne, nodeTwo);

                int edges = 0;
                double branchLengths = 0.0;
                while (MRCA != focalTree.getRoot()) {
                    edges++;
                    branchLengths += focalTree.getNodeHeight(focalTree.getParent(MRCA)) - focalTree.getNodeHeight(MRCA);
                    MRCA = focalTree.getParent(MRCA);
                }
                focalSmallM[index] = edges;
                focalLargeM[index] = branchLengths;
                index++;
            }
        }

        //fill out arrays further
        index = 0;
        for (int i = (externalNodeCount-1)*(externalNodeCount-2); i < dim; i++) {
            focalSmallM[i] = 1.0;
            focalLargeM[i] = focalTree.getNodeHeight(focalTree.getParent(focalTree.getExternalNode(index))) - focalTree.getNodeHeight(focalTree.getExternalNode(index));
            index++;
        }
    }

    public ArrayList<Double> getMetric(Tree tree, ArrayList<Double> lambda) {

        //check if taxon lists are in the same order!!
        if (focalTree.getExternalNodeCount() != tree.getExternalNodeCount()) {
            throw new RuntimeException("Different number of taxa in both trees.");
        } else {
            for (int i = 0; i < focalTree.getExternalNodeCount(); i++) {
                if (!focalTree.getNodeTaxon(focalTree.getExternalNode(i)).getId().equals(tree.getNodeTaxon(tree.getExternalNode(i)).getId())) {
                    throw new RuntimeException("Mismatch between taxa in both trees: " + focalTree.getNodeTaxon(focalTree.getExternalNode(i)).getId() + " vs. " + tree.getNodeTaxon(tree.getExternalNode(i)).getId());
                }
            }
        }

        double[] smallMTwo = new double[dim];
        double[] largeMTwo = new double[dim];

        int index = 0;
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            for (int j = i+1; j < tree.getExternalNodeCount(); j++) {
                //get two leaf nodes
                NodeRef nodeOne = tree.getExternalNode(i);
                NodeRef nodeTwo = tree.getExternalNode(j);

                //get common ancestor of 2 leaf nodes
                NodeRef MRCA = Tree.Utils.getCommonAncestor(tree, nodeOne, nodeTwo);

                int edges = 0;
                double branchLengths = 0.0;
                while (MRCA != tree.getRoot()) {
                    edges++;
                    branchLengths += tree.getNodeHeight(tree.getParent(MRCA)) - tree.getNodeHeight(MRCA);
                    MRCA = tree.getParent(MRCA);
                }
                smallMTwo[index] = edges;
                largeMTwo[index] = branchLengths;
                index++;
            }
        }

        //fill out arrays further
        index = 0;
        for (int i = (tree.getExternalNodeCount()-1)*(tree.getExternalNodeCount()-2); i < dim; i++) {
            smallMTwo[i] = 1.0;
            largeMTwo[i] = tree.getNodeHeight(tree.getParent(tree.getExternalNode(index))) - tree.getNodeHeight(tree.getExternalNode(index));
            index++;
        }

        double[] vArrayOne = new double[dim];
        double[] vArrayTwo = new double[dim];

        ArrayList<Double> results = new ArrayList<Double>();

        for (Double l : lambda) {
            double distance = 0.0;
            //calculate Euclidean distance for this lambda value
            for (int i = 0; i < dim; i++) {
                vArrayOne[i] = (1.0 - l)*focalSmallM[i] + l*focalLargeM[i];
                vArrayTwo[i] = (1.0 - l)*smallMTwo[i] + l*largeMTwo[i];
                distance += Math.pow(vArrayOne[i] - vArrayTwo[i],2);
            }
            distance = Math.sqrt(distance);
            results.add(distance);
        }

        return results;

    }

    /**
     * This method bypasses the constructor entirely, computing the metric on the two provided trees
     * and ignoring the internally stored tree.
     * @param tree1 Focal tree that will be used for computing the metric
     * @param tree2 Provided tree that will be compared to the focal tree
     * @param lambda Collection of lambda values for which to compute the metric
     * @return
     */
    public ArrayList<Double> getMetric(Tree tree1, Tree tree2, ArrayList<Double> lambda) {

        int dim = (tree1.getExternalNodeCount()-2)*(tree1.getExternalNodeCount()-1)+tree1.getExternalNodeCount();

        double[] smallMOne = new double[dim];
        double[] largeMOne = new double[dim];
        double[] smallMTwo = new double[dim];
        double[] largeMTwo = new double[dim];

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

                int edges = 0;
                double branchLengths = 0.0;
                while (MRCA != tree1.getRoot()) {
                    edges++;
                    branchLengths += tree1.getNodeHeight(tree1.getParent(MRCA)) - tree1.getNodeHeight(MRCA);
                    MRCA = tree1.getParent(MRCA);
                }
                smallMOne[index] = edges;
                largeMOne[index] = branchLengths;
                index++;
            }
        }

        //fill out arrays further
        index = 0;
        for (int i = (tree1.getExternalNodeCount()-1)*(tree1.getExternalNodeCount()-2); i < dim; i++) {
            smallMOne[i] = 1.0;
            largeMOne[i] = tree1.getNodeHeight(tree1.getParent(tree1.getExternalNode(index))) - tree1.getNodeHeight(tree1.getExternalNode(index));
            index++;
        }

        /*for (int i = 0; i < smallMOne.length; i++) {
            System.out.print(smallMOne[i] + " ");
        }
        System.out.println();
        for (int i = 0; i < largeMOne.length; i++) {
            System.out.print(largeMOne[i] + " ");
        }
        System.out.println("\n");*/

        index = 0;
        for (int i = 0; i < tree2.getExternalNodeCount(); i++) {
            for (int j = i+1; j < tree2.getExternalNodeCount(); j++) {
                //get two leaf nodes
                NodeRef nodeOne = tree2.getExternalNode(i);
                NodeRef nodeTwo = tree2.getExternalNode(j);

                //get common ancestor of 2 leaf nodes
                NodeRef MRCA = Tree.Utils.getCommonAncestor(tree2, nodeOne, nodeTwo);

                int edges = 0;
                double branchLengths = 0.0;
                while (MRCA != tree2.getRoot()) {
                    edges++;
                    branchLengths += tree2.getNodeHeight(tree2.getParent(MRCA)) - tree2.getNodeHeight(MRCA);
                    MRCA = tree2.getParent(MRCA);
                }
                smallMTwo[index] = edges;
                largeMTwo[index] = branchLengths;
                index++;
            }
        }

        //fill out arrays further
        index = 0;
        for (int i = (tree2.getExternalNodeCount()-1)*(tree2.getExternalNodeCount()-2); i < dim; i++) {
            smallMTwo[i] = 1.0;
            largeMTwo[i] = tree2.getNodeHeight(tree2.getParent(tree2.getExternalNode(index))) - tree2.getNodeHeight(tree2.getExternalNode(index));
            index++;
        }

        /*for (int i = 0; i < smallMTwo.length; i++) {
            System.out.print(smallMTwo[i] + " ");
        }
        System.out.println();
        for (int i = 0; i < largeMTwo.length; i++) {
            System.out.print(largeMTwo[i] + " ");
        }
        System.out.println("\n");*/

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

    public static void main(String[] args) {

        try {
            NewickImporter importer = new NewickImporter("(('A':1.2,'B':0.8):0.5,('C':0.8,'D':1.0):1.1)");
            Tree treeOne = importer.importNextTree();
            System.out.println("tree 1: " + treeOne);

            importer = new NewickImporter("((('A':0.8,'B':1.4):0.3,'C':0.7):0.9,'D':1.0)");
            Tree treeTwo = importer.importNextTree();
            System.out.println("tree 2: " + treeTwo + "\n");

            ArrayList<Double> lambdaValues = new ArrayList<Double>();
            lambdaValues.add(0.0);
            lambdaValues.add(0.5);
            lambdaValues.add(1.0);
            ArrayList<Double> metric = (new KCPathDifferenceMetric().getMetric(treeOne, treeTwo, lambdaValues));

            System.out.println("lambda (0.0) = " + metric.get(0));
            System.out.println("lambda (0.5) = " + metric.get(1));
            System.out.println("lambda (1.0) = " + metric.get(2));


            //Additional test for comparing a collection of trees against a (fixed) focal tree
            KCPathDifferenceMetric focalMetric = new KCPathDifferenceMetric(treeOne);
            metric = focalMetric.getMetric(treeTwo, lambdaValues);

            System.out.println("lambda (0.0) = " + metric.get(0));
            System.out.println("lambda (0.5) = " + metric.get(1));
            System.out.println("lambda (1.0) = " + metric.get(2));


        } catch(Importer.ImportException ie) {
            System.err.println(ie);
        } catch(IOException ioe) {
            System.err.println(ioe);
        }

    }

}
