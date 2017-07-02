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

package dr.evolution.tree.treemetrics;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;

import java.io.IOException;
import java.util.*;

import static dr.evolution.tree.treemetrics.TreeMetric.Utils.checkTreeTaxa;

/**
 * @author Guy Baele
 * Path difference metric according to Kendall & Colijn (2015)
 */
public class KendallColijnPathDifferenceMetric implements TreeMetric {

    public static Type TYPE = Type.KENDALL_COLIJN;

    private Tree focalTree;
    private int dim;
    private double[] focalSmallM, focalLargeM;
    private final boolean fixedFocalTree;
    private final double lambda;

    public KendallColijnPathDifferenceMetric(double lambda) {
        this.lambda = lambda;

        this.fixedFocalTree = false;
    }

    public KendallColijnPathDifferenceMetric(double lambda, Tree focalTree) {
        this.lambda = lambda;

        this.focalTree = focalTree;
        this.fixedFocalTree = true;
        //this.dim = (externalNodeCount-2)*(externalNodeCount-1)+externalNodeCount;
        this.dim = focalTree.getExternalNodeCount() * focalTree.getExternalNodeCount();
        this.focalSmallM = new double[dim];
        this.focalLargeM = new double[dim];

        traverse(focalTree, focalTree.getRoot(), 0.0, 0, focalLargeM, focalSmallM);
    }

    /**
     * Compute the metric between two trees. If tree1 is not the focal tree provided to the
     * constructor then it will store this as the new focal tree. If the focal tree is constant
     * from call to call then a cached set of precomputation will be used, increasing efficiency.
     * @param tree1
     * @param tree2
     * @return
     */
    @Override
    public double getMetric(Tree tree1, Tree tree2) {

        checkTreeTaxa(tree1, tree2);

        if (tree1 != focalTree) {
            if (fixedFocalTree) {
                // If we set a focal tree in the constructor then it makes sense to check it is the same
                // as the one set here.
                throw new RuntimeException("Focal tree is different from that set in the constructor.");
            }

            // cache tree1 and the pre-computed path for future calls
            focalTree = tree1;
            if (focalSmallM == null) {
                dim = focalTree.getExternalNodeCount() * focalTree.getExternalNodeCount();
                focalSmallM = new double[dim];
                focalLargeM = new double[dim];
            }
            traverse(focalTree, focalTree.getRoot(), 0.0, 0, focalLargeM, focalSmallM);
        }

        double[] smallMTwo = new double[dim];
        double[] largeMTwo = new double[dim];

        traverse(tree2, tree2.getRoot(), 0.0, 0, largeMTwo, smallMTwo);

        List<Double> results = new ArrayList<Double>();

        int n = tree1.getExternalNodeCount();

        return calculateMetric(focalSmallM, focalLargeM, smallMTwo, largeMTwo, n, lambda);
    }

    private double calculateMetric(double[] smallMOne, double[] largeMOne, double[] smallMTwo, double[] largeMTwo, int n, double l) {
        double distance = 0.0;
        //calculate Euclidean distance for this lambda value
        int k = 0;
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) { // include diagonal
                int index = (i * n) + j;
                double vOne = (1.0 - l) * smallMOne[index] + l * largeMOne[index];
                double vTwo = (1.0 - l) * smallMTwo[index] + l * largeMTwo[index];
                distance += Math.pow(vOne - vTwo, 2);
            }
        }
        return Math.sqrt(distance);
    }

    private Set<NodeRef> traverse(Tree tree, NodeRef node, double lengthFromRoot, int edgesFromRoot, double[] lengths, double[] edges) {
        NodeRef left = tree.getChild(node, 0);
        NodeRef right = tree.getChild(node, 1);


        Set<NodeRef> leftSet = null;
        Set<NodeRef> rightSet = null;

        if (!tree.isExternal(left)) {
            leftSet = traverse(tree, left, lengthFromRoot + tree.getBranchLength(left), edgesFromRoot + 1, lengths, edges);
        } else {
            leftSet = Collections.singleton(left);
            int index = (left.getNumber() * tree.getExternalNodeCount()) + left.getNumber();
            lengths[index] = tree.getBranchLength(left);
            edges[index] = 1;
        }
        if (!tree.isExternal(right)) {
            rightSet = traverse(tree, right, lengthFromRoot + tree.getBranchLength(right), edgesFromRoot + 1, lengths, edges);
        } else {
            rightSet = Collections.singleton(right);
            int index = (right.getNumber() * tree.getExternalNodeCount()) + right.getNumber();
            lengths[index] = tree.getBranchLength(right);
            edges[index] = 1;
        }

        for (NodeRef tip1 : leftSet) {
            for (NodeRef tip2 : rightSet) {
                int index;
                if (tip1.getNumber() < tip2.getNumber()) {
                    index = (tip1.getNumber() * tree.getExternalNodeCount()) + tip2.getNumber();
                } else {
                    index = (tip2.getNumber() * tree.getExternalNodeCount()) + tip1.getNumber();
                }
                lengths[index] = lengthFromRoot;
                edges[index] = edgesFromRoot;

            }
        }

        Set<NodeRef> tips = new HashSet<NodeRef>();
        tips.addAll(leftSet);
        tips.addAll(rightSet);

        return tips;
    }

    @Deprecated
    public ArrayList<Double> getMetric_old(Tree tree, ArrayList<Double> lambda) {

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
                NodeRef MRCA = TreeUtils.getCommonAncestor(tree, nodeOne, nodeTwo);

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

        int externalNodeCount = tree.getExternalNodeCount();
        int d = (externalNodeCount-2)*(externalNodeCount-1)+externalNodeCount;

        //fill out arrays further
        index = 0;
        for (int i = (externalNodeCount-1)*(externalNodeCount-2); i < d; i++) {
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

    @Deprecated
    public ArrayList<Double> getMetric_old(Tree tree1, Tree tree2, ArrayList<Double> lambda) {

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
                NodeRef MRCA = TreeUtils.getCommonAncestor(tree1, nodeOne, nodeTwo);

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

        int externalNodeCount = tree2.getExternalNodeCount();
        int d = (externalNodeCount-2)*(externalNodeCount-1)+externalNodeCount;

        //fill out arrays further
        index = 0;
        for (int i = (externalNodeCount-1)*(externalNodeCount-2); i < d; i++) {
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
                NodeRef MRCA = TreeUtils.getCommonAncestor(tree2, nodeOne, nodeTwo);

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
        for (int i = (externalNodeCount-1)*(externalNodeCount-2); i < d; i++) {
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
//        tree 1: ((A:1.2,B:0.8):0.5,(C:0.8,D:1.0):1.1);
//        tree 2: (((A:0.8,B:1.3999999999999997):0.30000000000000004,C:0.7000000000000002):0.8999999999999999,D:1.0);
//
//        lambda (0.0) = 2.0
//        lambda (0.5) = 1.9397164741270823
//        lambda (1.0) = 1.962141687034858
//        lambda (0.0) = 2.0
//        lambda (0.5) = 1.9397164741270823
//        lambda (1.0) = 1.962141687034858


        try {

            //4-taxa example
            NewickImporter importer = new NewickImporter("(('A':1.2,'B':0.8):0.5,('C':0.8,'D':1.0):1.1)");
            Tree treeOne = importer.importNextTree();
            System.out.println("4-taxa tree 1: " + treeOne);

            importer = new NewickImporter("((('A':0.8,'B':1.4):0.3,'C':0.7):0.9,'D':1.0)");
            Tree treeTwo = importer.importNextTree();
            System.out.println("4-taxa tree 2: " + treeTwo);
            System.out.println();

            double metrics[] = new double[] {
                    (new KendallColijnPathDifferenceMetric(0.0).getMetric(treeOne, treeTwo)),
                    (new KendallColijnPathDifferenceMetric(0.5).getMetric(treeOne, treeTwo)),
                    (new KendallColijnPathDifferenceMetric(1.0).getMetric(treeOne, treeTwo)),
            };

            System.out.println("Paired trees:");
            System.out.println("lambda (0.0) = " + metrics[0]);
            System.out.println("lambda (0.5) = " + metrics[1]);
            System.out.println("lambda (1.0) = " + metrics[2]);
            System.out.println();

            //Additional test for comparing a collection of trees against a (fixed) focal tree
            metrics = new double[] {
                    (new KendallColijnPathDifferenceMetric(0.0, treeOne).getMetric(treeOne, treeTwo)),
                    (new KendallColijnPathDifferenceMetric(0.5, treeOne).getMetric(treeOne, treeTwo)),
                    (new KendallColijnPathDifferenceMetric(1.0, treeOne).getMetric(treeOne, treeTwo)),
            };

            System.out.println("Focal trees:");
            System.out.println("lambda (0.0) = " + metrics[0]);
            System.out.println("lambda (0.5) = " + metrics[1]);
            System.out.println("lambda (1.0) = " + metrics[2]);
            System.out.println();
            System.out.println();

            //5-taxa example
            importer = new NewickImporter("(((('A':0.6,'B':0.6):0.1,'C':0.5):0.4,'D':0.7):0.1,'E':1.3)");
            treeOne = importer.importNextTree();
            System.out.println("5-taxa tree 1: " + treeOne);

            importer = new NewickImporter("((('A':0.8,'B':1.4):0.1,'C':0.7):0.2,('D':1.0,'E':0.9):1.3)");
            treeTwo = importer.importNextTree();
            System.out.println("5-taxa tree 2: " + treeTwo);
            System.out.println();

            //lambda = 0.0 should yield: sqrt(7) = 2.6457513110645907162
            //lambda = 1.0 should yield: sqrt(2.96) = 1.7204650534085252911

            metrics = new double[] {
                    (new KendallColijnPathDifferenceMetric(0.0, treeOne).getMetric(treeOne, treeTwo)),
                    (new KendallColijnPathDifferenceMetric(0.5, treeOne).getMetric(treeOne, treeTwo)),
                    (new KendallColijnPathDifferenceMetric(1.0, treeOne).getMetric(treeOne, treeTwo)),
            };

            System.out.println("Paired trees:");
            System.out.println("lambda (0.0) = " + metrics[0]);
            System.out.println("lambda (0.5) = " + metrics[1]);
            System.out.println("lambda (1.0) = " + metrics[2]);
            System.out.println();

            //Additional test for comparing a collection of trees against a (fixed) focal tree
            metrics = new double[] {
                    (new KendallColijnPathDifferenceMetric(0.0, treeOne).getMetric(treeOne, treeTwo)),
                    (new KendallColijnPathDifferenceMetric(0.5, treeOne).getMetric(treeOne, treeTwo)),
                    (new KendallColijnPathDifferenceMetric(1.0, treeOne).getMetric(treeOne, treeTwo)),
            };

            System.out.println("Focal trees:");
            System.out.println("lambda (0.0) = " + metrics[0]);
            System.out.println("lambda (0.5) = " + metrics[1]);
            System.out.println("lambda (1.0) = " + metrics[2]);
            System.out.println();


            //timings
//            long startTime = System.currentTimeMillis();
//            for (int i = 0; i < 1000000; i++) {
//                new KCPathDifferenceMetric().getMetric_old(treeOne, treeTwo, lambdaValues);
//            }
//            System.out.println("Old algorithm: " + (System.currentTimeMillis() - startTime) + " ms");

            long startTime = System.currentTimeMillis();
            for (int i = 0; i < 1000000; i++) {
                new KendallColijnPathDifferenceMetric(0.5).getMetric(treeOne, treeTwo);
            }
            System.out.println("New algorithm, 1M reps: " + (System.currentTimeMillis() - startTime) + " ms");

        } catch(Importer.ImportException ie) {
            System.err.println(ie);
        } catch(IOException ioe) {
            System.err.println(ioe);
        }

    }

    @Override
    public Type getType() {
        return TYPE;
    }

    @Override
    public String toString() {
        return getType().getShortName() + "(" + lambda + ")";
    }


}
