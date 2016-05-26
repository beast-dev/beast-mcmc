/*
 * SankoffParsimony.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evolution.parsimony;

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.*;
import dr.evolution.util.Taxon;

import java.util.Set;
import java.util.TreeSet;
import java.util.Iterator;

/**
 * Class for reconstructing characters using the Sankoff generalized parsimony methods. This will be
 * slower than the Fitch algorithm but it allows Weighted Parsimony.
 *
 * @version $Id: SankoffParsimony.java,v 1.7 2005/06/29 16:54:18 beth Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class SankoffParsimony implements ParsimonyCriterion {

    private final int stateCount;

    private int[][] stateSets;
    private double[][][] nodeScores;
    private int[][] nodeStates;

    private Tree tree = null;
    private final PatternList patterns;
    private final double[][] costMatrix;

    private final boolean compressStates = true;

    private boolean hasCalculatedSteps = false;
    private boolean hasRecontructedStates = false;

    private final double[] siteScores;

    public SankoffParsimony(PatternList patterns) {
        if (patterns == null) {
            throw new IllegalArgumentException("The patterns cannot be null");
        }
        stateCount = patterns.getDataType().getStateCount();
        this.costMatrix = new double[stateCount][stateCount];
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                if (i == j) {
                    costMatrix[i][j] = 0.0;
                } else {
                    costMatrix[i][j] = 1.0;
                }
            }
        }

        this.patterns = patterns;
        this.siteScores = new double[patterns.getPatternCount()];
    }

    public SankoffParsimony(PatternList patterns, double[][] costMatrix) {
        if (patterns == null) {
            throw new IllegalArgumentException("The patterns cannot be null");
        }
        stateCount = patterns.getDataType().getStateCount();
        if (costMatrix.length != stateCount || costMatrix[0].length != stateCount) {
            throw new IllegalArgumentException("The cost matrix is of the wrong dimension: expecting " + stateCount + " square");
        }
        this.costMatrix = costMatrix;

        this.patterns = patterns;
        this.siteScores = new double[patterns.getPatternCount()];
    }

    /**
     * Calculates the minimum number of siteScores for the parsimony reconstruction of a
     * a set of character patterns on a tree.
     * @param tree a tree object to reconstruct the characters on
     * @return number of parsimony siteScores
     */
    public double[] getSiteScores(Tree tree) {

        if (tree == null) {
            throw new IllegalArgumentException("The tree cannot be null");
        }

        if (this.tree == null || this.tree != tree) {
            this.tree = tree;

            initialize();
        }

        if (!hasCalculatedSteps) {
            calculateSteps(tree, tree.getRoot(), patterns);
            if (compressStates) {
                for (int i = 0; i < siteScores.length; i++) {
                    double[] Sr = nodeScores[tree.getRoot().getNumber()][i];
                    siteScores[i] = minScore(Sr, stateSets[i]);
                }
            } else {
                for (int i = 0; i < siteScores.length; i++) {
                    double[] Sr = nodeScores[tree.getRoot().getNumber()][i];
                    siteScores[i] = minScore(Sr);
                }
            }
            hasCalculatedSteps = true;
        }


        return siteScores;
    }

    public double getScore(Tree tree) {

        getSiteScores(tree);

        double score = 0;

        for (int i = 0; i < patterns.getPatternCount(); i++) {
            score += siteScores[i] * patterns.getPatternWeight(i);
        }
        return score;
    }

    /**
     * Returns the reconstructed character nodeStates for a given node in the tree. If this method is repeatedly
     * called with the same tree and patterns then only the first call will reconstruct the nodeStates and each
     * subsequent call will return the stored nodeStates.
     * @param tree a tree object to reconstruct the characters on
     * @param node the node of the tree
     * @return an array containing the reconstructed nodeStates for this node
     */
    public int[] getStates(Tree tree, NodeRef node) {

        getSiteScores(tree);

        if (!hasRecontructedStates) {
            if (compressStates) {
                for (int i = 0; i < patterns.getPatternCount(); i++) {
                    nodeStates[tree.getRoot().getNumber()][i] = minState(nodeScores[tree.getRoot().getNumber()][i], stateSets[i]);
                }
            } else {
                for (int i = 0; i < patterns.getPatternCount(); i++) {
                    nodeStates[tree.getRoot().getNumber()][i] = minState(nodeScores[tree.getRoot().getNumber()][i]);
                }
            }
            reconstructStates(tree, tree.getRoot(), nodeStates[tree.getRoot().getNumber()]);
            hasRecontructedStates = true;
        }

        return nodeStates[node.getNumber()];
    }

    private void initialize() {
        hasCalculatedSteps = false;
        hasRecontructedStates = false;

        if (compressStates) {
            stateSets = new int[patterns.getPatternCount()][];
        }

        nodeScores = new double[tree.getNodeCount()][patterns.getPatternCount()][];
        nodeStates = new int[tree.getNodeCount()][patterns.getPatternCount()];

        for (int i = 0; i < patterns.getPatternCount(); i++) {
            int[] pattern = patterns.getPattern(i);

            if (compressStates) {
                Set observedStates = new TreeSet();
                for (int j = 0; j < pattern.length; j++) {
                    boolean[] stateSet = patterns.getDataType().getStateSet(pattern[j]);
                    for (int k = 0; k < stateSet.length; k++) {
                        if (stateSet[k]) {
                            observedStates.add(new Integer(k));
                        }
                    }
                }

                stateSets[i] = new int[observedStates.size()];

                Iterator iter = observedStates.iterator();
                int j = 0;
                while (iter.hasNext()) {
                    stateSets[i][j] = ((Integer)iter.next()).intValue();
                    j++;
                }
            }

            for (int j = 0; j < tree.getExternalNodeCount(); j++) {
                NodeRef node = tree.getExternalNode(j);
                int state = pattern[patterns.getTaxonIndex(tree.getNodeTaxon(node).getId())];
                boolean[] stateSet = patterns.getDataType().getStateSet(state);

                nodeScores[j][i] = new double[stateCount];
                for (int k = 0; k < stateCount; k++) {
                    if (stateSet[k]) {
                        nodeScores[j][i][k] = 0.0;
                    } else {
                        nodeScores[j][i][k] = Double.POSITIVE_INFINITY;
                    }
                }
            }

            for (int j = 0; j < tree.getInternalNodeCount(); j++) {
                nodeScores[j + tree.getExternalNodeCount()][i] = new double[stateCount];
            }

        }
    }

    /**
     * This is the first pass of the Fitch algorithm. This calculates the set of nodeStates
     * at each node and counts the total number of siteScores (the score). If that is all that
     * is required then the second pass is not necessary.
     * @param tree
     * @param node
     * @param patterns
     */
    private void calculateSteps(Tree tree, NodeRef node, PatternList patterns) {

        if (!tree.isExternal(node)) {

            for (int i = 0; i < tree.getChildCount(node); i++) {
                calculateSteps(tree, tree.getChild(node, i), patterns);
            }

            if (compressStates) {
                for (int i = 0; i < patterns.getPatternCount(); i++) {
                    double[] Sc = nodeScores[tree.getChild(node, 0).getNumber()][i];
                    double[] Sa = nodeScores[node.getNumber()][i];

                    int[] set = stateSets[i];
                    for (int k = 0; k < set.length; k++) {
                        Sa[set[k]] = minCost(k, Sc, costMatrix, set);
                    }

                    for (int j = 1; j < tree.getChildCount(node); j++) {
                        Sc = nodeScores[tree.getChild(node, j).getNumber()][i];
                        for (int k = 0; k < set.length; k++) {
                            Sa[set[k]] += minCost(k, Sc, costMatrix, set);
                        }

                    }


                }
            } else {
                for (int i = 0; i < patterns.getPatternCount(); i++) {
                    double[] Sc = nodeScores[tree.getChild(node, 0).getNumber()][i];
                    double[] Sa = nodeScores[node.getNumber()][i];

                    for (int k = 0; k < stateCount; k++) {
                        Sa[k] = minCost(k, Sc, costMatrix);
                    }

                    for (int j = 1; j < tree.getChildCount(node); j++) {
                        Sc = nodeScores[tree.getChild(node, j).getNumber()][i];
                        for (int k = 0; k < stateCount; k++) {
                            Sa[k] += minCost(k, Sc, costMatrix);
                        }

                    }


                }

            }
        }
    }

    /**
     * The second pass of the algorithm. This reconstructs the ancestral nodeStates at
     * each node.
     * @param tree
     * @param node
     * @param parentStates
     */
    private void reconstructStates(Tree tree, NodeRef node, int[] parentStates) {

        for (int i = 0; i < patterns.getPatternCount(); i++) {

            double[] Sa = nodeScores[node.getNumber()][i];

            if (compressStates) {
                int[] set = stateSets[i];

                int minState = set[0];
                double minCost = Sa[minState] + costMatrix[parentStates[i]][minState];

                for (int j = 1; j < set.length; j++) {
                    double c = Sa[set[j]] + costMatrix[parentStates[i]][set[j]];
                    if (c < minCost) {
                        minState = set[j];
                        minCost = c;
                    }
                }
                nodeStates[node.getNumber()][i] = minState;
            } else {
                int minState = 0;
                double minCost = Sa[minState] + costMatrix[parentStates[i]][minState];

                for (int j = 1; j < Sa.length; j++) {
                    double c = Sa[j] + costMatrix[parentStates[i]][j];
                    if (c < minCost) {
                        minState = j;
                        minCost = c;
                    }
                }
                nodeStates[node.getNumber()][i] = minState;
            }

        }

        for (int i = 0; i < tree.getChildCount(node); i++) {
            reconstructStates(tree, tree.getChild(node, i), nodeStates[node.getNumber()]);
        }

    }


    private int minState(double[] s1) {

        int minState = 0;

        for (int j = 1; j < s1.length; j++) {
            if (s1[j] < s1[minState]) minState = j;
        }
        return minState;
    }

    private double minScore(double[] s1) {

        double minScore = s1[0];

        for (int j = 1; j < s1.length; j++) {
            if (s1[j] < minScore) minScore = s1[j];
        }
        return minScore;
    }

    private double minCost(int i, double[] s1, double[][] costMatrix) {

        double[] costRow = costMatrix[i];
        double minCost = costRow[0] + s1[0];

        for (int j = 1; j < s1.length; j++) {
            double cost = costRow[j] + s1[j];
            if (cost < minCost) minCost = cost;
        }
        return minCost;
    }

    private int minState(double[] s1, int[] set) {

        int minState = set[0];

        for (int j = 1; j < set.length; j++) {
            if (s1[set[j]] < s1[minState]) minState = set[j];
        }
        return minState;
    }

    private double minScore(double[] s1, int[] set) {

        double minScore = s1[set[0]];

        for (int j = 1; j < set.length; j++) {
            if (s1[set[j]] < minScore) minScore = s1[set[j]];
        }
        return minScore;
    }

    private double minCost(int i, double[] s1, double[][] costMatrix, int[] set) {

        double[] costRow = costMatrix[set[i]];
        double minCost = costRow[set[0]] + s1[set[0]];

        for (int j = 1; j < set.length; j++) {
            double cost = costRow[set[j]] + s1[set[j]];
            if (cost < minCost) minCost = cost;
        }
        return minCost;
    }

    public static void main(String[] argv) {
        FlexibleNode tip1 = new FlexibleNode(new Taxon("tip1"));
        FlexibleNode tip2 = new FlexibleNode(new Taxon("tip2"));
        FlexibleNode tip3 = new FlexibleNode(new Taxon("tip3"));
        FlexibleNode tip4 = new FlexibleNode(new Taxon("tip4"));
        FlexibleNode tip5 = new FlexibleNode(new Taxon("tip5"));

        FlexibleNode node1 = new FlexibleNode();
        node1.addChild(tip1);
        node1.addChild(tip2);

        FlexibleNode node2 = new FlexibleNode();
        node2.addChild(tip4);
        node2.addChild(tip5);

        FlexibleNode node3 = new FlexibleNode();
        node3.addChild(tip3);
        node3.addChild(node2);

        FlexibleNode root = new FlexibleNode();
        root.addChild(node1);
        root.addChild(node3);

        FlexibleTree tree = new FlexibleTree(root);

        Patterns patterns = new Patterns(Nucleotides.INSTANCE, tree);
        //patterns.addPattern(new int[] {1, 0, 1, 2, 2});
        //patterns.addPattern(new int[] {2, 1, 1, 1, 2});
        patterns.addPattern(new int[] {2, 3, 1, 3, 3});

        FitchParsimony fitch = new FitchParsimony(patterns, false);
        SankoffParsimony sankoff = new SankoffParsimony(patterns);

        for (int i = 0; i < patterns.getPatternCount(); i++) {
            double[] scores = fitch.getSiteScores(tree);
            System.out.println("Pattern = " + i);

            System.out.println("Fitch:");
            System.out.println("  No. Steps = " + scores[i]);
            System.out.println("    state(node1) = " + fitch.getStates(tree, node1)[i]);
            System.out.println("    state(node2) = " + fitch.getStates(tree, node2)[i]);
            System.out.println("    state(node3) = " + fitch.getStates(tree, node3)[i]);
            System.out.println("    state(root) = " + fitch.getStates(tree, root)[i]);

            scores = sankoff.getSiteScores(tree);
            System.out.println("Sankoff:");
            System.out.println("  No. Steps = " + scores[i]);
            System.out.println("    state(node1) = " + sankoff.getStates(tree, node1)[i]);
            System.out.println("    state(node2) = " + sankoff.getStates(tree, node2)[i]);
            System.out.println("    state(node3) = " + sankoff.getStates(tree, node3)[i]);
            System.out.println("    state(root) = " + sankoff.getStates(tree, root)[i]);

            System.out.println();
        }
    }

}