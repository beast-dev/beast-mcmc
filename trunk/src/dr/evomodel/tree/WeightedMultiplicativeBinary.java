/*
 * WeightedMultiplicativeBinary.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

/**
 *
 */
package dr.evomodel.tree;

import dr.evolution.io.Importer;
import dr.evolution.io.TreeTrace;
import dr.evolution.tree.Clade;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.inference.model.Likelihood;
import dr.inference.prior.Prior;
import dr.math.MathUtils;

import java.io.IOException;
import java.io.Reader;
import java.util.*;

/**
 * @author Sebastian Hoehna
 */
public class WeightedMultiplicativeBinary extends AbstractCladeImportanceDistribution {

    private final int TAXA_COUNT;

    private double EPSILON;

    private long samples = 0;

    private HashMap<BitSet, Clade> cladeProbabilities;

    private TreeTrace[] traces;

    private int burnin;

    /**
     * @param epsilon - the default number of occurences for each clade which wasn't
     *                observed to guarantee non-zero probabilities
     */
    public WeightedMultiplicativeBinary(Tree tree, double epsilon) {

        // initializing global variables
        cladeProbabilities = new HashMap<BitSet, Clade>();

        // setting global variables
        EPSILON = epsilon;

        TAXA_COUNT = tree.getExternalNodeCount();

    }

    /**
     * @param traces  - samples of trees in a tree traces array.
     * @param epsilon - the default number of occurences for each clade which wasn't
     *                observed to guarantee non-zero probabilities
     * @param burnIn  - number of trees discarded from the trace
     * @param verbose - hide the runtime status and outputs
     */
    public WeightedMultiplicativeBinary(TreeTrace[] traces, double epsilon,
                                        int burnIn, boolean verbose) {

        // initializing global variables
        cladeProbabilities = new HashMap<BitSet, Clade>();

        // setting global variables
        EPSILON = epsilon;
        this.traces = traces;

        // calculates the burn-in to 10% if it was set out of the boundaries
        int minMaxState = Integer.MAX_VALUE;
        for (TreeTrace trace : traces) {
            if (trace.getMaximumState() < minMaxState) {
                minMaxState = trace.getMaximumState();
            }
        }


        Tree tree = traces[0].getTree(0, burnIn);
        TAXA_COUNT = tree.getExternalNodeCount();

        if (burnIn < 0 || burnIn >= minMaxState) {
            this.burnin = minMaxState / (10 * traces[0].getStepSize());
            if (verbose)
                System.out
                        .println("WARNING: Burn-in larger than total number of states - using 10% of smallest trace");
        } else {
            this.burnin = burnIn;
        }

        // analyzing the whole trace -> reading the trees
        analyzeTrace(verbose);
    }

    /**
     * Actually analyzes the trace given the burn-in. Each tree from the trace
     * is read and the conditional clade frequencies incremented.
     *
     * @param verbose if true then progress is logged to stdout
     */
    public void analyzeTrace(boolean verbose) {

        if (verbose) {
            if (traces.length > 1)
                System.out.println("Combining " + traces.length + " traces.");
        }

        // get first tree to extract the taxon
        Tree tree = getTree(0);

        // read every tree from the trace
        for (TreeTrace trace : traces) {
            // do some output stuff
            int treeCount = trace.getTreeCount(burnin * trace.getStepSize());
            double stepSize = treeCount / 60.0;
            int counter = 1;

            if (verbose) {
                System.out.println("Analyzing " + treeCount + " trees...");
                System.out
                        .println("0              25             50             75            100");
                System.out
                        .println("|--------------|--------------|--------------|--------------|");
                System.out.print("*");
            }
            for (int i = 1; i < treeCount; i++) {
                // get the next tree
                tree = trace.getTree(i, burnin * trace.getStepSize());

                // add the tree and its clades to the frequencies
                addTree(tree);

                // some more output stuff
                if (i >= (int) Math.round(counter * stepSize) && counter <= 60) {
                    if (verbose) {
                        System.out.print("*");
                        System.out.flush();
                    }
                    counter += 1;
                }
            }
            if (verbose) {
                System.out.println("*");
            }
        }
    }

    /**
     * Creates the report. The estimated posterior of the given tree is printed.
     *
     * @throws IOException if general I/O error occurs
     */
    public void report(Tree tree) throws IOException {

        System.err.println("making report");

        SimpleTree sTree = new SimpleTree(tree);
        System.out
                .println("Estimated marginal posterior by condiational clade frequencies:");
        System.out.println(getTreeProbability(sTree));

        System.out.flush();

    }

    /**
     * Calculates the probability of a given tree.
     *
     * @param tree - the tree to be analyzed
     * @return estimated posterior probability in log
     */
    public double getTreeProbability(SimpleTree tree) {
        return calculateTreeProbabilityLog(tree);
//		return calculateTreeProbabilityLogRecursive(tree, tree.getRoot());
    }

    /**
     * Calculates the probability of a given tree.
     *
     * @param tree - the tree to be analyzed
     * @return estimated posterior probability in log
     */
    public double getTreeProbability(SimpleTree tree, HashMap<String, Integer> taxonMap) {
        return calculateTreeProbabilityLog(tree, taxonMap);
//		return calculateTreeProbabilityLogRecursive(tree, tree.getRoot());
    }

    /**
     * Calculates the probability of a given tree.
     *
     * @param tree - the tree to be analyzed
     * @return estimated posterior probability in log
     */
    private double calculateTreeProbabilityLog(Tree tree) {
        double prob = 0.0;

        // calculate the number of possible splits
        final double splits = Math.pow(2, tree.getExternalNodeCount() - 1) - 1;

        List<Clade> clades = new ArrayList<Clade>();
        List<Clade> parentClades = new ArrayList<Clade>();
        // get clades contained in the tree
        getClades(tree, tree.getRoot(), parentClades, clades);

        // for every clade multiply its probability to the
        // tree probability
        for (Clade c : clades) {

            // set the occurrences to epsilon
            double occurrences = EPSILON;
            if (cladeProbabilities.containsKey(c.getBits())) {
                // if we observed this clade in the trace, add the occurrences
                // to epsilon
                occurrences += cladeProbabilities.get(c.getBits())
                        .getSampleCount();
            }

            // multiply the conditional clade probability to the tree
            // probability
            prob += Math.log(occurrences / (samples + (splits * EPSILON)));

        }

        return prob;
    }

    /**
     * Calculates the probability of a given tree.
     *
     * @param tree - the tree to be analyzed
     * @return estimated posterior probability in log
     */
    private double calculateTreeProbabilityLog(Tree tree, HashMap<String, Integer> taxonMap) {
        double prob = 0.0;

        // calculate the number of possible splits
        final double splits = Math.pow(2, tree.getExternalNodeCount() - 1) - 1;

        List<Clade> clades = new ArrayList<Clade>();
        List<Clade> parentClades = new ArrayList<Clade>();
        // get clades contained in the tree
        getClades(tree, tree.getRoot(), parentClades, clades, taxonMap);

        // for every clade multiply its probability to the
        // tree probability
        for (Clade c : clades) {

            // set the occurrences to epsilon
            double occurrences = EPSILON;
            if (cladeProbabilities.containsKey(c.getBits())) {
                // if we observed this clade in the trace, add the occurrences
                // to epsilon
                double cladesInTreeSpace = getTrees(c.getBits().cardinality()) * getTrees(TAXA_COUNT - c.getBits().cardinality() + 1);
                occurrences += (cladeProbabilities.get(c.getBits())
                        .getSampleCount() / cladesInTreeSpace);
            }

            // multiply the conditional clade probability to the tree
            // probability
            prob += Math.log(occurrences / (samples + (splits * EPSILON)));

        }

        return prob;
    }

    private double getTrees(int n) {
        double trees = 1;
        for (int i = 3; i <= n; i++) {
            trees *= (2 * i - 3);
        }

        return trees;
    }

    /**
     * Calculates the probability of a given tree recursively.
     *
     * @param tree - the tree to be analyzed
     * @param node - the node at which the subtree is rooted for which the
     *             probability has to be calculated
     * @return estimated posterior probability in log
     */
    private double calculateTreeProbabilityLogRecursive(Tree tree, NodeRef node) {
        double prob = 0.0;

        NodeRef leftChild = tree.getChild(node, 0);
        NodeRef rightChild = tree.getChild(node, 1);

        if (tree.isExternal(leftChild) && tree.isExternal(rightChild)) {
            // both children are external nodes
            return 0.0;
        } else if (!tree.isExternal(leftChild) && !tree.isExternal(rightChild)) {
            // both children are internal nodes
            Clade leftSubclade = getClade(tree, leftChild);
            Clade rightSubclade = getClade(tree, rightChild);

            double sum = 0.0;
            if (cladeProbabilities.containsKey(leftSubclade.getBits())) {
                sum += (cladeProbabilities.get(leftSubclade.getBits())
                        .getSampleCount() + EPSILON)
                        / samples;
            } else {
                sum += EPSILON / samples;
            }
            if (cladeProbabilities.containsKey(rightSubclade.getBits())) {
                sum += (cladeProbabilities.get(rightSubclade.getBits())
                        .getSampleCount() + EPSILON)
                        / samples;
            } else {
                sum += EPSILON / samples;
            }

            prob += Math.log(sum / 2.0);

            prob += calculateTreeProbabilityLogRecursive(tree, leftChild);
            prob += calculateTreeProbabilityLogRecursive(tree, rightChild);

            return prob;
        } else {

            Clade leftSubclade = getClade(tree, leftChild);
            Clade rightSubclade = getClade(tree, rightChild);

            double sum = 0.0;
            if (leftSubclade.getSize() > 1) {
                if (cladeProbabilities.containsKey(leftSubclade.getBits())) {
                    sum += (cladeProbabilities.get(leftSubclade.getBits())
                            .getSampleCount() + EPSILON)
                            / samples;
                } else {
                    sum += EPSILON / samples;
                }
            }
            if (rightSubclade.getSize() > 1) {
                if (cladeProbabilities.containsKey(rightSubclade.getBits())) {
                    sum += (cladeProbabilities.get(rightSubclade.getBits())
                            .getSampleCount() + EPSILON)
                            / samples;
                } else {
                    sum += EPSILON / samples;
                }
            }

            prob += Math.log(sum);

            if (!tree.isExternal(leftChild)) {
                prob += calculateTreeProbabilityLogRecursive(tree, leftChild);
            }
            if (!tree.isExternal(rightChild)) {
                prob += calculateTreeProbabilityLogRecursive(tree, rightChild);
            }

            return prob;
        }
    }

    /* (non-Javadoc)
      * @see dr.evomodel.tree.AbstractCladeImportanceDistribution#getChanceForNodeHeights(dr.evomodel.tree.TreeModel, dr.inference.model.Likelihood, dr.inference.prior.Prior)
      */
    @Override
    public double getChanceForNodeHeights(TreeModel tree,
                                          Likelihood likelihood, Prior prior) {
        // TODO Auto-generated method stub
        return 0;
    }

    /* (non-Javadoc)
      * @see dr.evomodel.tree.AbstractCladeImportanceDistribution#setNodeHeights(dr.evomodel.tree.TreeModel, dr.inference.model.Likelihood, dr.inference.prior.Prior)
      */
    @Override
    public double setNodeHeights(TreeModel tree, Likelihood likelihood,
                                 Prior prior) {
        // TODO Auto-generated method stub
        return 0;
    }

    /**
     * get the i'th tree of the trace
     *
     * @param index
     * @return the i'th tree of the trace
     */
    public final Tree getTree(int index) {

        int oldTreeCount = 0;
        int newTreeCount = 0;
        for (TreeTrace trace : traces) {
            newTreeCount += trace.getTreeCount(burnin * trace.getStepSize());

            if (index < newTreeCount) {
                return trace.getTree(index - oldTreeCount, burnin
                        * trace.getStepSize());
            }
            oldTreeCount = newTreeCount;
        }
        throw new RuntimeException("Couldn't find tree " + index);
    }

    /**
     * increments the number of occurrences for all conditional clades
     *
     * @param tree - the tree to be added
     */
    public void addTree(Tree tree) {

        samples++;


        List<Clade> clades = new ArrayList<Clade>();
        List<Clade> parentClades = new ArrayList<Clade>();
        // get clades contained in the tree
        getClades(tree, tree.getRoot(), parentClades, clades);

        // increment the occurrences of the clade and the conditional clade
        for (Clade c : clades) {

            // increment the clade occurrences
            if (cladeProbabilities.containsKey(c.getBits())) {
                Clade tmp = cladeProbabilities.get(c.getBits());
                tmp.addHeight(c.getHeight());
                // frequency += cladeProbabilities.get(c);
            } else {
                // just to set the first value of the height value list
                c.addHeight(c.getHeight());
                cladeProbabilities.put(c.getBits(), c);
            }
        }
    }

    /**
     * increments the number of occurrences for all conditional clades
     *
     * @param tree - the tree to be added
     */
    public void addTree(Tree tree, HashMap<String, Integer> taxonMap) {

        samples++;


        List<Clade> clades = new ArrayList<Clade>();
        List<Clade> parentClades = new ArrayList<Clade>();
        // get clades contained in the tree
        getClades(tree, tree.getRoot(), parentClades, clades, taxonMap);

        // increment the occurrences of the clade and the conditional clade
        for (Clade c : clades) {

            // increment the clade occurrences
            if (cladeProbabilities.containsKey(c.getBits())) {
                Clade tmp = cladeProbabilities.get(c.getBits());
                tmp.addHeight(c.getHeight());
                // frequency += cladeProbabilities.get(c);
            } else {
                // just to set the first value of the height value list
                c.addHeight(c.getHeight());
                cladeProbabilities.put(c.getBits(), c);
            }
        }
    }

    /**
     * @param reader  the readers to be analyzed
     * @param burnin  the burnin in states
     * @param verbose true if progress should be logged to stdout
     * @return an analyses of the trees in a log file.
     * @throws java.io.IOException if general I/O error occurs
     */
    public static ConditionalCladeFrequency analyzeLogFile(Reader[] reader,
                                                           double e, int burnin, boolean verbose) throws IOException {

        TreeTrace[] trace = new TreeTrace[reader.length];
        for (int i = 0; i < reader.length; i++) {
            try {
                trace[i] = TreeTrace.loadTreeTrace(reader[i]);
            } catch (Importer.ImportException ie) {
                throw new RuntimeException(ie.toString());
            }
            reader[i].close();

        }

        return new ConditionalCladeFrequency(trace, e, burnin, verbose);
    }

    /*
      * (non-Javadoc)
      *
      * @see
      * dr.evolution.tree.ImportanceDistribution#getTreeProbability(dr.evolution
      * .tree.Tree)
      */
    public double getTreeProbability(Tree tree) {
        return calculateTreeProbabilityLogRecursive(tree, tree.getRoot());
    }

    /*
      * (non-Javadoc)
      *
      * @see
      * dr.evolution.tree.ImportanceDistribution#splitClade(dr.evolution.tree
      * .Clade, dr.evolution.tree.Clade[])
      */
    public double splitClade(Clade parent, Clade[] children) {
        // the number of all possible clades is 2^n with n the number of tips
        // reduced by 2 because we wont consider the clades with all or no tips
        // contained
        // note: this time we consider each clade of a split separately with its
        // own probability because every clade has a different chance for
        // itself.
        // #splits = 2^(n) - 1
        final double splits = Math.pow(2, parent.getSize()) - 1;

        double prob = 0;
        double sum = 0.0;

        List<Clade> childClades = getPossibleChildren(parent);

        for (Clade child : childClades) {
            sum += child.getSampleCount();
        }
        sum += EPSILON * splits;

        double randomNumber = Math.random() * sum;
        for (Clade child : childClades) {
            randomNumber -= (child.getSampleCount() + EPSILON);
            if (randomNumber < 0) {
                children[0] = child;
                double chance = (child.getSampleCount() + EPSILON) / samples;
                // the other clade which would have resulted into the same split
                BitSet secondChild = (BitSet) children[0].getBits().clone();
                secondChild.xor(parent.getBits());

                if (secondChild.cardinality() > 1) {
                    Clade counterClade = cladeProbabilities.get(secondChild);
                    if (counterClade != null) {
                        chance += (counterClade.getSampleCount() + EPSILON)
                                / samples;
                    } else {
                        chance += EPSILON / samples;
                    }
                    prob = chance / 2.0;
                } else {
                    prob = chance;
                }

                break;
            }
        }

        // we take a clade which we haven't seen so far
        if (randomNumber >= 0) {
//			System.out.println("Random Clade");
            BitSet newChild;
            do {
                do {
                    newChild = (BitSet) parent.getBits().clone();
                    int index = -1;
                    do {
                        index = newChild.nextSetBit(index + 1);
                        if (index > -1 && MathUtils.nextBoolean()) {
                            newChild.clear(index);
                        }
                    } while (index > -1);
                } while (newChild.cardinality() == 0
                        || newChild.cardinality() == parent.getSize());
            } while (cladeProbabilities.containsKey(newChild));

            Clade randomClade = new Clade(newChild, 0.5);
            children[0] = randomClade;

            BitSet secondChild = (BitSet) children[0].getBits().clone();
            secondChild.xor(parent.getBits());

            if (cladeProbabilities.containsKey(secondChild)) {
                children[1] = cladeProbabilities.get(secondChild);
            } else {
                children[1] = new Clade(secondChild, 0.5);
            }

            if (children[0].getSize() > 1 && children[1].getSize() > 1) {
                prob = (children[0].getSampleCount()
                        + children[1].getSampleCount() + (2.0 * EPSILON))
                        / (samples * 2.0);
            } else {
                if (children[0].getSize() > 1) {
                    prob = (children[0].getSampleCount() + EPSILON) / samples;
                } else {
                    prob = (children[1].getSampleCount() + EPSILON) / samples;
                }
            }
        } else {
            BitSet secondChild = (BitSet) children[0].getBits().clone();
            secondChild.xor(parent.getBits());
            children[1] = cladeProbabilities.get(secondChild);
            // children[1] = childClades.get(secondChild);
            if (children[1] == null) {
                children[1] = new Clade(secondChild, 0.5);
                children[1].addHeight(0.5);
            }
        }

        return Math.log(prob);

    }

    /**
     * Finds all possible children clades which we have observed already. A
     * clade is a possible child clade if is a subset of taxa of the parent
     *
     * @param parent - the parent clade of which we want to find the possible
     *               children
     * @return a List<Clade> of the possible child clades
     */
    private List<Clade> getPossibleChildren(Clade parent) {
        List<Clade> children = new ArrayList<Clade>();

        Set<BitSet> keys = cladeProbabilities.keySet();
        for (BitSet key : keys) {
            if (key.cardinality() < parent.getSize()) {
                if (containsClade(parent.getBits(), key)) {
                    children.add(cladeProbabilities.get(key));
                }
            }
        }

        return children;
    }
}
