/*
 * ImportanceSubtreeSwap.java
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
package dr.evomodel.operators;

import dr.evolution.tree.MutableTree.InvalidTreeException;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.ConditionalCladeFrequency;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.ImportanceSubtreeSwapParser;
import dr.inference.operators.*;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Sebastian Hoehna
 *         <p/>
 *         This class implements a subtree swap operator. The first subtree is
 *         chosen randomly and the second one is chosen according to the
 *         importance of the new tree. The importance are calculated by the
 *         multiplied clade probabilities.
 */
public class ImportanceSubtreeSwap extends AbstractTreeOperator {

    public final int SAMPLE_EVERY = 10;

    private final TreeModel tree;

    private final int samples;

    private int sampleCount = 0;

    private boolean burnin = false;

    private final ConditionalCladeFrequency probabilityEstimater;

    private final OperatorSchedule schedule;

    /**
     *
     */
    public ImportanceSubtreeSwap(TreeModel tree, double weight, int samples, int epsilon) {
        this.tree = tree;
        setWeight(weight);
        this.samples = samples;
        sampleCount = 0;
        probabilityEstimater = new ConditionalCladeFrequency(tree, epsilon);
        schedule = getOperatorSchedule(tree);
    }

    /**
     *
     */
    public ImportanceSubtreeSwap(TreeModel tree, double weight, int samples) {
        this.tree = tree;
        setWeight(weight);
        this.samples = samples;
        sampleCount = 0;
        double epsilon = 1 - Math.pow(0.5, 1.0 / samples);
        probabilityEstimater = new ConditionalCladeFrequency(tree, epsilon);
        schedule = getOperatorSchedule(tree);
    }

    private OperatorSchedule getOperatorSchedule(TreeModel treeModel) {

        ExchangeOperator narrowExchange = new ExchangeOperator(
                ExchangeOperator.NARROW, treeModel, 10);
        ExchangeOperator wideExchange = new ExchangeOperator(
                ExchangeOperator.WIDE, treeModel, 3);
        SubtreeSlideOperator subtreeSlide = new SubtreeSlideOperator(treeModel,
                10.0, 1.0, true, false, false, false, CoercionMode.COERCION_ON);
        NNI nni = new NNI(treeModel, 10.0);
        WilsonBalding wilsonBalding = new WilsonBalding(treeModel, 3.0);
        FNPR fnpr = new FNPR(treeModel, 5.0);

        OperatorSchedule schedule = new SimpleOperatorSchedule();
        schedule.addOperator(narrowExchange);
        schedule.addOperator(wideExchange);
        schedule.addOperator(subtreeSlide);
        schedule.addOperator(nni);
        schedule.addOperator(wilsonBalding);
        schedule.addOperator(fnpr);

        return schedule;
    }

    /*
     * (non-Javadoc)
     *
     * @see dr.inference.operators.SimpleMCMCOperator#doOperation()
     */
    @Override
    public double doOperation() throws OperatorFailedException {
        if (!burnin) {
            if (sampleCount < samples * SAMPLE_EVERY) {
                sampleCount++;
                if (sampleCount % SAMPLE_EVERY == 0) {
                    probabilityEstimater.addTree(tree);
                }
                setAcceptCount(0);
                setRejectCount(0);
                setTransitions(0);

                return doUnguidedOperation();

            } else {
                return importanceExchange();
            }
        } else {

            return doUnguidedOperation();

        }
    }

    private double doUnguidedOperation() throws OperatorFailedException {
        int index = schedule.getNextOperatorIndex();
        SimpleMCMCOperator operator = (SimpleMCMCOperator) schedule.getOperator(index);

        return operator.doOperation();
    }

    /**
     * WARNING: Assumes strictly bifurcating tree.
     *
     * @throws InvalidTreeException
     */
    private double importanceExchange() throws OperatorFailedException {

        final int nodeCount = tree.getNodeCount();
        final NodeRef root = tree.getRoot();

        NodeRef i;
        int indexI;
        int indexJ;

        do {
            indexI = MathUtils.nextInt(nodeCount);
            i = tree.getNode(indexI);
        } while (root == i
                || (tree.getParent(i) == root &&
                tree.getNodeHeight(i) > tree.getNodeHeight(getOtherChild(tree, tree.getParent(i), i))));

        List<Integer> secondNodeIndices = new ArrayList<Integer>();
        List<Double> probabilities = new ArrayList<Double>();
        NodeRef j, iP, jP;
        iP = tree.getParent(i);
        double sum = 0.0;
        double backward = calculateTreeProbability(tree);
        int offset = (int) -backward;
        backward = Math.exp(backward + offset);

        tree.beginTreeEdit();
        for (int n = 0; n < nodeCount; n++) {
            j = tree.getNode(n);
            if (j != root) {
                jP = tree.getParent(j);

                if ((iP != jP) && (i != jP) && (j != iP)
                        && (tree.getNodeHeight(j) < tree.getNodeHeight(iP))
                        && (tree.getNodeHeight(i) < tree.getNodeHeight(jP))) {
                    secondNodeIndices.add(n);

                    swap(tree, tree.getNode(indexI), tree.getNode(n));
                    double prob = Math.exp(calculateTreeProbability(tree)
                            + offset);
                    probabilities.add(prob);
                    swap(tree, tree.getNode(indexI), tree.getNode(n));
                    sum += prob;
                }
            }
        }

        double ran = Math.random() * sum;
        int index = 0;
        while (ran > 0.0) {
            ran -= probabilities.get(index);
            index++;
        }
        index--;

        j = tree.getNode(secondNodeIndices.get(index));
        jP = tree.getParent(j);

        // *******************************************
        // assuming we would have chosen j first
        double sumForward2 = 0.0;
        NodeRef k, kP;
        indexJ = secondNodeIndices.get(index);
        for (int n = 0; n < nodeCount; n++) {
            k = tree.getNode(n);
            if (k != root) {
                kP = tree.getParent(k);

                if ((jP != kP) && (j != kP) && (k != jP)
                        && (tree.getNodeHeight(k) < tree.getNodeHeight(jP))
                        && (tree.getNodeHeight(j) < tree.getNodeHeight(kP))) {

                    swap(tree, tree.getNode(indexJ), tree.getNode(n));
                    double prob = Math.exp(calculateTreeProbability(tree)
                            + offset);
                    sumForward2 += prob;
                    swap(tree, tree.getNode(indexJ), tree.getNode(n));
                }
            }
        }

        swap(tree, i, j);
        double forward = probabilities.get(index);

        iP = tree.getParent(i);
        double sumBackward = 0.0;
        for (int n = 0; n < nodeCount; n++) {
            j = tree.getNode(n);
            if (j != root) {
                jP = tree.getParent(j);

                if ((iP != jP) && (i != jP) && (j != iP)
                        && (tree.getNodeHeight(j) < tree.getNodeHeight(iP))
                        && (tree.getNodeHeight(i) < tree.getNodeHeight(jP))) {

                    swap(tree, tree.getNode(indexI), tree.getNode(n));
                    double prob = Math.exp(calculateTreeProbability(tree)
                            + offset);
                    sumBackward += prob;
                    swap(tree, tree.getNode(indexI), tree.getNode(n));

                }
            }
        }

        // *******************************************
        // assuming we would have chosen j first
        double sumBackward2 = 0.0;
        j = tree.getNode(secondNodeIndices.get(index));
        jP = tree.getParent(j);
        for (int n = 0; n < nodeCount; n++) {
            k = tree.getNode(n);
            if (k != root) {
                kP = tree.getParent(k);

                if ((jP != kP) && (j != kP) && (k != jP)
                        && (tree.getNodeHeight(k) < tree.getNodeHeight(jP))
                        && (tree.getNodeHeight(j) < tree.getNodeHeight(kP))) {

                    swap(tree, tree.getNode(indexJ), tree.getNode(n));
                    double prob = Math.exp(calculateTreeProbability(tree)
                            + offset);
                    sumBackward2 += prob;
                    swap(tree, tree.getNode(indexJ), tree.getNode(n));
                }
            }
        }

        tree.endTreeEdit();

        // AR - not sure whether this check is necessary
        try {
            tree.checkTreeIsValid();
        } catch (InvalidTreeException e) {
            throw new OperatorFailedException(e.getMessage());
        }


        double forwardProb = (forward / sum) + (forward / sumForward2);
        double backwardProb = (backward / sumBackward)
                + (backward / sumBackward2);

        double hastingsRatio = Math.log(backwardProb / forwardProb);

        // throw new OperatorFailedException(
        // "Couldn't find valid wide move on this tree!");

        return hastingsRatio;
    }

    /* exchange subtrees whose root are i and j */
    private void swap(TreeModel tree, NodeRef i, NodeRef j)
            throws OperatorFailedException {

        NodeRef iP = tree.getParent(i);
        NodeRef jP = tree.getParent(j);

        tree.removeChild(iP, i);
        tree.removeChild(jP, j);
        tree.addChild(jP, i);
        tree.addChild(iP, j);
    }

    private double calculateTreeProbability(Tree tree) {
        // return calculateTreeProbabilityMult(tree);
        // return calculateTreeProbabilityLog(tree);
        return probabilityEstimater.getTreeProbability(tree);
        // return 0.0;
    }

    public void setBurnin(boolean burnin) {
        this.burnin = burnin;
    }

    /*
     * (non-Javadoc)
     *
     * @see dr.inference.operators.SimpleMCMCOperator#getOperatorName()
     */
    @Override
    public String getOperatorName() {
        return ImportanceSubtreeSwapParser.IMPORTANCE_SUBTREE_SWAP;
    }

    /*
     * (non-Javadoc)
     *
     * @see dr.inference.operators.MCMCOperator#getPerformanceSuggestion()
     */
    public String getPerformanceSuggestion() {
        // TODO Auto-generated method stub
        return "";
    }

}
