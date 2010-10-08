/*
 * ImportancePruneAndRegraft.java
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
import dr.evomodelxml.operators.ImportancePruneAndRegraftParser;
import dr.inference.operators.*;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Sebastian Hoehna
 */
public class ImportancePruneAndRegraft extends AbstractTreeOperator {

    public final int SAMPLE_EVERY = 10;

    private final TreeModel tree;

    private final int samples;

    private int sampleCount;

    private boolean burnin = false;

    private final ConditionalCladeFrequency probabilityEstimater;

    private final OperatorSchedule schedule;

    /**
     *
     */
    public ImportancePruneAndRegraft(TreeModel tree, double weight,
                                     int samples, int epsilon) {
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
    public ImportancePruneAndRegraft(TreeModel tree, double weight, int samples) {
        this.tree = tree;
        setWeight(weight);
        this.samples = samples;
        sampleCount = 0;
        // double epsilon = 1 - Math.pow(0.5, samples);
        double epsilon = 1 - Math.pow(0.5, 1.0 / samples);
        // double epsilon = 1;
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
                return importancePruneAndRegraft();
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

    private double importancePruneAndRegraft() throws OperatorFailedException {
        final int nodeCount = tree.getNodeCount();
        final NodeRef root = tree.getRoot();

        NodeRef i;

        do {
            int indexI = MathUtils.nextInt(nodeCount);
            i = tree.getNode(indexI);
        } while (root == i || tree.getParent(i) == root);

        List<Integer> secondNodeIndices = new ArrayList<Integer>();
        List<Double> probabilities = new ArrayList<Double>();
        NodeRef j, iP, jP;
        iP = tree.getParent(i);
        double iParentHeight = tree.getNodeHeight(iP);
        double sum = 0.0;
        double backwardLikelihood = calculateTreeProbability(tree);
        int offset = (int) -backwardLikelihood;
        double backward = Math.exp(backwardLikelihood + offset);
        final NodeRef oldBrother = getOtherChild(tree, iP, i);
        final NodeRef oldGrandfather = tree.getParent(iP);

        tree.beginTreeEdit();
        for (int n = 0; n < nodeCount; n++) {
            j = tree.getNode(n);
            if (j != root) {
                jP = tree.getParent(j);

                if ((iP != jP)
                        && (tree.getNodeHeight(j) < iParentHeight && iParentHeight < tree
                        .getNodeHeight(jP))) {
                    secondNodeIndices.add(n);

                    pruneAndRegraft(tree, i, iP, j, jP);
                    double prob = Math.exp(calculateTreeProbability(tree)
                            + offset);
                    probabilities.add(prob);
                    sum += prob;

                    pruneAndRegraft(tree, i, iP, oldBrother, oldGrandfather);
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

        if (iP != jP) {
            pruneAndRegraft(tree, i, iP, j, jP);
            tree.pushTreeChangedEvent(i);
        }

        tree.endTreeEdit();

        // AR - not sure whether this check is necessary
        try {
            tree.checkTreeIsValid();
        } catch (InvalidTreeException e) {
            throw new OperatorFailedException(e.getMessage());
        }

        double forward = probabilities.get(index);

        // tree.pushTreeChangedEvent(jP);
        // tree.pushTreeChangedEvent(oldGrandfather);
        tree.pushTreeChangedEvent(i);

        double forwardProb = (forward / sum);
        double backwardProb = (backward / (sum - forward + backward));
        final double hastingsRatio = Math.log(backwardProb / forwardProb);

        return hastingsRatio;
    }

    private void pruneAndRegraft(TreeModel tree, NodeRef i, NodeRef iP,
                                 NodeRef j, NodeRef jP) throws OperatorFailedException {
        // tree.beginTreeEdit();

        // the grandfather
        NodeRef iG = tree.getParent(iP);
        // the brother
        NodeRef iB = getOtherChild(tree, iP, i);
        // prune
        tree.removeChild(iP, iB);
        tree.removeChild(iG, iP);
        tree.addChild(iG, iB);

        // reattach
        tree.removeChild(jP, j);
        tree.addChild(iP, j);
        tree.addChild(jP, iP);

    }

    private double calculateTreeProbability(Tree tree) {
        return probabilityEstimater.getTreeProbability(tree);
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
        return ImportancePruneAndRegraftParser.IMPORTANCE_PRUNE_AND_REGRAFT;
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
