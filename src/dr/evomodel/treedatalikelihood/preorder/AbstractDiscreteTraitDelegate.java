/*
 * DataSimulationDelegate.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.preorder;

import beagle.Beagle;
import dr.evolution.alignment.PatternList;
import dr.evolution.tree.*;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.treedatalikelihood.*;
import dr.inference.model.Model;

import java.util.Arrays;
import java.util.List;

import static dr.evolution.tree.TreeTrait.DA.factory;

/**
 * AbstractDiscreteTraitDelegate - interface for a plugin delegate for data simulation on a tree.
 *
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class AbstractDiscreteTraitDelegate extends ProcessSimulationDelegate.AbstractDelegate
        implements TreeTrait.TraitInfo<double[]> {

    public AbstractDiscreteTraitDelegate(String name,
                                         Tree tree,
                                         BeagleDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree);
        this.likelihoodDelegate = likelihoodDelegate;
        this.beagle = likelihoodDelegate.getBeagleInstance();

        assert (this.likelihoodDelegate.isUsePreOrder()); /// TODO: reinitialize beagle instance if usePreOrder = false

        this.evolutionaryProcessDelegate = likelihoodDelegate.getEvolutionaryProcessDelegate();
        this.siteRateModel = likelihoodDelegate.getSiteRateModel();

        this.patternCount = likelihoodDelegate.getPatternList().getPatternCount();
        this.stateCount = likelihoodDelegate.getPatternList().getDataType().getStateCount();
        this.categoryCount = siteRateModel.getCategoryCount();

        // put preOrder partials right after postOrder partials
        this.preOrderPartialOffset = likelihoodDelegate.getPartialBufferCount();

        // put scaleBuffers for preOrder partials right after those for postOrder partials
        this.preOrderScaleBufferOffset = tree.getNodeCount() - tree.getExternalNodeCount() - 1;

        this.patternList = likelihoodDelegate.getPatternList();

        this.postOrderPartial = new double[stateCount * patternCount * categoryCount];
        this.preOrderPartial = new double[stateCount * patternCount * categoryCount];
        this.rootPostOrderPartials = new double[stateCount * patternCount * categoryCount];
        this.Q = new double[stateCount * stateCount];
        this.cLikelihood = new double[categoryCount * patternCount];

        this.grandDenominator = new double[patternCount];
        this.grandNumerator = new double[patternCount];
        this.grandNumeratorIncrementLowerBound = new double[patternCount];
        this.grandNumeratorIncrementUpperBound = new double[patternCount];

    }

    @Override
    public void simulate(final int[] operations, final int operationCount,
                         final int rootNodeNumber) {
        //This function updates preOrder Partials for all nodes
        this.simulateRoot(rootNodeNumber); //set up pre-order partials at root node first

        // TODO Use something like this: likelihoodDelegate.getPartialBufferIndex(0);
        int[] beagleOperations = new int[operationCount * Beagle.OPERATION_TUPLE_SIZE];
        for (int i = 0; i < operationCount; ++i) {
            beagleOperations[i * Beagle.OPERATION_TUPLE_SIZE] = getPreOrderPartialIndex(operations[i * 5]);
            beagleOperations[i * Beagle.OPERATION_TUPLE_SIZE + 1] = getPreOrderScaleBufferIndex(operations[i * 5]);
            beagleOperations[i * Beagle.OPERATION_TUPLE_SIZE + 2] = Beagle.NONE;
            beagleOperations[i * Beagle.OPERATION_TUPLE_SIZE + 3] = getPreOrderPartialIndex(operations[i * 5 + 1]);
            beagleOperations[i * Beagle.OPERATION_TUPLE_SIZE + 4] = evolutionaryProcessDelegate.getMatrixIndex(operations[i * 5 + 2]);
            beagleOperations[i * Beagle.OPERATION_TUPLE_SIZE + 5] = getPostOrderPartialIndex(operations[i * 5 + 3]);
            beagleOperations[i * Beagle.OPERATION_TUPLE_SIZE + 6] = evolutionaryProcessDelegate.getMatrixIndex(operations[i * 5 + 4]);
            // TODO Shouldn't this transformation happen in vectorizeNodeOperations() such that we only transform once?
        }

        beagle.updatePrePartials(beagleOperations, operationCount, Beagle.NONE);  // Update all nodes with no rescaling
        //TODO:error control
    }

    @Override
    public void setupStatistics() {
        throw new RuntimeException("Not used (?) with BEAGLE");
    }

    @Override
    protected void simulateRoot(int rootNumber) {
        //This function sets preOrderPartials at Root for now.
        if (DEBUG) {
            System.err.println("Simulate root node " + rootNumber);
        }
        final double[] frequencies = evolutionaryProcessDelegate.getRootStateFrequencies();
        double[] rootPreOrderPartial = new double[stateCount * patternCount * categoryCount];
        for (int i = 0; i < patternCount * categoryCount; ++i) {
            System.arraycopy(frequencies, 0, rootPreOrderPartial, i * stateCount, stateCount);
        }
//        InstanceDetails instanceDetails = beagle.getDetails();

        beagle.setPartials(getPreOrderPartialIndex(rootNumber), rootPreOrderPartial); // Call beagle.setPartials()
        //TODO: find the right error message for control
    }

    @Override
    protected void simulateNode(int v0, int v1, int v2, int v3, int v4) {
        throw new RuntimeException("Not used with BEAGLE");
    }

    @Override
    protected void constructTraits(Helper treeTraitHelper) {
        treeTraitHelper.addTrait(factory(this));
    }

    public static String getName(String name) {
        return "derivative." + name;
    }

    @Override
    public String getTraitName() {
        return getName(name);
    }

    @Override
    public TreeTrait.Intent getTraitIntent() {
        return TreeTrait.Intent.NODE;
    }

    @Override
    public Class getTraitClass() {
        return double[].class;
    }

    @Override
    public double[] getTrait(Tree tree, NodeRef node) {

        assert (tree == this.tree);
        assert (node == null); // Implies: get trait for all nodes at same time

        //update all preOrder partials first
        simulationProcess.cacheSimulatedTraits(node);

        final double[] frequencies = evolutionaryProcessDelegate.getRootStateFrequencies();
        final double[] patternWeights = patternList.getPatternWeights();
        final double[] categoryWeights = siteRateModel.getCategoryProportions();
        final double[] categoryRates = siteRateModel.getCategoryRates();

        // TODO Should get from BeagleDataLikelihoodDelegate (already computed)
        beagle.getPartials(getPostOrderPartialIndex(tree.getRoot().getNumber()), Beagle.NONE, rootPostOrderPartials);

        double[] gradient = new double[tree.getNodeCount() - 1];

//        beagle.getSiteLogLikelihoods(cLikelihood);
        for (int category = 0; category < categoryCount; category++) {
            for (int pattern = 0; pattern < patternCount; pattern++) {

                final int patternIndex = category * patternCount + pattern;

                double sumOverEndState = 0.0;
                for (int k = 0; k < stateCount; k++) {
                    sumOverEndState += frequencies[k] * rootPostOrderPartials[patternIndex * stateCount + k];
                }

                cLikelihood[category * patternCount + pattern] = sumOverEndState;
//                if (sumOverEndState < 1E-20){ // underflow occurred in postOrderPartials
//                    System.err.println("Underflow error occurred in postOrder Partials, try turn on scalingScheme=\"always\"");
//                    if(DEBUG){
//                        System.err.println("underflow occurred in postOrder Partials");
//                    }
//                }
            }
        }

        int index = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); ++nodeNum) {

            if (!tree.isRoot(tree.getNode(nodeNum))) {

                Arrays.fill(grandDenominator, 0.0);
                Arrays.fill(grandNumerator, 0.0);

                beagle.getPartials(getPostOrderPartialIndex(nodeNum), Beagle.NONE, postOrderPartial);
                beagle.getPartials(getPreOrderPartialIndex(nodeNum), Beagle.NONE, preOrderPartial);

                evolutionaryProcessDelegate.getSubstitutionModelForBranch(nodeNum).getInfinitesimalMatrix(Q);

                for (int category = 0; category < categoryCount; category++) {

                    final double rate = categoryRates[category];
                    final double weight = categoryWeights[category];

                    for (int pattern = 0; pattern < patternCount; pattern++) {

                        final int patternOffset = category * patternCount + pattern;

                        double numerator = 0.0;
                        double denominator = 0.0;

                        for (int k = 0; k < stateCount; k++) {

                            double sumOverEndState = 0;
                            for (int j = 0; j < stateCount; j++) {
                                sumOverEndState += Q[k * stateCount + j]
                                        * postOrderPartial[patternOffset * stateCount + j];
                            }

                            numerator += sumOverEndState * preOrderPartial[patternOffset * stateCount + k];
                            denominator += postOrderPartial[patternOffset * stateCount + k]
                                    * preOrderPartial[patternOffset * stateCount + k];
                        }

                        double grandNumeratorIncrement = weight * rate * numerator / denominator * cLikelihood[patternOffset];
                        if (denominator == 0) {  // Now instead evaluate the bound of the gradient
                            grandNumeratorIncrement = 0.0;  // if numerator == 0, it is 0
                            if (numerator != 0.0) {
                                grandNumeratorIncrementLowerBound[pattern] += weight * rate * (numerator > 0 ? 0.0 : numerator);
                                grandNumeratorIncrementUpperBound[pattern] += weight * rate * (numerator > 0 ? numerator : 0.0);
                            }
                        } else if (Double.isNaN(denominator)) {
                            System.err.println("something wrong with preOrder partial calculation.");
                        }
                        grandNumerator[pattern] += grandNumeratorIncrement;
                        grandDenominator[pattern] += weight * cLikelihood[patternOffset];
                    }
                }

                for (int pattern = 0; pattern < patternCount; pattern++) {

                    gradient[index] += (grandNumerator[pattern] + (grandNumeratorIncrementLowerBound[pattern] + grandNumeratorIncrementUpperBound[pattern]) / 2.0)
                            / grandDenominator[pattern] * patternWeights[pattern];

                    if (Double.isNaN(gradient[index])) {
                        System.err.println("bad");
                    }
                }

                final double branchLength = tree.getBranchLength(tree.getNode(nodeNum)); // TODO Delegate for rate models with multiple multipliers
                gradient[index] *= branchLength;
                index++;
            }
        }

        return gradient;
    }

    @Override
    public boolean isTraitLoggable() {
        return false;
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        // Do nothing
    }

    @Override
    public void modelRestored(Model model) {
        // Do nothing
    }

    @Override
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int[] operations) {
        int k = 0;
        for (NodeOperation tmpNodeOperation : nodeOperations) {
            //nodeNumber = ParentNodeNumber, leftChild = nodeNumber, rightChild = siblingNodeNumber
            operations[k++] = tmpNodeOperation.getLeftChild();
            operations[k++] = tmpNodeOperation.getNodeNumber();
            operations[k++] = tmpNodeOperation.getLeftChild();
            operations[k++] = tmpNodeOperation.getRightChild();
            operations[k++] = tmpNodeOperation.getRightChild();
            // TODO Transform to Beagle operations here
        }
        return nodeOperations.size();
    }

    private int getPostOrderPartialIndex(final int nodeNumber) {
        return likelihoodDelegate.getPartialBufferIndex(nodeNumber);
    }

    private int getPreOrderPartialIndex(final int nodeNumber) {
        return preOrderPartialOffset + nodeNumber;
    }

    private int getPreOrderScaleBufferIndex(final int nodeNumber) {
        return preOrderScaleBufferOffset + nodeNumber;
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    private final BeagleDataLikelihoodDelegate likelihoodDelegate;
    private final Beagle beagle;
    private final EvolutionaryProcessDelegate evolutionaryProcessDelegate;
    private final SiteRateModel siteRateModel;
    private final PatternList patternList;

    private final int patternCount;
    private final int stateCount;
    private final int categoryCount;
    private final int preOrderPartialOffset;
    private final int preOrderScaleBufferOffset;

    private final double[] postOrderPartial;
    private final double[] preOrderPartial;
    private final double[] rootPostOrderPartials;
    private final double[] Q;
    private final double[] cLikelihood;

    private final double[] grandDenominator;
    private final double[] grandNumerator;
    private final double[] grandNumeratorIncrementLowerBound;
    private final double[] grandNumeratorIncrementUpperBound;

    private static final boolean DEBUG = false;
}

