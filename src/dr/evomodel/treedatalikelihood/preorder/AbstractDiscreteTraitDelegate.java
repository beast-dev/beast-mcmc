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
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.*;
import dr.inference.model.Model;

import java.util.Arrays;
import java.util.List;

import static dr.evolution.tree.TreeTrait.DA.factory;
import static dr.evomodel.treedatalikelihood.preorder.AbstractDiscreteTraitDelegate.MatrixChoice.GRADIENT;
import static dr.evomodel.treedatalikelihood.preorder.AbstractDiscreteTraitDelegate.MatrixChoice.HESSIAN;

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
        beagle.updatePrePartials(operations, operationCount, Beagle.NONE);  // Update all nodes with no rescaling
        //TODO:error control
    }

    @Override
    public void setupStatistics() {
        throw new RuntimeException("Not used (?) with BEAGLE");
    }

    @Override
    protected void simulateRoot(int rootNumber) {
        //This function sets preOrderPartials at Root for now.
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
        treeTraitHelper.addTrait(new TreeTrait.DA() {
            @Override
            public String getTraitName() {
                return "Hessian";
            }

            @Override
            public Intent getIntent() {
                return null;
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                return getHessian(tree, node);
            }
        });
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

    public enum MatrixChoice {
        GRADIENT {
            @Override
            public void getMatrix(SubstitutionModel model, double[] matrix) {
                model.getInfinitesimalMatrix(matrix);
            }

            @Override
            public double getBranchRateDifferential(ArbitraryBranchRates branchRates, double rate) {
                return branchRates.getBranchRateDifferential(rate);
            }

            @Override
            public double getRateDifferential(double rate) {
                return rate;
            }
        },
        HESSIAN {
            @Override
            public void getMatrix(SubstitutionModel model, double[] matrix) {
                double[] tmp = new double[matrix.length];
                model.getInfinitesimalMatrix(tmp);
                Arrays.fill(matrix, 0.0);

                final int stateCount = model.getDataType().getStateCount();
                for (int i = 0; i < stateCount; ++i){
                    for ( int j = 0; j < stateCount; ++j){
                        for ( int k = 0; k < stateCount; ++k){
                            matrix[i * stateCount + j] += tmp[i * stateCount + k] * tmp[k * stateCount + j];
                        }
                    }
                }
            }

            @Override
            public double getBranchRateDifferential(ArbitraryBranchRates branchRates, double rate) {
                // TODO
                return branchRates.getBranchRateSecondDifferential(rate);
            }

            @Override
            public double getRateDifferential(double rate) {
                return rate * rate;
            }
        };

        public abstract void getMatrix(SubstitutionModel model, double[] matrix);

        public abstract double getBranchRateDifferential(ArbitraryBranchRates branchRates, double rate);

        public abstract double getRateDifferential(double rate);

    }

    public double[] getHessian(Tree tree, NodeRef node) {
        final double[] patternGradient = getTrait(tree, node, GRADIENT); // this is 1/L * firstDerivative(L)
        final double[] patternDiagonalHessian  = getTrait(tree, node, HESSIAN); // this is 1/L * secondDerivative(L)
        final double[] patternWeights = patternList.getPatternWeights();
        double[] patternDiagonalLogHessian = new double[patternGradient.length];
        double[] diagonalLogHessian = new double[tree.getNodeCount() - 1];
        Arrays.fill(diagonalLogHessian, 0.0);

        for (int i = 0; i < patternGradient.length; ++i){
            patternDiagonalLogHessian[i] = patternDiagonalHessian[i] - patternGradient[i] * patternGradient[i];
        }

        int v = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); ++nodeNum){
            if (!tree.isRoot(tree.getNode(nodeNum))) {
                for (int pattern = 0; pattern < patternCount; pattern++) {
                    diagonalLogHessian[nodeNum] += patternDiagonalLogHessian[v * patternCount + pattern] * patternWeights[pattern];
                }
                v++;
            }
        }
        return diagonalLogHessian;
    }

    @Override
    public double[] getTrait(Tree tree, NodeRef node) {
        final double[] patternGradient = getTrait(tree, node, GRADIENT);
        final double[] patternWeights = patternList.getPatternWeights();
        double[] gradient = new double[tree.getNodeCount() - 1];
        Arrays.fill(gradient, 0.0);
        int v = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); ++nodeNum){
            if (!tree.isRoot(tree.getNode(nodeNum))) {
                for (int pattern = 0; pattern < patternCount; pattern++) {
                    gradient[nodeNum] += patternGradient[v * patternCount + pattern] * patternWeights[pattern];
                }
                v++;
            }
        }
        return gradient;
    }

    private double[] getTrait(Tree tree, NodeRef node, MatrixChoice matrixChoice) {

        assert (tree == this.tree);
        assert (node == null); // Implies: get trait for all nodes at same time

        //update all preOrder partials first
        simulationProcess.cacheSimulatedTraits(node);

        final double[] frequencies = evolutionaryProcessDelegate.getRootStateFrequencies();
        final double[] categoryWeights = siteRateModel.getCategoryProportions();
        final double[] categoryRates = siteRateModel.getCategoryRates();

        beagle.getPartials(getPostOrderPartialIndex(tree.getRoot().getNumber()), Beagle.NONE, rootPostOrderPartials);

        double[] derivative = new double[(tree.getNodeCount() - 1) * patternCount];

//        beagle.getSiteLogLikelihoods(cLikelihood);
        for (int category = 0; category < categoryCount; category++) {
            for (int pattern = 0; pattern < patternCount; pattern++) {

                final int patternIndex = category * patternCount + pattern;

                double sumOverEndState = 0.0;
                for (int k = 0; k < stateCount; k++) {
                    sumOverEndState += frequencies[k] * rootPostOrderPartials[patternIndex * stateCount + k];
                }

                cLikelihood[patternIndex] = sumOverEndState;
//                if (sumOverEndState < 1E-20){ // underflow occurred in postOrderPartials
//                    System.err.println("Underflow error occurred in postOrder Partials, try turn on scalingScheme=\"always\"");
//                }
            }
        }

        int index = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); ++nodeNum) {

            if (!tree.isRoot(tree.getNode(nodeNum))) {

                Arrays.fill(grandDenominator, 0.0);
                Arrays.fill(grandNumerator, 0.0);
                Arrays.fill(grandNumeratorIncrementLowerBound, 0.0);
                Arrays.fill(grandNumeratorIncrementUpperBound, 0.0);

                beagle.getPartials(getPostOrderPartialIndex(nodeNum), Beagle.NONE, postOrderPartial);
                beagle.getPartials(getPreOrderPartialIndex(nodeNum), Beagle.NONE, preOrderPartial);

                matrixChoice.getMatrix(evolutionaryProcessDelegate.getSubstitutionModelForBranch(nodeNum), Q);

                for (int category = 0; category < categoryCount; category++) {

                    final double weight = categoryWeights[category];
                    final double weightedRate = weight * matrixChoice.getRateDifferential(categoryRates[category]);

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

                        if (numerator != 0.0) {
                            if (denominator == 0.0) {
                                grandNumeratorIncrementLowerBound[pattern] += weightedRate * Math.min(numerator, 0.0);
                                grandNumeratorIncrementUpperBound[pattern] += weightedRate * Math.max(numerator, 0.0);
                            } else {
                                grandNumerator[pattern] += weightedRate * numerator *
                                        cLikelihood[patternOffset] / denominator;
                            }
                        }

                        grandDenominator[pattern] += weight * cLikelihood[patternOffset];
                    }
                }

                for (int pattern = 0; pattern < patternCount; pattern++) {

                    final double numerator = clampGradientNumerator(grandNumerator[pattern],
                            grandNumeratorIncrementLowerBound[pattern], grandNumeratorIncrementUpperBound[pattern]);

                    derivative[index * patternCount + pattern] += numerator / grandDenominator[pattern];

                    if (Double.isNaN(derivative[index]) && DEBUG) {
                        System.err.println("bad"); // OK, this should be invoked by underflow in lnL only now.
                    }
                }

//                final double branchLength = tree.getBranchLength(tree.getNode(nodeNum)); // TODO Delegate for rate models with multiple multiplier
//                gradient[index] *= branchLength;
                index++;
            }
        }

        return derivative;
    }

    // TODO Delegate for alternative behavior, like least value in magnitude
    private double clampGradientNumerator(double unbounded, double lowerBound, double upperBound) {
        return unbounded + (lowerBound + upperBound) / 2;
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
            operations[k++] = getPreOrderPartialIndex(tmpNodeOperation.getLeftChild());
            operations[k++] = Beagle.NONE;//getPreOrderScaleBufferIndex(tmpNodeOperation.getLeftChild());
            operations[k++] = Beagle.NONE;
            operations[k++] = getPreOrderPartialIndex(tmpNodeOperation.getNodeNumber());
            operations[k++] = evolutionaryProcessDelegate.getMatrixIndex(tmpNodeOperation.getLeftChild());
            operations[k++] = getPostOrderPartialIndex(tmpNodeOperation.getRightChild());
            operations[k++] = evolutionaryProcessDelegate.getMatrixIndex(tmpNodeOperation.getRightChild());
            // TODO Transform to Beagle operations here
        }
        return nodeOperations.size();
    }

    @Override
    public int getSingleOperationSize() {
        return Beagle.OPERATION_TUPLE_SIZE;
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

