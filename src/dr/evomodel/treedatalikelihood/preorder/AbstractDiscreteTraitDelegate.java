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
import beagle.InstanceDetails;
import dr.evolution.tree.*;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.treedatalikelihood.*;
import dr.inference.model.Model;

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

    private final BeagleDataLikelihoodDelegate likelihoodDelegate;
    private final Beagle beagle;
    private final EvolutionaryProcessDelegate evolutionaryProcessDelegate;
    private final SiteRateModel siteRateModel;
    private final int preOrderPartialOffset;

    // TODO Please use auto-format everywhere

    public AbstractDiscreteTraitDelegate(String name,
                                         Tree tree,
                                         BeagleDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree);
        this.likelihoodDelegate = likelihoodDelegate;
        this.beagle = likelihoodDelegate.getBeagleInstance();
        assert (this.likelihoodDelegate.isUsePreOrder()); /// TODO: reinitialize beagle instance if usePreOrder = false
        evolutionaryProcessDelegate = this.likelihoodDelegate.getEvolutionaryProcessDelegate();
        siteRateModel = this.likelihoodDelegate.getSiteRateModel();

        patternCount = likelihoodDelegate.getPatternList().getPatternCount();
        stateCount = likelihoodDelegate.getPatternList().getDataType().getStateCount();
        categoryCount = siteRateModel.getCategoryCount();

        // put preOrder partials right after postOrder partials
        preOrderPartialOffset = likelihoodDelegate.getPartialBufferCount();

    }

    @Override
    public void simulate(final int[] operations, final int operationCount,
                         final int rootNodeNumber) {
        //This function updates preOrder Partials for all nodes
        if (DEBUG) {
            System.err.println("Setting Root preOrder partial.");
        }

        this.simulateRoot(rootNodeNumber); //set up pre-order partials at root node first

        if (DEBUG) {
            System.err.println("Now update preOrder partials at all other nodes");
        }

        // TODO Use something like this: likelihoodDelegate.getPartialBufferIndex(0);

        int[] beagleoperations = new int[operationCount * Beagle.OPERATION_TUPLE_SIZE];
        for (int i = 0; i < operationCount; ++i) {
            beagleoperations[i * Beagle.OPERATION_TUPLE_SIZE] = getPreOrderPartialIndex(operations[i * 5]);
            beagleoperations[i * Beagle.OPERATION_TUPLE_SIZE + 1] = Beagle.NONE;
            beagleoperations[i * Beagle.OPERATION_TUPLE_SIZE + 2] = Beagle.NONE;
            beagleoperations[i * Beagle.OPERATION_TUPLE_SIZE + 3] = getPreOrderPartialIndex(operations[i * 5 + 1]);
            beagleoperations[i * Beagle.OPERATION_TUPLE_SIZE + 4] = evolutionaryProcessDelegate.getMatrixIndex(operations[i * 5 + 2]);
            beagleoperations[i * Beagle.OPERATION_TUPLE_SIZE + 5] = getPostOrderPartialIndex(operations[i * 5 + 3]);
            beagleoperations[i * Beagle.OPERATION_TUPLE_SIZE + 6] = evolutionaryProcessDelegate.getMatrixIndex(operations[i * 5 + 4]);
        }

        beagle.updatePrePartials(beagleoperations, operationCount, Beagle.NONE);  // Update all nodes with no rescaling
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

        final double[] postOrderPartial = new double[stateCount * patternCount * categoryCount];
        final double[] preOrderPartial = new double[stateCount * patternCount * categoryCount];
        final double[] frequencies = evolutionaryProcessDelegate.getRootStateFrequencies();
        final double[] rootPostOrderPartials = new double[stateCount * patternCount * categoryCount];

        //create a matrix for fetching the infinitesimal matrix Q
        double[] Q = new double[stateCount * stateCount];

        // TODO Should get from BeagleDataLikelihoodDelegate (already computed)
        double[] clikelihood = new double[categoryCount * patternCount];  // likelihood for each category doesn't come in free.

        beagle.getPartials(getPostOrderPartialIndex(tree.getRoot().getNumber()), Beagle.NONE, rootPostOrderPartials);

        double[] grandDenominator = new double[patternCount];
        double[] grandNumerator = new double[patternCount];

        double[] gradient = new double[tree.getNodeCount() - 1];
        double branchLength;

//        beagle.getSiteLogLikelihoods(clikelihood);
        for (int s = 0; s < categoryCount; s++) {
            for (int m = 0; m < patternCount; m++) {
                double clikelihoodTmp = 0;
                for (int k = 0; k < stateCount; k++) {
                    clikelihoodTmp += frequencies[k] * rootPostOrderPartials[s * patternCount * stateCount + m * stateCount + k];
                }
                clikelihood[s * patternCount + m] = clikelihoodTmp;
            }
        }

        int v = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); ++nodeNum) {
            if (!tree.isRoot(tree.getNode(nodeNum))) {
                for (int m = 0; m < patternCount; ++m) {
                    grandDenominator[m] = 0;
                    grandNumerator[m] = 0;
                }
                branchLength = tree.getBranchLength(tree.getNode(nodeNum));

                beagle.getPartials(getPostOrderPartialIndex(nodeNum), Beagle.NONE, postOrderPartial);
                beagle.getPartials(getPreOrderPartialIndex(nodeNum), Beagle.NONE, preOrderPartial);

                evolutionaryProcessDelegate.getSubstitutionModelForBranch(nodeNum).getInfinitesimalMatrix(Q);

                double[] tmpNumerator = new double[patternCount * categoryCount];

                //now calculate weights
                double denominator = 0;
                double numerator = 0;
                double tmp = 0;
                double[] weights = siteRateModel.getCategoryProportions();

                for (int s = 0; s < categoryCount; s++) {
                    double rs = siteRateModel.getRateForCategory(s);
                    double ws = weights[s];
                    for (int m = 0; m < patternCount; m++) {
                        numerator = 0;
                        denominator = 0;
                        for (int k = 0; k < stateCount; k++) {
                            tmp = 0;
                            for (int j = 0; j < stateCount; j++) {
                                tmp += Q[k * stateCount + j] * postOrderPartial[(s * patternCount + m) * stateCount + j];
                            }
                            numerator += tmp * preOrderPartial[(s * patternCount + m) * stateCount + k];
                            denominator += postOrderPartial[(s * patternCount + m) * stateCount + k] * preOrderPartial[(s * patternCount + m) * stateCount + k];
                        }
                        if (Double.isNaN(denominator)) {
                            System.err.println("bad bad");
                        }
                        tmpNumerator[s * patternCount + m] = ws * rs * numerator / denominator * clikelihood[s * patternCount + m];
                        grandNumerator[m] += tmpNumerator[s * patternCount + m];
                        grandDenominator[m] += ws * clikelihood[s * patternCount + m];
                    }
                }

                for (int m = 0; m < patternCount; m++) {
                    gradient[v] += grandNumerator[m] / grandDenominator[m] * likelihoodDelegate.getPatternList().getPatternWeight(m); // this assumes that root node is always at the end

                    if (Double.isNaN(gradient[v])) {
                        System.err.println("bad");
                    }
                }
                gradient[v] *= branchLength;
                v++;
            }
        }
        // See TipGradientViaFullConditionalDelegate.getTrait() as an example of using post- and pre-order partials together
        return gradient;
    }

    @Override
    public boolean isTraitLoggable() {
        return false;
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        // TODO When we start to cache intermediate calculation
    }

    @Override
    public void modelRestored(Model model) {

    }

    @Override
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int[] operations) {
        int k = 0;
        for (int i = 0; i < nodeOperations.size(); ++i) {
            NodeOperation tmpNodeOperation = nodeOperations.get(i);
            //nodeNumber = ParentNodeNumber, leftChild = nodeNumber, rightChild = siblingNodeNumber
            operations[k++] = tmpNodeOperation.getLeftChild();
            operations[k++] = tmpNodeOperation.getNodeNumber();
            operations[k++] = tmpNodeOperation.getLeftChild();
            operations[k++] = tmpNodeOperation.getRightChild();
            operations[k++] = tmpNodeOperation.getRightChild();
        }
        return nodeOperations.size();
    }

    private final int getPostOrderPartialIndex(final int nodeNumber) {
        return likelihoodDelegate.getPartialBufferIndex(nodeNumber);
    }

    private final int getPreOrderPartialIndex(final int nodeNumber) {
        return preOrderPartialOffset + nodeNumber;
    }


    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    private final int patternCount;
    private final int stateCount;
    private final int categoryCount;

    private static final boolean DEBUG = false;
}

