/*
 * DataSimulationDelegate.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.math.matrixAlgebra.WrappedVector;

import java.util.List;

/**
 * AbstractDiscreteTraitDelegate - interface for a plugin delegate for data simulation on a tree.
 *
 * @author Xiang Ji
 * @author Marc Suchard
 */
public abstract class AbstractDiscreteTraitDelegate extends ProcessSimulationDelegate.AbstractDelegate {

    private static final String GRADIENT_TRAIT_NAME = "Gradient";
    private static final String HESSIAN_TRAIT_NAME = "Hessian";

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

        this.patternList = likelihoodDelegate.getPatternList();
        this.gradient = new double[tree.getNodeCount() - 1];

        likelihoodDelegate.addModelListener(this);
        likelihoodDelegate.addModelRestoreListener(this);

        this.substitutionProcessKnown = false;
    }

    private void printMatrix(double[] matrix) {
        for (int rate = 0; rate < siteRateModel.getCategoryCount(); ++rate) {

            System.err.println("\nRate = " + rate);
            for (int i = 0; i < stateCount; ++i) {
                double[] row = new double[stateCount];
                System.arraycopy(matrix, rate * stateCount * stateCount + i * stateCount,
                        row, 0, stateCount);
                System.err.println(new WrappedVector.Raw(row));
            }
        }
    }

    private void debugMatrixTranspose(final int[] operations) {

        double[] matrix = new double[stateCount * stateCount * siteRateModel.getCategoryCount()];
        int matrixIndex = operations[4];
        beagle.getTransitionMatrix(matrixIndex, matrix);

        printMatrix(matrix);

        int destination = 1;

        beagle.transposeTransitionMatrices(new int[]{ matrixIndex }, new int[]{ destination }, 1);
        beagle.getTransitionMatrix(destination, matrix);

        printMatrix(matrix);
    }

    private static final boolean DEBUG_TRANSPOSE = false;

    @Override
    public void simulate(final int[] operations, final int operationCount,
                         final int rootNodeNumber) {
        //This function updates preOrder Partials for all nodes
        this.simulateRoot(rootNodeNumber);

        if (DEBUG_TRANSPOSE) { debugMatrixTranspose(operations); }

        beagle.updatePrePartials(operations, operationCount, Beagle.NONE);

        getNodeDerivatives(tree, gradient, null);

        if (COUNT_TOTAL_OPERATIONS) {
            ++simulateCount;
            updatePrePartialCount += operationCount;
        }
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

        beagle.setPartials(getPreOrderPartialIndex(rootNumber), rootPreOrderPartial);
    }

    @Override
    protected void simulateNode(int v0, int v1, int v2, int v3, int v4) {
        throw new RuntimeException("Not used with BEAGLE");
    }

    protected String getGradientTraitName() {
        return GRADIENT_TRAIT_NAME;
    }

    protected String getHessianTraitName() {
        return HESSIAN_TRAIT_NAME;
    }

    @Override
    protected void constructTraits(Helper treeTraitHelper) {

        treeTraitHelper.addTrait(new TreeTrait.DA() {
            @Override
            public String getTraitName() {
                return getGradientTraitName();
            }

            @Override
            public Intent getIntent() {
                return Intent.BRANCH;
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                return getGradient(node);
            }
        });

        treeTraitHelper.addTrait(new TreeTrait.DA() {
            @Override
            public String getTraitName() {
                return getHessianTraitName();
            }

            @Override
            public Intent getIntent() {
                return Intent.BRANCH;
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                return getHessian(tree, node);
            }
        });
    }

    private double[] getHessian(Tree tree, NodeRef node) {

        //update all preOrder partials first
        simulationProcess.cacheSimulatedTraits(node);

        double[] second = new double[tree.getNodeCount() - 1];
        getNodeDerivatives(tree, null, second);

        return second;
    }

    private double[] getGradient(NodeRef node) {

        if (COUNT_TOTAL_OPERATIONS) {
            ++getTraitCount;
        }

        //update all preOrder partials first
        simulationProcess.cacheSimulatedTraits(node);
        return gradient.clone();
    }

    abstract protected void cacheDifferentialMassMatrix(Tree tree, boolean cacheSquaredMatrix);

    private static final boolean USE_CACHE = true;

    private void getNodeDerivatives(Tree tree, double[] first, double[] second) {

        final int[] postBufferIndices = new int[tree.getNodeCount() - 1];
        final int[] preBufferIndices = new int[tree.getNodeCount() - 1];
        final int[] firstDervIndices = new int[tree.getNodeCount() - 1];
        final int[] secondDeriveIndices = new int[tree.getNodeCount() - 1];

        boolean needsUpdate = !substitutionProcessKnown || second != null;
        if (!USE_CACHE || needsUpdate) {
            cacheDifferentialMassMatrix(tree, second != null);
            substitutionProcessKnown = true;
        }

        int u = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); nodeNum++) {
            if (!tree.isRoot(tree.getNode(nodeNum))) {
                postBufferIndices[u] = getPostOrderPartialIndex(nodeNum);
                preBufferIndices[u]  = getPreOrderPartialIndex(nodeNum);
                firstDervIndices[u]  = getFirstDerivativeMatrixBufferIndex(nodeNum);
                secondDeriveIndices[u] = getSecondDerivativeMatrixBufferIndex(nodeNum);
                u++;
            }
        }

        double[] firstSquared = (second != null) ? new double[second.length] : null;

        beagle.calculateEdgeDifferentials(postBufferIndices, preBufferIndices,
                firstDervIndices, new int[] { 0 }, tree.getNodeCount() - 1,
                null, first, firstSquared);

        if (second != null) {
            beagle.calculateEdgeDifferentials(postBufferIndices, preBufferIndices,
                    secondDeriveIndices, new int[] { 0 }, tree.getNodeCount() - 1,
                    null, second, null);

            for (int i = 0; i < second.length; ++i) {
                second[i] -= firstSquared[i];
            }
        }
    }

    protected int getFirstDerivativeMatrixBufferIndex(int nodeNum) {
        return evolutionaryProcessDelegate.getInfinitesimalMatrixBufferIndex(nodeNum);
    }

    protected int getSecondDerivativeMatrixBufferIndex(int nodeNum) {
        return evolutionaryProcessDelegate.getInfinitesimalSquaredMatrixBufferIndex(nodeNum);
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        substitutionProcessKnown = false;
    }

    @Override
    public void modelRestored(Model model) {
        substitutionProcessKnown = false;
    }

    @Override
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int[] operations) {
        int k = 0;
        for (NodeOperation tmpNodeOperation : nodeOperations) {
            //nodeNumber = ParentNodeNumber, leftChild = nodeNumber, rightChild = siblingNodeNumber
            operations[k++] = getPreOrderPartialIndex(tmpNodeOperation.getLeftChild());
            operations[k++] = Beagle.NONE;
            operations[k++] = Beagle.NONE;
            operations[k++] = getPreOrderPartialIndex(tmpNodeOperation.getNodeNumber());
            operations[k++] = evolutionaryProcessDelegate.getMatrixIndex(tmpNodeOperation.getLeftChild());
            operations[k++] = getPostOrderPartialIndex(tmpNodeOperation.getRightChild());
            operations[k++] = evolutionaryProcessDelegate.getMatrixIndex(tmpNodeOperation.getRightChild());
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

    @Override
    public String toString() {
        if (COUNT_TOTAL_OPERATIONS) {

            return "\tsimulateCount = " + simulateCount + "\n" +
                    "\tgetTraitCount = " + getTraitCount + "\n" +
                    "\tupPrePartialCount = " + updatePrePartialCount + "\n";

        } else {
            return super.toString();
        }
    }
    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    protected final BeagleDataLikelihoodDelegate likelihoodDelegate;
    protected final Beagle beagle;
    protected EvolutionaryProcessDelegate evolutionaryProcessDelegate;
    protected final SiteRateModel siteRateModel;
    protected final PatternList patternList;

    protected final int patternCount;
    protected final int stateCount;
    protected final int categoryCount;
    private int preOrderPartialOffset;

    protected final double[] gradient;

    private boolean substitutionProcessKnown;

    private static final boolean COUNT_TOTAL_OPERATIONS = true;
    private long simulateCount = 0;
    private long getTraitCount = 0;
    private long updatePrePartialCount = 0;
}
