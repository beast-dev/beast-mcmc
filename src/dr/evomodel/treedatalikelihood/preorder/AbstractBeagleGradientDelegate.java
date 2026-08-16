/*
 * AbstractBeagleGradientDelegate.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.treedatalikelihood.preorder;

import beagle.Beagle;
import dr.evolution.alignment.PatternList;
import dr.evolution.tree.*;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.treedatalikelihood.*;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.PreOrderMessageProvider;
import dr.inference.model.Model;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * AbstractBeagleGradientDelegate - interface for a plugin delegate for data simulation on a tree.
 *
 * @author Xiang Ji
 * @author Marc Suchard
 */
public abstract class AbstractBeagleGradientDelegate extends ProcessSimulationDelegate.AbstractDelegate
        implements PreOrderMessageProvider {

    private static final String GRADIENT_TRAIT_NAME = "Gradient";
    private static final String HESSIAN_TRAIT_NAME = "Hessian";

    protected AbstractBeagleGradientDelegate(String name,
                                             Tree tree,
                                             BeagleDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree);
        this.tree = tree;
        this.likelihoodDelegate = likelihoodDelegate;
        this.beagle = likelihoodDelegate.getBeagleInstance();
        this.preOrderType = getPreOrderType();

        assert (this.likelihoodDelegate.isUsePreOrder()); /// TODO: reinitialize beagle instance if usePreOrder = false

        this.evolutionaryProcessDelegate = likelihoodDelegate.getEvolutionaryProcessDelegate();
        this.siteRateModel = likelihoodDelegate.getSiteRateModel();

        this.patternCount = likelihoodDelegate.getPatternList().getPatternCount();
        this.stateCount = likelihoodDelegate.getPatternList().getDataType().getStateCount();
        this.categoryCount = siteRateModel.getCategoryCount();

        // put preOrder partials right after postOrder partials
        this.preOrderPartialOffset = likelihoodDelegate.getPartialBufferCount();

        this.patternList = likelihoodDelegate.getPatternList();

        likelihoodDelegate.addModelListener(this);
        likelihoodDelegate.addModelRestoreListener(this);

        this.substitutionProcessKnown = false;
    }

    @Override
    public int getStateCount() { return stateCount; }

    @Override
    public int getPatternCount() { return patternCount; }

    @Override
    public int getCategoryCount() { return categoryCount; }

    @Override
    public void getPreorderPartials(int nodeNumber, DiscretePartialsType type, double[] out) {

        assert out.length >= categoryCount * patternCount * stateCount;

        if (type == preOrderType) {
            beagle.getPartials(getPreOrderPartialIndex(nodeNumber), Beagle.NONE, out);
        } else {
            Arrays.fill(out, 0.0);
        }
    }

    protected DiscretePartialsType getPreOrderType() {
        return DiscretePartialsType.BOTTOM;
    }

    abstract protected int getGradientLength();

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

        beagle.updatePrePartials_v5(operations, operationCount, Beagle.NONE, preOrderType.getBeagleType());

        if (gradient == null) {
            gradient = new double[getGradientLength()];
        }

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

    double[] getHessian(Tree tree, NodeRef node) {

        //update all preOrder partials first
        simulationProcess.cacheSimulatedTraits(node);

        double[] second = new double[getGradientLength()];
        getNodeDerivatives(tree, null, second);

        return second;
    }

    protected double[] getGradient(NodeRef node) {

        if (COUNT_TOTAL_OPERATIONS) {
            ++getTraitCount;
        }

        //update all preOrder partials first
        simulationProcess.cacheSimulatedTraits(node);
        return gradient.clone();
    }

    abstract protected void getNodeDerivatives(Tree tree, double[] first, double[] second);

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
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int rootNodeNumber, int[] operations) {
        if (preOrderType == DiscretePartialsType.BOTTOM) {
            return vectorizeNodeOperationsBottom(nodeOperations, operations);
        } else if (preOrderType == DiscretePartialsType.TOP) {
            return vectorizeNodeOperationsTop(nodeOperations, rootNodeNumber, operations);
        } else if (preOrderType == DiscretePartialsType.BOTTOM_SPECTRAL) {
            return vectorizeNodeOperationsBottom(nodeOperations, operations);
        } else if (preOrderType == DiscretePartialsType.TOP_SPECTRAL) {
            return vectorizeNodeOperationsTop(nodeOperations, rootNodeNumber, operations);
        } else {
            throw new RuntimeException("Not yet implemented");
        }
    }

    private int vectorizeNodeOperationsBottom(List<NodeOperation> nodeOperations, int[] operations) {
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

    private int vectorizeNodeOperationsTop(List<NodeOperation> nodeOperations, int rootNodeNumber, int[] operations) {
        final boolean scaling = likelihoodDelegate.isUsingScaleFactors();
        if (scaling) {
            initializeScaleBufferTracking();
        }

        int k = 0;
        for (NodeOperation tmpNodeOperation : nodeOperations) {
            //nodeNumber = ParentNodeNumber, leftChild = nodeNumber, rightChild = siblingNodeNumber
            final int nodeNumber = tmpNodeOperation.getLeftChild();
            final int parentNumber = tmpNodeOperation.getNodeNumber();
            final int siblingNumber = tmpNodeOperation.getRightChild();

            operations[k++] = getPreOrderPartialIndex(nodeNumber);
            operations[k++] = Beagle.NONE;
            operations[k++] = Beagle.NONE;
            operations[k++] = getPreOrderPartialIndex(parentNumber);
            operations[k++] = (parentNumber == rootNodeNumber) ? Beagle.NONE :
                    evolutionaryProcessDelegate.getMatrixIndex(parentNumber);
            operations[k++] = getPostOrderPartialIndex(siblingNumber);
            operations[k++] = evolutionaryProcessDelegate.getMatrixIndex(siblingNumber);

            if (scaling) {
                buildPreCumulativeScaleBuffer(nodeNumber, parentNumber, siblingNumber, rootNodeNumber);
            }
        }
        return nodeOperations.size();
    }

    /**
     * Builds, for every node touched by a pre-order traversal, the two
     * rescaling-correction buffers beagleCalculateAdjointCrossProductDerivative
     * needs (see its contract in beagle.h) -- a node's CUMULATIVE post-order
     * scale (its own scale-write plus every rescaled descendant's) and the
     * CUMULATIVE pre-order scale already baked into its pre-order partial by
     * every sibling subtree incorporated while constructing it.
     *
     * Both are tracked as flat lists of RAW individual scale-write buffers
     * (never as a previously-built cumulative buffer fed back into another
     * beagleAccumulateScaleFactors call) because accumulateScaleFactors
     * always LOG-transforms its inputs and treats its target as already-LOG:
     * feeding it an already-cumulative buffer as an *input* would silently
     * re-apply log() to an already-logged value and corrupt the sum.
     */
    private void initializeScaleBufferTracking() {
        final int nodeCount = tree.getNodeCount();
        if (scaleBufferBase < 0) {
            // BeagleDataLikelihoodDelegate reserves 2*nodeCount extra scale
            // buffers (beyond its own doubled/store-restore pool) specifically
            // for this: [base, base+nodeCount) for post-cumulative buffers
            // indexed directly by node number, [base+nodeCount, base+2*nodeCount)
            // for pre-cumulative buffers. Root's slots in each half go unused.
            scaleBufferBase = likelihoodDelegate.getScaleBufferCount();
        }
        postCumulativeRawList = new List[nodeCount];
        preCumulativeRawList = new List[nodeCount];
        postCumulativeScaleBufferIndex = new int[nodeCount];
        preCumulativeScaleBufferIndex = new int[nodeCount];
        Arrays.fill(postCumulativeScaleBufferIndex, Beagle.NONE);
        Arrays.fill(preCumulativeScaleBufferIndex, Beagle.NONE);
    }

    private List<Integer> getPostCumulativeRawList(int nodeNumber) {
        if (postCumulativeRawList[nodeNumber] == null) {
            List<Integer> list = new ArrayList<>();
            collectRescaledIndividualBuffers(nodeNumber, list);
            postCumulativeRawList[nodeNumber] = list;
            if (!list.isEmpty()) {
                int slot = scaleBufferBase + nodeNumber;
                beagle.resetScaleFactors(slot);
                beagle.accumulateScaleFactors(toIntArray(list), list.size(), slot);
                postCumulativeScaleBufferIndex[nodeNumber] = slot;
            }
        }
        return postCumulativeRawList[nodeNumber];
    }

    private void collectRescaledIndividualBuffers(int nodeNumber, List<Integer> list) {
        NodeRef node = tree.getNode(nodeNumber);
        if (tree.isExternal(node)) {
            return; // tips are never independently rescaled
        }
        int individual = likelihoodDelegate.getNodeIndividualScaleBufferIndex(nodeNumber);
        if (individual != Beagle.NONE) {
            list.add(individual);
        }
        for (int i = 0; i < tree.getChildCount(node); ++i) {
            collectRescaledIndividualBuffers(tree.getChild(node, i).getNumber(), list);
        }
    }

    private void buildPreCumulativeScaleBuffer(int nodeNumber, int parentNumber, int siblingNumber, int rootNodeNumber) {
        List<Integer> list = new ArrayList<>();
        if (parentNumber != rootNodeNumber) {
            // parentNumber was already processed earlier in this same top-down
            // traversal (its pre-order buffer must already exist for THIS
            // node's beagle.updatePrePartials_v5 operation to reference it).
            list.addAll(preCumulativeRawList[parentNumber]);
        }
        list.addAll(getPostCumulativeRawList(siblingNumber));
        preCumulativeRawList[nodeNumber] = list;
        if (!list.isEmpty()) {
            int slot = scaleBufferBase + tree.getNodeCount() + nodeNumber;
            beagle.resetScaleFactors(slot);
            beagle.accumulateScaleFactors(toIntArray(list), list.size(), slot);
            preCumulativeScaleBufferIndex[nodeNumber] = slot;
        }
    }

    private static int[] toIntArray(List<Integer> list) {
        int[] array = new int[list.size()];
        for (int i = 0; i < array.length; ++i) {
            array[i] = list.get(i);
        }
        return array;
    }

    /** @return the cumulative (LOG) post-order scale buffer for nodeNumber's
     * OWN partial (its own scale-write plus every rescaled descendant's), or
     * Beagle.NONE if nothing at-or-below it was ever independently rescaled.
     * Only valid after simulate() has run for the current tree/model state. */
    protected int getPostCumulativeScaleBufferIndex(int nodeNumber) {
        return (postCumulativeScaleBufferIndex == null) ? Beagle.NONE : postCumulativeScaleBufferIndex[nodeNumber];
    }

    /** @return the cumulative (LOG) rescaling already baked into nodeNumber's
     * pre-order partial, or Beagle.NONE if it carries no such scaling. Only
     * valid after simulate() has run for the current tree/model state. */
    protected int getPreCumulativeScaleBufferIndex(int nodeNumber) {
        return (preCumulativeScaleBufferIndex == null) ? Beagle.NONE : preCumulativeScaleBufferIndex[nodeNumber];
    }

    /** @return the cumulative (LOG) scale buffer used to correct the last
     * likelihood evaluation's logL, or Beagle.NONE if rescaling was not applied. */
    protected int getRootCumulativeScaleBufferIndex() {
        return likelihoodDelegate.getRootCumulativeScaleBufferIndex();
    }

    @Override
    public int getSingleOperationSize() {
        return Beagle.OPERATION_TUPLE_SIZE;
    }

    protected int getPostOrderPartialIndex(final int nodeNumber) {
        return likelihoodDelegate.getPartialBufferIndex(nodeNumber);
    }

    protected int getPreOrderPartialIndex(final int nodeNumber) {
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
    protected final DiscretePartialsType preOrderType;

    protected final int patternCount;
    protected final int stateCount;
    protected final int categoryCount;
    private final int preOrderPartialOffset;

    protected double[] gradient;

    protected boolean substitutionProcessKnown;
    protected Tree tree;

    // Rescaling-correction bookkeeping for beagleCalculateAdjointCrossProductDerivative
    // (see initializeScaleBufferTracking); all by node number, rebuilt every simulate() call.
    private int scaleBufferBase = -1;
    private List<Integer>[] postCumulativeRawList;
    private List<Integer>[] preCumulativeRawList;
    private int[] postCumulativeScaleBufferIndex;
    private int[] preCumulativeScaleBufferIndex;

    private static final boolean COUNT_TOTAL_OPERATIONS = true;
    final boolean DEBUG = false;
    private long simulateCount = 0;
    private long getTraitCount = 0;
    private long updatePrePartialCount = 0;
}
