/*
 * SubstitutionModelDelegate.java
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

package dr.evomodel.treedatalikelihood;

import beagle.Beagle;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.SubstitutionModel;

import java.io.Serializable;

/**
 * A simple substitution model delegate with the same substitution model over the whole tree
 * @author Andrew Rambaut
 * @version $Id$
 */
public final class HomogenousSubstitutionModelDelegate implements EvolutionaryProcessDelegate, Serializable {

    private final SubstitutionModel substitutionModel;

    private final int eigenCount = 1;
    private final int nodeCount;

    private final BufferIndexHelper eigenBufferHelper;
    private final BufferIndexHelper matrixBufferHelper;
    private final boolean cacheQMatrices;

    /**
     * A class which handles substitution models including epoch models where multiple
     * substitution models on a branch are convolved.
     * @param tree
     * @param branchModel Describes which substitution models use on each branch
     */
    public HomogenousSubstitutionModelDelegate(Tree tree, BranchModel branchModel) {
        this(tree, branchModel, 0, false);
    }

    public HomogenousSubstitutionModelDelegate(Tree tree, BranchModel branchModel, int partitionNumber) {
        this(tree, branchModel, partitionNumber, false);
    }

    /**
     * A class which handles substitution models including epoch models where multiple
     * substitution models on a branch are convolved.
     * @param tree
     * @param branchModel Describes which substitution models use on each branch
     * @param partitionNumber which data partition is this (used to offset eigen and matrix buffer numbers)
     */
    public HomogenousSubstitutionModelDelegate(Tree tree, BranchModel branchModel, int partitionNumber, boolean cacheQMatrices) {

        assert(branchModel.getSubstitutionModels().size() == 1) : "this delegate should only be used with simple branch models";

        this.substitutionModel = branchModel.getRootSubstitutionModel();
        
        nodeCount = tree.getNodeCount();

        // two eigen buffers for each decomposition for store and restore.
        eigenBufferHelper = new BufferIndexHelper(eigenCount, 0, partitionNumber);

        // two matrices for each node less the root
        matrixBufferHelper = new BufferIndexHelper(nodeCount, 0, partitionNumber);

        this.cacheQMatrices = cacheQMatrices;

    }// END: Constructor

    @Override
    public boolean canReturnComplexDiagonalization() {
        return substitutionModel.canReturnComplexDiagonalization();
    }

    @Override
    public int getEigenBufferCount() {
        return eigenBufferHelper.getBufferCount();
    }

    @Override
    public int getMatrixBufferCount() {
        return matrixBufferHelper.getBufferCount();
    }

    @Override
    public int getFirstOrderDifferentialMatrixBufferIndex(int branchIndex) {
        return matrixBufferHelper.getBufferCount() + getEigenIndex(0);
    }

    @Override
    public int getSecondOrderDifferentialMatrixBufferIndex(int branchIndex) {
        return matrixBufferHelper.getBufferCount() + getEigenBufferCount() + getEigenIndex(0);
    }

    @Override
    public void cacheFirstOrderDifferentialMatrix(Beagle beagle, int branchIndex, double[][] differentialMatrix) {
        double[] cacheMatrix = convertDifferentialMatrixToArray(differentialMatrix);
        beagle.setTransitionMatrix(getFirstOrderDifferentialMatrixBufferIndex(branchIndex), cacheMatrix, 0.0);
    }

    private double[] convertDifferentialMatrixToArray(double[][] differentialMatrix) {
        final int stateCount = substitutionModel.getDataType().getStateCount();
        assert(differentialMatrix.length == stateCount * stateCount);
        double[] cacheMatrix = new double[stateCount * stateCount];
        for (int i = 0; i < stateCount; i++) {
            System.arraycopy(differentialMatrix[i], 0, cacheMatrix, i * stateCount, stateCount);
        }
        return cacheMatrix;
    }

    @Override
    public void cacheSecondOrderDifferentialMatrix(Beagle beagle, int branchIndex, double[][] differentialMatrix) {
        double[] cacheMatrix = convertDifferentialMatrixToArray(differentialMatrix);
        beagle.setTransitionMatrix(getSecondOrderDifferentialMatrixBufferIndex(branchIndex), cacheMatrix, 0.0);
    }

    @Override
    public int getCachedMatrixBufferCount() {
        return 2 * getEigenBufferCount();
    }

    @Override
    public int getSubstitutionModelCount() {
        return 1;
    }

    @Override
    public SubstitutionModel getSubstitutionModel(int index) {
        assert(index == 0);
        return substitutionModel;
    }

    @Override
    public int getEigenIndex(int bufferIndex) {
        return eigenBufferHelper.getOffsetIndex(bufferIndex);
    }

    @Override
    public int getMatrixIndex(int branchIndex) {
        return matrixBufferHelper.getOffsetIndex(branchIndex);
    }

    @Override
    public double[] getRootStateFrequencies() {
        return substitutionModel.getFrequencyModel().getFrequencies();
    }

    @Override
    public void updateSubstitutionModels(Beagle beagle, boolean flip) {
        if (flip) {
            eigenBufferHelper.flipOffset(0);
        }
        EigenDecomposition ed = substitutionModel.getEigenDecomposition();

        beagle.setEigenDecomposition(
                eigenBufferHelper.getOffsetIndex(0),
                ed.getEigenVectors(),
                ed.getInverseEigenVectors(),
                ed.getEigenValues());

        // TODO Why is this here?  Are these matrices only needed right before `beagle.calculateEdgeDerivative()` is called?
        if (cacheQMatrices) {
            final int stateCount = substitutionModel.getDataType().getStateCount();
            double[] infinitesimalMatrix = new double[stateCount * stateCount];
            double[] infinitesimalMatrixSquared = new double[stateCount * stateCount];
            substitutionModel.getInfinitesimalMatrix(infinitesimalMatrix);
            beagle.setTransitionMatrix(getFirstOrderDifferentialMatrixBufferIndex(0), infinitesimalMatrix, 0.0);
            for (int l = 0; l < stateCount; l++) {
                for (int j = 0; j < stateCount; j++) {
                    double sumOverState = 0.0;
                    for (int k = 0; k < stateCount; k++) {
                        sumOverState += infinitesimalMatrix[l * stateCount + k] * infinitesimalMatrix[k * stateCount + j];
                    }
                    infinitesimalMatrixSquared[l * stateCount + j] = sumOverState;
                }
            }
            beagle.setTransitionMatrix(getSecondOrderDifferentialMatrixBufferIndex(0), infinitesimalMatrixSquared, 0.0);
        }
    }

    @Override
    public SubstitutionModel getSubstitutionModelForBranch(int branchIndex) {
        return substitutionModel;
    }

    @Override
    public void updateTransitionMatrices(Beagle beagle, int[] branchIndices, double[] edgeLengths, int updateCount, boolean flip) {

        int[] probabilityIndices = new int[updateCount];

        for (int i = 0; i < updateCount; i++) {
            if (flip) {
                matrixBufferHelper.flipOffset(branchIndices[i]);
            }
            probabilityIndices[i] = matrixBufferHelper.getOffsetIndex(branchIndices[i]);
        }// END: i loop

        beagle.updateTransitionMatrices(eigenBufferHelper.getOffsetIndex(0),
                probabilityIndices,
                null, // firstDerivativeIndices
                null, // secondDerivativeIndices
                edgeLengths,
                updateCount);

    }

    @Override
    public void flipTransitionMatrices(int[] branchIndices, int updateCount) {
        for (int i = 0; i < updateCount; i++) {
            matrixBufferHelper.flipOffset(branchIndices[i]);
        }
    }


    @Override
    public void storeState() {
        eigenBufferHelper.storeState();
        matrixBufferHelper.storeState();
    }

    @Override
    public void restoreState() {
        eigenBufferHelper.restoreState();
        matrixBufferHelper.restoreState();
    }

}// END: class
