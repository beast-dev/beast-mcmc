/*
 * ActionSubstitutionModelDelegate.java
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

package dr.evomodel.treedatalikelihood.action;

import beagle.Beagle;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import org.newejml.data.DMatrixSparseCSC;
import org.newejml.sparse.csc.CommonOps_DSCC;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class OldActionSubstitutionModelDelegate implements ActionEvolutionaryProcessDelegate{

    private final double[] branchLengths;
    private DMatrixSparseCSC[] sparseQs;
    private DMatrixSparseCSC[] scaledQs;
    private Tree tree;
    private final BranchModel branchModel;
    private final int nodeCount;
    private final int stateCount;

    public OldActionSubstitutionModelDelegate(Tree tree,
                                              BranchModel branchModel,
                                              int nodeCount) {
        this.branchLengths = new double[nodeCount];
        this.tree = tree;
        this.branchModel = branchModel;
        this.nodeCount = nodeCount;
        this.stateCount = branchModel.getRootSubstitutionModel().getFrequencyModel().getFrequencyCount();
        this.sparseQs = new DMatrixSparseCSC[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            sparseQs[i] = new DMatrixSparseCSC(stateCount, stateCount, 10 * stateCount);
        }
        this.scaledQs = new DMatrixSparseCSC[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            scaledQs[i] = new DMatrixSparseCSC(stateCount, stateCount, 10 * stateCount);
        }
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return false;
    }

    @Override
    public int getEigenBufferCount() {
        return 0;
    }

    @Override
    public int getMatrixBufferCount() {
        return 2 * nodeCount - 2;
    }

    @Override
    public int getInfinitesimalMatrixBufferIndex(int branchIndex) {
        return branchIndex;
    }

    @Override
    public int getInfinitesimalSquaredMatrixBufferIndex(int branchIndex) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public int getFirstOrderDifferentialMatrixBufferIndex(int branchIndex) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public int getSecondOrderDifferentialMatrixBufferIndex(int branchIndex) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public void cacheInfinitesimalMatrix(Beagle beagle, int bufferIndex, double[] differentialMatrix) {

    }

    @Override
    public void cacheInfinitesimalSquaredMatrix(Beagle beagle, int bufferIndex, double[] differentialMatrix) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public void cacheFirstOrderDifferentialMatrix(Beagle beagle, int branchIndex, double[] differentialMassMatrix) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public int getCachedMatrixBufferCount(PreOrderSettings settings) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public int getSubstitutionModelCount() {
        return branchModel.getSubstitutionModels().size();
    }

    @Override
    public SubstitutionModel getSubstitutionModel(int index) {
        return branchModel.getSubstitutionModels().get(index);
    }

    @Override
    public SubstitutionModel getSubstitutionModelForBranch(int branchIndex) {
        BranchModel.Mapping mapping = branchModel.getBranchModelMapping(tree.getNode(branchIndex));
        int[] order = mapping.getOrder();

        if (order.length > 1) {
            throw new RuntimeException("Not yet implemented");
        }

        return getSubstitutionModel(order[0]);
    }

    @Override
    public int getEigenIndex(int bufferIndex) {
        throw new RuntimeException("Eigen index should not be called for actions.");
    }

    @Override
    public int getMatrixIndex(int branchIndex) {
        return branchIndex;
    }

    @Override
    public double[] getRootStateFrequencies() {
        return branchModel.getRootFrequencyModel().getFrequencies();
    }

    @Override
    public void updateSubstitutionModels(Beagle beagle, boolean flipBuffers) {
        final int stateCount = branchModel.getRootSubstitutionModel().getFrequencyModel().getFrequencyCount();
        for (int i = 0; i < nodeCount; i++) {
            DMatrixSparseCSC sparseQ = sparseQs[i];
            SubstitutionModel substitutionModel = getSubstitutionModelForBranch(i);
            sparseQ.zero();
            double[] Q = new double[stateCount * stateCount];
            substitutionModel.getInfinitesimalMatrix(Q);
            for (int j = 0; j < stateCount; j++) {
                for (int k = 0; k < stateCount; k++) {
                    final double entryValue = Q[j * stateCount + k];
                    if (entryValue != 0.0) {
                        sparseQ.set(j, k, Q[j * stateCount + k]);
                    }
                }
            }
        }
    }

    @Override
    public void updateTransitionMatrices(Beagle beagle, int[] branchIndices, double[] edgeLengths, int updateCount, boolean flipBuffers) {
        for (int i = 0; i < updateCount; i++) {
            branchLengths[branchIndices[i]] = edgeLengths[i];
        }
    }

    @Override
    public void flipTransitionMatrices(int[] branchIndices, int updateCount) {
        throw new RuntimeException("Transition matrix is not used for actions.");
    }

    @Override
    public void storeState() {

    }

    @Override
    public void restoreState() {

    }

    @Override
    public DMatrixSparseCSC getScaledInstantaneousMatrix(int nodeIndex, double categoryRate) {
        CommonOps_DSCC.scale(branchLengths[nodeIndex] * categoryRate, sparseQs[nodeIndex], scaledQs[nodeIndex]);
        return scaledQs[nodeIndex];
    }
}
