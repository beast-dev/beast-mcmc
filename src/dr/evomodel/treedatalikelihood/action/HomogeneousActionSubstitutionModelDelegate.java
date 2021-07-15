/*
 * HomogeneousActionSubstitutionModelDelegate.java
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

package dr.evomodel.treedatalikelihood.action;

import beagle.Beagle;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import org.newejml.data.DMatrixSparseCSC;
import org.newejml.data.DMatrixSparseTriplet;
import org.newejml.ops.DConvertMatrixStruct;
import org.newejml.sparse.csc.CommonOps_DSCC;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class HomogeneousActionSubstitutionModelDelegate implements ActionEvolutionaryProcessDelegate {

    private final SubstitutionModel substitutionModel;
    private final double[] branchLengths;
    private DMatrixSparseTriplet sparseQ;

    public HomogeneousActionSubstitutionModelDelegate(SubstitutionModel substitutionModel,
                                                      int nodeCount) {
//        assert(branchModel.getSubstitutionModels().size() == 1) : "this delegate should only be used with simple branch models";
//        this.substitutionModel = branchModel.getRootSubstitutionModel();
        this.substitutionModel = substitutionModel;
        this.branchLengths = new double[nodeCount];
    }


    @Override
    public DMatrixSparseCSC getScaledInstantaneousMatrix(int nodeIndex) {
        DMatrixSparseCSC scaledQ = DConvertMatrixStruct.convert(sparseQ, (DMatrixSparseCSC) null);
        CommonOps_DSCC.scale(branchLengths[nodeIndex], scaledQ, scaledQ);
        return scaledQ;
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
        return 1;
    }

    @Override
    public int getInfinitesimalMatrixBufferIndex(int branchIndex) {
        return 1;
    }

    @Override
    public int getInfinitesimalSquaredMatrixBufferIndex(int branchIndex) {
        return 1;
    }

    @Override
    public int getFirstOrderDifferentialMatrixBufferIndex(int branchIndex) {
        return 1;
    }

    @Override
    public int getSecondOrderDifferentialMatrixBufferIndex(int branchIndex) {
        return 1;
    }

    @Override
    public void cacheInfinitesimalMatrix(Beagle beagle, int bufferIndex, double[] differentialMatrix) {

    }

    @Override
    public void cacheInfinitesimalSquaredMatrix(Beagle beagle, int bufferIndex, double[] differentialMatrix) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public void cacheFirstOrderDifferentialMatrix(Beagle beagle, int branchIndex, double[] differentialMassMatrix) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public int getCachedMatrixBufferCount(PreOrderSettings settings) {
        throw new RuntimeException("Not yet implemented.");
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
    public SubstitutionModel getSubstitutionModelForBranch(int branchIndex) {
        return substitutionModel;
    }

    @Override
    public int getEigenIndex(int bufferIndex) {
        throw new RuntimeException("Eigen index should not be called for actions.");
    }

    @Override
    public int getMatrixIndex(int branchIndex) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public double[] getRootStateFrequencies() {
        return substitutionModel.getFrequencyModel().getFrequencies();
    }

    @Override
    public void updateSubstitutionModels(Beagle beagle, boolean flipBuffers) {
        final int stateCount = substitutionModel.getFrequencyModel().getFrequencyCount();
        sparseQ = new DMatrixSparseTriplet(stateCount, stateCount, stateCount * stateCount);
        double[] Q = new double[stateCount * stateCount];
        substitutionModel.getInfinitesimalMatrix(Q);
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                sparseQ.addItem(i, j, Q[i * stateCount + j]);
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
}
