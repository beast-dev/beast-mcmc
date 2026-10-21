/*
 * HierarchicalSubstitutionEvolutionaryProcess.java
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

package dr.evomodel.treedatalikelihood;

import beagle.Beagle;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;

/**
 * @author Jeff Thorne
 * @author Tae-Kun Seo
 * @author Sangchul Choi
 * @author Xiang Ji
 */
public final class HierarchicalSubstitutionEvolutionaryProcess implements EvolutionaryProcessDelegate {

    private final HierarchicalSubstitutionModel substitutionModel;
    private final FrequencyModel rootFrequencyModel;

    // two model slots (BEAGLE eigen buffers), flipped on a model change so a rejected move re-uploads nothing
    private final BufferIndexHelper eigenBufferHelper;
    private final BufferIndexHelper matrixBufferHelper;

    public HierarchicalSubstitutionEvolutionaryProcess(Tree tree, BranchModel branchModel, SiteRateModel siteRateModel) {

        if (branchModel.getSubstitutionModels().size() != 1 ||
                !(branchModel.getRootSubstitutionModel() instanceof HierarchicalSubstitutionModel)) {
            throw new IllegalArgumentException(
                    "The hierarchical representation needs a single hierarchicalSubstitutionModel on every branch");
        }

        this.substitutionModel = (HierarchicalSubstitutionModel) branchModel.getRootSubstitutionModel();
        if (!substitutionModel.isBinary()) {
            throw new IllegalArgumentException("The hierarchical representation needs a binary character tree");
        }

        this.rootFrequencyModel = branchModel.getRootFrequencyModel();

        eigenBufferHelper = new BufferIndexHelper(1, 0);
        matrixBufferHelper = new BufferIndexHelper(tree.getNodeCount(), 0, 0);
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return false;
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
    public int getInfinitesimalMatrixBufferIndex(int branchIndex) {
        throw new UnsupportedOperationException("Gradients are not yet available with the hierarchical representation");
    }

    @Override
    public int getInfinitesimalSquaredMatrixBufferIndex(int branchIndex) {
        throw new UnsupportedOperationException("Gradients are not yet available with the hierarchical representation");
    }

    @Override
    public int getFirstOrderDifferentialMatrixBufferIndex(int branchIndex) {
        throw new UnsupportedOperationException("Gradients are not yet available with the hierarchical representation");
    }

    @Override
    public int getSecondOrderDifferentialMatrixBufferIndex(int branchIndex) {
        throw new UnsupportedOperationException("Gradients are not yet available with the hierarchical representation");
    }

    @Override
    public void cacheInfinitesimalMatrix(Beagle beagle, int bufferIndex, double[] differentialMatrix) {
        throw new UnsupportedOperationException("Gradients are not yet available with the hierarchical representation");
    }

    @Override
    public void cacheInfinitesimalSquaredMatrix(Beagle beagle, int bufferIndex, double[] differentialMatrix) {
        throw new UnsupportedOperationException("Gradients are not yet available with the hierarchical representation");
    }

    @Override
    public void cacheFirstOrderDifferentialMatrix(Beagle beagle, int branchIndex, double[] differentialMassMatrix) {
        throw new UnsupportedOperationException("Gradients are not yet available with the hierarchical representation");
    }

    @Override
    public int getCachedMatrixBufferCount(PreOrderSettings settings) {
        return 0;
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
        return eigenBufferHelper.getOffsetIndex(bufferIndex);
    }

    @Override
    public int getMatrixIndex(int branchIndex) {
        return matrixBufferHelper.getOffsetIndex(branchIndex);
    }

    @Override
    public double[] getRootStateFrequencies() {
        return rootFrequencyModel.getFrequencies();
    }

    @Override
    public void updateSubstitutionModels(Beagle beagle, boolean flipBuffers) {
        if (flipBuffers) {
            eigenBufferHelper.flipOffset(0);
        }
        // The model goes into an eigen buffer; see BEAGLE_FLAG_HIERARCHICAL_REPRESENTATION in beagle.h.
        beagle.setEigenDecomposition(eigenBufferHelper.getOffsetIndex(0),
                substitutionModel.getBeagleModelTriples(),
                substitutionModel.getBeagleModelUnused(),
                substitutionModel.getBeagleModelLambdas());
    }

    @Override
    public void updateTransitionMatrices(Beagle beagle, int[] branchIndices, double[] edgeLengths, int updateCount,
                                         boolean flipBuffers) {
        int[] probabilityIndices = new int[updateCount];
        for (int i = 0; i < updateCount; i++) {
            if (flipBuffers) {
                matrixBufferHelper.flipOffset(branchIndices[i]);
            }
            probabilityIndices[i] = matrixBufferHelper.getOffsetIndex(branchIndices[i]);
        }
        // One call for all flagged branches: BEAGLE builds the bundles from the model, the category rates
        // and state frequencies buffer 0 it already holds.
        beagle.updateTransitionMatrices(eigenBufferHelper.getOffsetIndex(0), probabilityIndices,
                null, null, edgeLengths, updateCount);
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
}
