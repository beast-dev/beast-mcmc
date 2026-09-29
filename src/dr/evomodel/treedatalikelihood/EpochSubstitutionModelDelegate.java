/*
 * EpochSubstitutionModelDelegate.java
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
import dr.evomodel.branchratemodel.TimeVaryingBranchRateModel.EpochTimeProvider;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;

import java.util.List;

/**
 * Handles one substitution model for each epoch. The branches must be augmented with degree-2 nodes at the epoch
 * transition times (see EpochLikelihoodTraversal), so that the matrix of each branch segment comes from a single
 * eigen decomposition and no convolution of matrices is needed.
 * <p>
 * The matrix buffers cover the original nodes and the augmented nodes. The epoch of the matrix of a node is looked
 * up in the AugmentedNodeRegistry, which the traversal fills.
 * <p>
 * One further matrix buffer provides an identity matrix. A branch length of zero makes every eigenvalue term
 * exp(0 t) = 1, and so gives the identity, with any eigen decomposition, in the standard and in the spectral
 * representation.
 *
 * @author Marc A Suchard
 */
public class EpochSubstitutionModelDelegate implements EpochEvolutionaryProcessDelegate {

    private final List<SubstitutionModel> substitutionModelList;
    private final FrequencyModel rootFrequencyModel;
    private final EpochTimeProvider epochTimeProvider;

    private final int eigenCount;
    private final double[] transitionTimes;

    private final AugmentedNodeRegistry registry;

    private final BufferIndexHelper eigenBufferHelper;
    private final BufferIndexHelper matrixBufferHelper;

    private final int identityMatrixIndex;
    private boolean identityMatrixSet = false;

    /**
     * @param tree              the tree
     * @param branchModel       supplies the substitution models, ordered from the most recent epoch, and the root
     *                          frequencies
     * @param epochTimeProvider supplies the start time of each epoch: the first entry is the start of the most
     *                          recent epoch and the remaining entries are the transition times
     */
    public EpochSubstitutionModelDelegate(Tree tree, BranchModel branchModel, EpochTimeProvider epochTimeProvider) {
        this(tree, branchModel, epochTimeProvider, 0);
    }

    /**
     * @param augmentedNodeCapacity  maximum number of augmented nodes in use at once; if not positive then
     *                               taxon count times number of transition times, which cannot be exceeded
     */
    public EpochSubstitutionModelDelegate(Tree tree, BranchModel branchModel, EpochTimeProvider epochTimeProvider,
                                          int augmentedNodeCapacity) {

        this.substitutionModelList = branchModel.getSubstitutionModels();
        this.rootFrequencyModel = branchModel.getRootFrequencyModel();
        this.epochTimeProvider = epochTimeProvider;

        eigenCount = substitutionModelList.size();

        final int boundaryCount = eigenCount - 1;
        if (epochTimeProvider.getEpochTimes().length != eigenCount) {
            throw new IllegalArgumentException("EpochSubstitutionModelDelegate needs one epoch start time for each of " +
                    eigenCount + " substitution models but found " + epochTimeProvider.getEpochTimes().length);
        }
        transitionTimes = new double[boundaryCount];

        final int capacity = augmentedNodeCapacity > 0 ? augmentedNodeCapacity :
                tree.getExternalNodeCount() * boundaryCount;

        registry = new AugmentedNodeRegistry(tree.getNodeCount(), boundaryCount, capacity);

        // two eigen buffers for each decomposition for store and restore
        eigenBufferHelper = new BufferIndexHelper(eigenCount, 0);

        // two matrices for each original and augmented node for store and restore, then the identity
        matrixBufferHelper = new BufferIndexHelper(registry.getTotalNodeCount(), 0);
        identityMatrixIndex = matrixBufferHelper.getBufferCount();
    }

    @Override
    public double[] getEpochTransitionTimes() {
        final double[] times = epochTimeProvider.getEpochTimes();

        System.arraycopy(times, 1, transitionTimes, 0, transitionTimes.length);
        for (int i = 1; i < transitionTimes.length; ++i) {
            if (!(transitionTimes[i - 1] < transitionTimes[i])) {
                throw new IllegalArgumentException("Epoch transition times must be strictly ascending");
            }
        }

        return transitionTimes;
    }

    @Override
    public AugmentedNodeRegistry getAugmentedNodeRegistry() {
        return registry;
    }

    @Override
    public int getIdentityMatrixIndex() {
        return identityMatrixIndex;
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        for (SubstitutionModel model : substitutionModelList) {
            if (model.canReturnComplexDiagonalization()) {
                return true;
            }
        }
        return false;
    }

    @Override
    public int getEigenBufferCount() {
        return eigenBufferHelper.getBufferCount();
    }

    @Override
    public int getMatrixBufferCount() {
        return identityMatrixIndex + 1;
    }

    /**
     * The infinitesimal matrix of an epoch is cached once, and used by every segment in the epoch.
     *
     * @param branchIndex an original or augmented node
     * @return the buffer of the infinitesimal matrix of the substitution model of the segment above the node
     */
    @Override
    public int getInfinitesimalMatrixBufferIndex(int branchIndex) {
        return getInfinitesimalMatrixBufferIndexForEpoch(registry.getMatrixEpoch(branchIndex));
    }

    @Override
    public int getInfinitesimalSquaredMatrixBufferIndex(int branchIndex) {
        return getInfinitesimalSquaredMatrixBufferIndexForEpoch(registry.getMatrixEpoch(branchIndex));
    }

    private int getInfinitesimalMatrixBufferIndexForEpoch(int epoch) {
        return getMatrixBufferCount() + eigenBufferHelper.getOffsetIndex(epoch);
    }

    private int getInfinitesimalSquaredMatrixBufferIndexForEpoch(int epoch) {
        return getMatrixBufferCount() + getEigenBufferCount() + eigenBufferHelper.getOffsetIndex(epoch);
    }

    @Override
    public int getFirstOrderDifferentialMatrixBufferIndex(int branchIndex) {
        throw new UnsupportedOperationException(DIFFERENTIAL_MASS_MESSAGE);
    }

    @Override
    public int getSecondOrderDifferentialMatrixBufferIndex(int branchIndex) {
        throw new UnsupportedOperationException(DIFFERENTIAL_MASS_MESSAGE);
    }

    /**
     * @param bufferIndex the epoch, that is, the index of the substitution model
     */
    @Override
    public void cacheInfinitesimalMatrix(Beagle beagle, int bufferIndex, double[] differentialMatrix) {
        beagle.setDifferentialMatrix(getInfinitesimalMatrixBufferIndexForEpoch(bufferIndex), differentialMatrix);
    }

    /**
     * @param bufferIndex the epoch, that is, the index of the substitution model
     */
    @Override
    public void cacheInfinitesimalSquaredMatrix(Beagle beagle, int bufferIndex, double[] differentialMatrix) {
        beagle.setDifferentialMatrix(getInfinitesimalSquaredMatrixBufferIndexForEpoch(bufferIndex),
                differentialMatrix);
    }

    @Override
    public void cacheFirstOrderDifferentialMatrix(Beagle beagle, int branchIndex, double[] differentialMassMatrix) {
        throw new UnsupportedOperationException(DIFFERENTIAL_MASS_MESSAGE);
    }

    @Override
    public int getCachedMatrixBufferCount(PreOrderSettings settings) {
        if (settings.branchInfinitesimalDerivative) {
            throw new UnsupportedOperationException(DIFFERENTIAL_MASS_MESSAGE);
        }
        return settings.branchRateDerivative ? 2 * getEigenBufferCount() : 0;
    }

    @Override
    public int getSubstitutionModelCount() {
        return substitutionModelList.size();
    }

    @Override
    public SubstitutionModel getSubstitutionModel(int index) {
        return substitutionModelList.get(index);
    }

    /**
     * @param branchIndex an original or augmented node
     * @return the substitution model of the segment above the node
     */
    @Override
    public SubstitutionModel getSubstitutionModelForBranch(int branchIndex) {
        return getSubstitutionModel(registry.getMatrixEpoch(branchIndex));
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
        for (int i = 0; i < eigenCount; i++) {
            if (flipBuffers) {
                eigenBufferHelper.flipOffset(i);
            }

            EigenDecomposition ed = substitutionModelList.get(i).getEigenDecomposition();

            beagle.setEigenDecomposition(
                    eigenBufferHelper.getOffsetIndex(i),
                    ed.getEigenVectors(),
                    ed.getInverseEigenVectors(),
                    ed.getEigenValues());
        }
    }

    /**
     * The identity does not depend on the category rates, but BEAGLE needs them to have been set before it updates
     * any matrix, which is after updateSubstitutionModels() has been called for the first time.
     */
    private void setIdentityMatrix(Beagle beagle) {
        beagle.updateTransitionMatrices(eigenBufferHelper.getOffsetIndex(0),
                new int[]{identityMatrixIndex},
                null, // firstDerivativeIndices
                null, // secondDerivativeIndices
                new double[]{0.0},
                1);

        identityMatrixSet = true;
    }

    /**
     * @param branchIndices original and augmented nodes whose matrices to update
     */
    @Override
    public void updateTransitionMatrices(Beagle beagle, int[] branchIndices, double[] edgeLengths, int updateCount,
                                         boolean flipBuffers) {

        if (!identityMatrixSet) {
            setIdentityMatrix(beagle);
        }

        final int[] counts = new int[eigenCount];
        for (int i = 0; i < updateCount; i++) {
            ++counts[registry.getMatrixEpoch(branchIndices[i])];
        }

        final int[][] probabilityIndices = new int[eigenCount][];
        final double[][] lengths = new double[eigenCount][];
        for (int k = 0; k < eigenCount; k++) {
            if (counts[k] > 0) {
                probabilityIndices[k] = new int[counts[k]];
                lengths[k] = new double[counts[k]];
                counts[k] = 0;
            }
        }

        for (int i = 0; i < updateCount; i++) {
            final int branch = branchIndices[i];
            final int k = registry.getMatrixEpoch(branch);

            if (flipBuffers) {
                matrixBufferHelper.flipOffset(branch);
            }
            probabilityIndices[k][counts[k]] = matrixBufferHelper.getOffsetIndex(branch);
            lengths[k][counts[k]] = edgeLengths[i];
            ++counts[k];
        }

        for (int k = 0; k < eigenCount; k++) {
            if (counts[k] > 0) {
                beagle.updateTransitionMatrices(eigenBufferHelper.getOffsetIndex(k),
                        probabilityIndices[k],
                        null, // firstDerivativeIndices
                        null, // secondDerivativeIndices
                        lengths[k],
                        counts[k]);
            }
        }
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
        registry.storeState();
    }

    @Override
    public void restoreState() {
        eigenBufferHelper.restoreState();
        matrixBufferHelper.restoreState();
        registry.restoreState();
    }

    private static final String DIFFERENTIAL_MASS_MESSAGE =
            "Branch-specific differential mass matrices are not yet supported for epoch models with augmented nodes";
}
