/*
 * SequentialCanonicalOUMessagePasser.java
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

package dr.evomodel.treedatalikelihood.continuous.integration;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.adapter.HomogeneousCanonicalOUBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalRootPrior;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTreeMessagePasser;
import dr.evomodel.treedatalikelihood.continuous.framework.MatrixUtils;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.timeseries.engine.gaussian.CanonicalBranchMessageContribution;
import dr.inference.timeseries.engine.gaussian.CanonicalBranchMessageContributionUtils;
import dr.inference.timeseries.engine.gaussian.CanonicalGaussianMessageOps;
import dr.inference.timeseries.engine.gaussian.CanonicalLocalTransitionAdjoints;
import dr.inference.timeseries.engine.gaussian.CanonicalTransitionAdjointUtils;
import dr.inference.timeseries.gaussian.OrthogonalBlockDiagonalSelectionMatrixParameterization;
import dr.inference.timeseries.gaussian.OUProcessModel;
import dr.inference.timeseries.representation.CanonicalGaussianState;
import dr.inference.timeseries.representation.CanonicalGaussianTransition;
import dr.inference.timeseries.representation.CanonicalGaussianUtils;
import dr.util.TaskPool;

import java.util.Arrays;

/**
 * <p>This class implements the canonical OU tree traversals and exact/partially
 * observed tip elimination in canonical form.
 */
public final class SequentialCanonicalOUMessagePasser implements CanonicalTreeMessagePasser {

    private static final double BRANCH_LENGTH_FD_RELATIVE_STEP = 1.0e-6;
    private static final double BRANCH_LENGTH_FD_ABSOLUTE_STEP = 1.0e-8;
    private static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);
    private static final String BRANCH_GRADIENT_THREADS_PROPERTY =
            "beast.experimental.canonicalBranchGradientThreads";

    private final Tree tree;
    private final int dim;
    private final int nodeCount;
    private final int tipCount;

    private CanonicalTipObservation[] tipObservations;
    private CanonicalTipObservation[] storedTipObservations;
    private CanonicalGaussianState[] postOrder;
    private CanonicalGaussianState[] storedPostOrder;
    private CanonicalGaussianState[] preOrder;
    private CanonicalGaussianState[] storedPreOrder;
    private CanonicalGaussianState[] branchAboveParent;
    private CanonicalGaussianState[] storedBranchAboveParent;

    private boolean hasPostOrderState;
    private boolean storedHasPostOrderState;
    private boolean hasPreOrderState;
    private boolean storedHasPreOrderState;
    private double lastRootDiffusionScale;
    private double storedLastRootDiffusionScale;
    private boolean hasFixedRootValue;
    private boolean storedHasFixedRootValue;

    private double[] fixedRootValue;
    private double[] storedFixedRootValue;
    private final BranchGradientInputs preparedBranchGradientInputs;
    private final BranchGradientWorkspace mainWorkspace;
    private final TaskPool branchGradientTaskPool;
    private final BranchGradientWorkspace[] branchGradientWorkspaces;

    /**
     * Frozen branch-local adjoints prepared from the current post-order and
     * pre-order passes. These inputs are tree-blind: each active entry is just
     * one independent branch contribution plus the branch length needed by the
     * OU pullbacks.
     */
    public static final class BranchGradientInputs {
        private final int dimension;
        private final int capacity;
        private final int nodeSlotCount;
        private final int[] activeChildIndices;
        private final double[] activeBranchLengths;
        private final CanonicalLocalTransitionAdjoints[] activeAdjoints;
        private final OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis[] activeOrthogonalPreparedBasis;
        private final boolean[] stagedActiveByChild;
        private final double[] stagedBranchLengthsByChild;
        private final CanonicalLocalTransitionAdjoints[] stagedAdjointsByChild;
        private final OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis[] stagedOrthogonalPreparedBasisByChild;
        private final CanonicalGaussianState rootPreOrder;
        private final CanonicalGaussianState rootPostOrder;
        private int activeBranchCount;
        private double rootDiffusionScale;
        private boolean hasOrthogonalPreparedBasis;

        private BranchGradientInputs(final int capacity, final int dimension) {
            this.dimension = dimension;
            this.capacity = capacity;
            this.nodeSlotCount = capacity + 1;
            this.activeChildIndices = new int[capacity];
            this.activeBranchLengths = new double[capacity];
            this.activeAdjoints = new CanonicalLocalTransitionAdjoints[capacity];
            this.activeOrthogonalPreparedBasis =
                    new OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis[capacity];
            this.stagedActiveByChild = new boolean[nodeSlotCount];
            this.stagedBranchLengthsByChild = new double[nodeSlotCount];
            this.stagedAdjointsByChild = new CanonicalLocalTransitionAdjoints[nodeSlotCount];
            this.stagedOrthogonalPreparedBasisByChild =
                    new OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis[nodeSlotCount];
            for (int i = 0; i < capacity; ++i) {
                this.activeAdjoints[i] = new CanonicalLocalTransitionAdjoints(dimension);
            }
            for (int i = 0; i < nodeSlotCount; ++i) {
                this.stagedAdjointsByChild[i] = new CanonicalLocalTransitionAdjoints(dimension);
            }
            this.rootPreOrder = new CanonicalGaussianState(dimension);
            this.rootPostOrder = new CanonicalGaussianState(dimension);
            this.activeBranchCount = 0;
            this.rootDiffusionScale = 0.0;
            this.hasOrthogonalPreparedBasis = false;
        }

        public int getDimension() {
            return dimension;
        }

        public int getActiveBranchCount() {
            return activeBranchCount;
        }

        public int getActiveChildIndex(final int activeIndex) {
            checkActiveIndex(activeIndex);
            return activeChildIndices[activeIndex];
        }

        public double getBranchLength(final int activeIndex) {
            checkActiveIndex(activeIndex);
            return activeBranchLengths[activeIndex];
        }

        public CanonicalLocalTransitionAdjoints getLocalAdjoints(final int activeIndex) {
            checkActiveIndex(activeIndex);
            return activeAdjoints[activeIndex];
        }

        public double getRootDiffusionScale() {
            return rootDiffusionScale;
        }

        public CanonicalGaussianState getRootPreOrderState() {
            return rootPreOrder;
        }

        public CanonicalGaussianState getRootPostOrderState() {
            return rootPostOrder;
        }

        private void clear() {
            Arrays.fill(stagedActiveByChild, false);
            activeBranchCount = 0;
            rootDiffusionScale = 0.0;
            hasOrthogonalPreparedBasis = false;
        }

        private void checkCompatible(final int expectedCapacity,
                                     final int expectedDimension) {
            if (capacity < expectedCapacity) {
                throw new IllegalArgumentException(
                        "BranchGradientInputs capacity " + capacity
                                + " is too small for " + expectedCapacity + " branches.");
            }
            if (dimension != expectedDimension) {
                throw new IllegalArgumentException(
                        "BranchGradientInputs dimension mismatch: "
                                + dimension + " vs " + expectedDimension + ".");
            }
        }

        private void addBranch(final int childIndex,
                               final double branchLength,
                               final CanonicalLocalTransitionAdjoints source,
                               final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalSelection,
                               final double[] stationaryMean) {
            if (activeBranchCount >= capacity) {
                throw new IllegalStateException("BranchGradientInputs capacity exceeded.");
            }
            if (orthogonalSelection != null) {
                hasOrthogonalPreparedBasis = true;
                OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis prepared =
                        activeOrthogonalPreparedBasis[activeBranchCount];
                if (prepared == null) {
                    prepared = orthogonalSelection.createPreparedBranchBasis();
                    activeOrthogonalPreparedBasis[activeBranchCount] = prepared;
                }
                orthogonalSelection.prepareBranchBasis(branchLength, stationaryMean, prepared);
            }
            activeChildIndices[activeBranchCount] = childIndex;
            activeBranchLengths[activeBranchCount] = branchLength;
            copyAdjoints(source, activeAdjoints[activeBranchCount]);
            activeBranchCount++;
        }

        private void stageBranch(final int childIndex,
                                 final double branchLength,
                                 final CanonicalLocalTransitionAdjoints source,
                                 final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalSelection,
                                 final double[] stationaryMean) {
            checkChildIndex(childIndex);
            stagedBranchLengthsByChild[childIndex] = branchLength;
            copyAdjoints(source, stagedAdjointsByChild[childIndex]);
            if (orthogonalSelection != null) {
                OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis prepared =
                        stagedOrthogonalPreparedBasisByChild[childIndex];
                if (prepared == null) {
                    prepared = orthogonalSelection.createPreparedBranchBasis();
                    stagedOrthogonalPreparedBasisByChild[childIndex] = prepared;
                }
                orthogonalSelection.prepareBranchBasis(branchLength, stationaryMean, prepared);
            }
            stagedActiveByChild[childIndex] = true;
        }

        private void clearStagedBranch(final int childIndex) {
            checkChildIndex(childIndex);
            stagedActiveByChild[childIndex] = false;
        }

        private void compactStagedBranches(final int rootIndex,
                                           final boolean useOrthogonalPreparedBasis) {
            checkChildIndex(rootIndex);
            activeBranchCount = 0;
            for (int childIndex = 0; childIndex < nodeSlotCount; ++childIndex) {
                if (childIndex == rootIndex || !stagedActiveByChild[childIndex]) {
                    continue;
                }
                if (activeBranchCount >= capacity) {
                    throw new IllegalStateException("BranchGradientInputs capacity exceeded during compaction.");
                }
                activeChildIndices[activeBranchCount] = childIndex;
                activeBranchLengths[activeBranchCount] = stagedBranchLengthsByChild[childIndex];
                activeAdjoints[activeBranchCount] = stagedAdjointsByChild[childIndex];
                if (useOrthogonalPreparedBasis) {
                    activeOrthogonalPreparedBasis[activeBranchCount] =
                            stagedOrthogonalPreparedBasisByChild[childIndex];
                }
                activeBranchCount++;
            }
            hasOrthogonalPreparedBasis = useOrthogonalPreparedBasis && activeBranchCount > 0;
        }

        private void prepareOrthogonalBasisForActiveBranches(
                final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalSelection,
                final double[] stationaryMean) {
            if (orthogonalSelection == null) {
                hasOrthogonalPreparedBasis = false;
                return;
            }

            for (int activeIndex = 0; activeIndex < activeBranchCount; ++activeIndex) {
                OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis prepared =
                        activeOrthogonalPreparedBasis[activeIndex];
                if (prepared == null) {
                    prepared = orthogonalSelection.createPreparedBranchBasis();
                    activeOrthogonalPreparedBasis[activeIndex] = prepared;
                }
                orthogonalSelection.prepareBranchBasis(activeBranchLengths[activeIndex], stationaryMean, prepared);
            }
            hasOrthogonalPreparedBasis = activeBranchCount > 0;
        }

        private OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis
        getOrthogonalPreparedBasis(final int activeIndex) {
            checkActiveIndex(activeIndex);
            if (!hasOrthogonalPreparedBasis) {
                throw new IllegalStateException("No orthogonal prepared branch basis is available.");
            }
            return activeOrthogonalPreparedBasis[activeIndex];
        }

        private void checkActiveIndex(final int activeIndex) {
            if (activeIndex < 0 || activeIndex >= activeBranchCount) {
                throw new IndexOutOfBoundsException(
                        "Active branch index " + activeIndex
                                + " is outside [0, " + activeBranchCount + ").");
            }
        }

        private void checkChildIndex(final int childIndex) {
            if (childIndex < 0 || childIndex >= nodeSlotCount) {
                throw new IndexOutOfBoundsException(
                        "Child index " + childIndex
                                + " is outside [0, " + nodeSlotCount + ").");
            }
        }
    }

    private static final class BranchGradientWorkspace {
        final double[][] traitCovariance;
        final CanonicalGaussianTransition transition;
        final CanonicalGaussianState state;
        final CanonicalGaussianState siblingProduct;
        final CanonicalGaussianState downwardParentState;
        final CanonicalGaussianState combinedState;
        final CanonicalGaussianState parentPosterior;
        final CanonicalGaussianState pairState;
        final CanonicalGaussianMessageOps.Workspace gaussianWorkspace;
        final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
        final CanonicalTransitionAdjointUtils.Workspace transitionAdjointWorkspace;
        final CanonicalBranchMessageContribution contribution;
        final CanonicalLocalTransitionAdjoints adjoints;
        final int[] observedIndexScratch;
        final int[] missingIndexScratch;
        final int[] reducedIndexByTraitScratch;
        final double[] mean;
        final double[] mean2;
        final double[] orthogonalStationaryMeanScratch;
        final double[] orthogonalCompressedGradientScratch;
        final double[] orthogonalNativeGradientScratch;
        final double[][] orthogonalRotationGradientScratch;
        final double[][] covariance;
        final double[][] covariance2;
        final double[][] transitionMatrix;
        final double[][] covarianceAdjoint;
        final double[] varianceFlat;
        final double[] precisionFlat;
        final double[] reducedPrecisionFlatScratch;
        final double[] reducedCovarianceFlatScratch;
        final double[] reducedInformationScratch;
        final double[] reducedMeanScratch;
        final double[] localGradientA;
        final double[] localGradientQ;
        final double[] localGradientMuVector;
        final double[] localGradientMuScalar;
        OrthogonalBlockDiagonalSelectionMatrixParameterization.BranchGradientWorkspace orthogonalBranchWorkspace;

        BranchGradientWorkspace(final int dim) {
            this.traitCovariance = new double[dim][dim];
            this.transition = new CanonicalGaussianTransition(dim);
            this.state = new CanonicalGaussianState(dim);
            this.siblingProduct = new CanonicalGaussianState(dim);
            this.downwardParentState = new CanonicalGaussianState(dim);
            this.combinedState = new CanonicalGaussianState(dim);
            this.parentPosterior = new CanonicalGaussianState(dim);
            this.pairState = new CanonicalGaussianState(2 * dim);
            this.gaussianWorkspace = new CanonicalGaussianMessageOps.Workspace(dim);
            this.contributionWorkspace = new CanonicalBranchMessageContributionUtils.Workspace(dim);
            this.transitionAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(dim);
            this.contribution = new CanonicalBranchMessageContribution(dim);
            this.adjoints = new CanonicalLocalTransitionAdjoints(dim);
            this.observedIndexScratch = new int[dim];
            this.missingIndexScratch = new int[dim];
            this.reducedIndexByTraitScratch = new int[dim];
            this.mean = new double[dim];
            this.mean2 = new double[dim];
            this.orthogonalStationaryMeanScratch = new double[dim];
            this.orthogonalCompressedGradientScratch = new double[dim + 2 * (dim / 2)];
            this.orthogonalNativeGradientScratch = new double[((dim & 1) == 1 ? 1 : 0) + 3 * (dim / 2)];
            this.orthogonalRotationGradientScratch = new double[dim][dim];
            this.covariance = new double[dim][dim];
            this.covariance2 = new double[dim][dim];
            this.transitionMatrix = new double[dim][dim];
            this.covarianceAdjoint = new double[dim][dim];
            this.varianceFlat = new double[dim * dim];
            this.precisionFlat = new double[dim * dim];
            this.reducedPrecisionFlatScratch = new double[4 * dim * dim];
            this.reducedCovarianceFlatScratch = new double[4 * dim * dim];
            this.reducedInformationScratch = new double[2 * dim];
            this.reducedMeanScratch = new double[2 * dim];
            this.localGradientA = new double[dim * dim];
            this.localGradientQ = new double[dim * dim];
            this.localGradientMuVector = new double[dim];
            this.localGradientMuScalar = new double[1];
        }

        double[] localGradientMu(final int gradientLength, final int dim) {
            if (gradientLength == 1) {
                return localGradientMuScalar;
            }
            if (gradientLength == dim) {
                return localGradientMuVector;
            }
            throw new IllegalArgumentException(
                    "Stationary-mean gradient length must be 1 or " + dim + ", found " + gradientLength);
        }

        void clearLocalGradientBuffers(final int gradALength,
                                       final int gradMuLength,
                                       final int dim,
                                       final boolean orthogonalSelection,
                                       final int compressedGradientLength) {
            Arrays.fill(localGradientA, 0, gradALength, 0.0);
            Arrays.fill(localGradientQ, 0.0);
            Arrays.fill(localGradientMu(gradMuLength, dim), 0.0);
            if (orthogonalSelection) {
                Arrays.fill(orthogonalCompressedGradientScratch, 0, compressedGradientLength, 0.0);
                for (int i = 0; i < orthogonalRotationGradientScratch.length; ++i) {
                    Arrays.fill(orthogonalRotationGradientScratch[i], 0.0);
                }
            }
        }

        OrthogonalBlockDiagonalSelectionMatrixParameterization.BranchGradientWorkspace
        ensureOrthogonalBranchWorkspace(final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalSelection) {
            if (orthogonalBranchWorkspace == null) {
                orthogonalBranchWorkspace = orthogonalSelection.createBranchGradientWorkspace();
            }
            return orthogonalBranchWorkspace;
        }
    }

    public SequentialCanonicalOUMessagePasser(final Tree tree, final int dim) {
        this(tree, dim, resolveBranchGradientParallelism());
    }

    public SequentialCanonicalOUMessagePasser(final Tree tree,
                                              final int dim,
                                              final int branchGradientParallelism) {
        this.tree = tree;
        this.dim = dim;
        this.nodeCount = tree.getNodeCount();
        this.tipCount = tree.getExternalNodeCount();
        this.tipObservations = allocateTipObservations(nodeCount, dim);
        this.storedTipObservations = allocateTipObservations(nodeCount, dim);
        this.postOrder = allocateStates(nodeCount, dim);
        this.storedPostOrder = allocateStates(nodeCount, dim);
        this.preOrder = allocateStates(nodeCount, dim);
        this.storedPreOrder = allocateStates(nodeCount, dim);
        this.branchAboveParent = allocateStates(nodeCount, dim);
        this.storedBranchAboveParent = allocateStates(nodeCount, dim);
        this.preparedBranchGradientInputs = new BranchGradientInputs(Math.max(0, nodeCount - 1), dim);
        this.mainWorkspace = new BranchGradientWorkspace(dim);
        this.fixedRootValue = new double[dim];
        this.storedFixedRootValue = new double[dim];
        this.branchGradientTaskPool = new TaskPool(nodeCount, Math.max(1, branchGradientParallelism));
        final int branchGradientWorkspaceCount =
                branchGradientTaskPool.getNumThreads() <= 1 ? 1 : branchGradientTaskPool.getNumThreads() + 1;
        this.branchGradientWorkspaces = new BranchGradientWorkspace[branchGradientWorkspaceCount];
        for (int i = 0; i < branchGradientWorkspaces.length; ++i) {
            branchGradientWorkspaces[i] = new BranchGradientWorkspace(dim);
        }
        this.hasPostOrderState = false;
        this.storedHasPostOrderState = false;
        this.hasPreOrderState = false;
        this.storedHasPreOrderState = false;
        this.lastRootDiffusionScale = 0.0;
        this.storedLastRootDiffusionScale = 0.0;
        this.hasFixedRootValue = false;
        this.storedHasFixedRootValue = false;
    }

    private static int resolveBranchGradientParallelism() {
        final String propertyValue = System.getProperty(BRANCH_GRADIENT_THREADS_PROPERTY);
        if (propertyValue != null) {
            try {
                return Math.max(1, Integer.parseInt(propertyValue));
            } catch (NumberFormatException ignored) {
                // Fall back to the default hardware-concurrency choice below.
            }
        }
        return Math.max(1, Runtime.getRuntime().availableProcessors());
    }

    private int branchGradientChunkSize(final int taskLimit,
                                        final int targetChunksPerWorker,
                                        final int maxChunkSize) {
        final int workerCount = Math.max(1, branchGradientWorkspaces.length);
        final int suggested =
                (taskLimit + workerCount * targetChunksPerWorker - 1) / (workerCount * targetChunksPerWorker);
        return Math.max(1, Math.min(maxChunkSize, suggested));
    }

    private int branchGradientPreparationChunkSize(final int taskLimit) {
        return branchGradientChunkSize(taskLimit, 2, 16);
    }

    private int branchGradientJointChunkSize(final int taskLimit) {
        return branchGradientChunkSize(taskLimit, 3, 8);
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public int getTipCount() {
        return tipCount;
    }

    @Override
    public void setTipObservation(final int tipIndex, final CanonicalTipObservation observation) {
        tipObservations[tipIndex].copyFrom(observation);
    }

    @Override
    public double computePostOrderLogLikelihood(final CanonicalBranchTransitionProvider transitionProvider,
                                                final CanonicalRootPrior rootPrior) {
        if (tree.isExternal(tree.getRoot())) {
            throw new UnsupportedOperationException("Single-node trees are not yet supported by the canonical OU passer.");
        }

        for (int i = 0; i < nodeCount; i++) {
            if (!tree.isExternal(tree.getNode(i))) {
                computePostOrderAtInternalNode(i, transitionProvider, mainWorkspace);
            }
        }

        final int rootIndex = tree.getRoot().getNumber();
        transitionProvider.fillTraitCovariance(mainWorkspace.traitCovariance);
        hasPostOrderState = true;
        return rootPrior.computeLogMarginalLikelihood(postOrder[rootIndex], mainWorkspace.traitCovariance);
    }

    @Override
    public CanonicalGaussianState getPostOrderState(final int nodeIndex) {
        return postOrder[nodeIndex];
    }

    @Override
    public void computePreOrder(final CanonicalBranchTransitionProvider transitionProvider,
                                final CanonicalRootPrior rootPrior) {
        final int rootIndex = tree.getRoot().getNumber();
        transitionProvider.fillTraitCovariance(mainWorkspace.traitCovariance);
        CanonicalGaussianMessageOps.clearState(preOrder[rootIndex]);
        hasFixedRootValue = rootPrior.isFixedRoot();
        if (hasFixedRootValue) {
            rootPrior.fillFixedRootValue(fixedRootValue);
            lastRootDiffusionScale = 0.0;
        } else {
            rootPrior.fillRootPriorState(mainWorkspace.traitCovariance, preOrder[rootIndex]);
            lastRootDiffusionScale = rootPrior.getDiffusionScale();
        }
        clearState(branchAboveParent[rootIndex]);
        computePreOrderRecursive(tree.getRoot().getNumber(), transitionProvider, mainWorkspace);
        hasPreOrderState = true;
    }

    @Override
    public CanonicalGaussianState getPreOrderState(final int nodeIndex) {
        return preOrder[nodeIndex];
    }

    @Override
    public void computeGradientQ(final CanonicalBranchTransitionProvider transitionProvider, final double[] gradQ) {
        Arrays.fill(gradQ, 0.0);
        ensureGradientState();

        final HomogeneousCanonicalOUBranchTransitionProvider ouProvider = requireOUProvider(transitionProvider);
        final OUProcessModel processModel = ouProvider.getProcessModel();
        final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalSelection =
                processModel.getSelectionMatrixParameterization()
                        instanceof OrthogonalBlockDiagonalSelectionMatrixParameterization
                        ? (OrthogonalBlockDiagonalSelectionMatrixParameterization)
                        processModel.getSelectionMatrixParameterization()
                        : null;
        final int rootIndex = tree.getRoot().getNumber();
        final BranchGradientWorkspace workspace = mainWorkspace;

        for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
            if (childIndex == rootIndex) {
                continue;
            }
            if (!fillLocalAdjointsForBranch(childIndex, transitionProvider, workspace)) {
                continue;
            }
            if (orthogonalSelection != null) {
                transposeFromFlatInto(workspace.adjoints.dLogL_dOmega, workspace.covarianceAdjoint, dim);
                orthogonalSelection.accumulateDiffusionGradient(
                        processModel.getDiffusionMatrix(),
                        transitionProvider.getEffectiveBranchLength(childIndex),
                        workspace.covarianceAdjoint,
                        gradQ);
            } else {
                copyFlatToSquare(workspace.adjoints.dLogL_dOmega, workspace.covarianceAdjoint, dim);
                processModel.accumulateDiffusionGradient(
                        transitionProvider.getEffectiveBranchLength(childIndex),
                        workspace.covarianceAdjoint,
                        gradQ);
            }
        }

        if (lastRootDiffusionScale > 0.0) {
            accumulateRootDiffusionGradient(gradQ, workspace);
        }
    }

    @Override
    public void computeGradientBranchLengths(final CanonicalBranchTransitionProvider transitionProvider,
                                             final double[] gradT) {
        Arrays.fill(gradT, 0.0);
        ensureGradientState();

        final HomogeneousCanonicalOUBranchTransitionProvider ouProvider = requireOUProvider(transitionProvider);
        final int rootIndex = tree.getRoot().getNumber();
        for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
            if (childIndex == rootIndex) {
                continue;
            }
            final CanonicalTipObservation tipObservation =
                    tree.isExternal(tree.getNode(childIndex)) ? tipObservations[childIndex] : null;
            if (tipObservation != null && tipObservation.observedCount == 0) {
                gradT[childIndex] = 0.0;
                continue;
            }
            final double branchLength = transitionProvider.getEffectiveBranchLength(childIndex);
            gradT[childIndex] = finiteDifferenceBranchLengthGradient(childIndex, ouProvider, branchLength);
        }
    }

    @Override
    public void computeGradientA(final CanonicalBranchTransitionProvider transitionProvider, final double[] gradA) {
        Arrays.fill(gradA, 0.0);
        ensureGradientState();

        final HomogeneousCanonicalOUBranchTransitionProvider ouProvider = requireOUProvider(transitionProvider);
        final OUProcessModel processModel = ouProvider.getProcessModel();
        if (processModel.getSelectionMatrixParameterization()
                instanceof OrthogonalBlockDiagonalSelectionMatrixParameterization) {
            computeOrthogonalBlockGradientA(
                    transitionProvider,
                    processModel,
                    (OrthogonalBlockDiagonalSelectionMatrixParameterization)
                            processModel.getSelectionMatrixParameterization(),
                    gradA);
            return;
        }
        final int rootIndex = tree.getRoot().getNumber();
        final BranchGradientWorkspace workspace = mainWorkspace;

        for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
            if (childIndex == rootIndex) {
                continue;
            }
            if (!fillLocalAdjointsForBranch(childIndex, transitionProvider, workspace)) {
                continue;
            }

            final double branchLength = transitionProvider.getEffectiveBranchLength(childIndex);
            copyFlatToSquare(workspace.adjoints.dLogL_dF, workspace.transitionMatrix, dim);
            processModel.accumulateSelectionGradient(
                    branchLength,
                    workspace.transitionMatrix,
                    workspace.adjoints.dLogL_df,
                    gradA);

            transposeFromFlatInto(workspace.adjoints.dLogL_dOmega, workspace.covarianceAdjoint, dim);
            processModel.accumulateSelectionGradientFromCovariance(
                    branchLength,
                    workspace.covarianceAdjoint,
                    gradA);
        }

        transposeFlatSquareInPlace(gradA, dim);
    }

    private void computeOrthogonalBlockGradientA(final CanonicalBranchTransitionProvider transitionProvider,
                                                 final OUProcessModel processModel,
                                                 final OrthogonalBlockDiagonalSelectionMatrixParameterization parameterization,
                                                 final double[] gradA) {
        final BranchGradientWorkspace workspace = mainWorkspace;
        final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter =
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) parameterization.getMatrixParameter();
        if (!(blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider)) {
            throw new IllegalStateException(
                    "Orthogonal block native gradient requires an OrthogonalMatrixProvider rotation parameter.");
        }

        final int compressedBlockDim = blockParameter.getCompressedDDimension();
        final int nativeBlockDim = blockParameter.getBlockDiagonalNParameters();
        final OrthogonalMatrixProvider orthogonalRotation =
                (OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter();
        final int angleDim = orthogonalRotation.getOrthogonalParameter().getDimension();
        final int nativeDim = nativeBlockDim + angleDim;
        if (gradA.length != nativeDim) {
            throw new IllegalArgumentException(
                    "Orthogonal block selection gradient expects native parameter length "
                            + nativeDim + ", found " + gradA.length);
        }

        final int rootIndex = tree.getRoot().getNumber();
        if (compressedBlockDim > workspace.orthogonalCompressedGradientScratch.length
                || nativeBlockDim > workspace.orthogonalNativeGradientScratch.length) {
            throw new IllegalStateException(
                    "Orthogonal block scratch is too small for native gradient dimensions "
                            + compressedBlockDim + " and " + nativeBlockDim + ".");
        }
        final double[] stationaryMean = workspace.orthogonalStationaryMeanScratch;
        processModel.getInitialMean(stationaryMean);
        final double[] compressedBlockGradient = workspace.orthogonalCompressedGradientScratch;
        Arrays.fill(compressedBlockGradient, 0, compressedBlockDim, 0.0);
        final double[] nativeBlockGradient = workspace.orthogonalNativeGradientScratch;
        Arrays.fill(nativeBlockGradient, 0, nativeBlockDim, 0.0);
        final double[][] rotationGradient = workspace.orthogonalRotationGradientScratch;
        clearSquare(rotationGradient);

        for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
            if (childIndex == rootIndex) {
                continue;
            }
            if (!fillLocalAdjointsForBranch(childIndex, transitionProvider, workspace)) {
                continue;
            }

            parameterization.accumulateNativeGradientFromAdjoints(
                    processModel.getDiffusionMatrix(),
                    stationaryMean,
                    transitionProvider.getEffectiveBranchLength(childIndex),
                    workspace.adjoints,
                    compressedBlockGradient,
                    rotationGradient);
        }

        blockParameter.chainGradient(compressedBlockGradient, nativeBlockGradient);
        final double[] angleGradient = orthogonalRotation.pullBackGradient(rotationGradient);
        System.arraycopy(nativeBlockGradient, 0, gradA, 0, nativeBlockDim);
        System.arraycopy(angleGradient, 0, gradA, nativeBlockDim, angleGradient.length);
    }

    @Override
    public void computeGradientMu(final CanonicalBranchTransitionProvider transitionProvider, final double[] gradMu) {
        Arrays.fill(gradMu, 0.0);
        ensureGradientState();

        final HomogeneousCanonicalOUBranchTransitionProvider ouProvider = requireOUProvider(transitionProvider);
        final OUProcessModel processModel = ouProvider.getProcessModel();
        final int rootIndex = tree.getRoot().getNumber();
        final BranchGradientWorkspace workspace = mainWorkspace;

        for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
            if (childIndex == rootIndex) {
                continue;
            }
            if (!fillLocalAdjointsForBranch(childIndex, transitionProvider, workspace)) {
                continue;
            }

            final double branchLength = transitionProvider.getEffectiveBranchLength(childIndex);
            processModel.fillTransitionMatrix(branchLength, workspace.transitionMatrix);
            accumulateStationaryMeanGradient(
                    workspace.transitionMatrix,
                    workspace.adjoints.dLogL_df,
                    gradMu);
        }
    }

    @Override
    public void computeJointGradients(final CanonicalBranchTransitionProvider transitionProvider,
                                      final double[] gradA,
                                      final double[] gradQ,
                                      final double[] gradMu) {
        prepareBranchGradientInputs(transitionProvider, preparedBranchGradientInputs);
        computeJointGradients(transitionProvider, preparedBranchGradientInputs, gradA, gradQ, gradMu);
    }

    public BranchGradientInputs createBranchGradientInputs() {
        return new BranchGradientInputs(Math.max(0, nodeCount - 1), dim);
    }

    public BranchGradientInputs prepareBranchGradientInputs(
            final CanonicalBranchTransitionProvider transitionProvider) {
        final BranchGradientInputs out = createBranchGradientInputs();
        prepareBranchGradientInputs(transitionProvider, out);
        return out;
    }

    public void prepareBranchGradientInputs(final CanonicalBranchTransitionProvider transitionProvider,
                                            final BranchGradientInputs out) {
        ensureGradientState();
        out.checkCompatible(Math.max(0, nodeCount - 1), dim);
        out.clear();

        final HomogeneousCanonicalOUBranchTransitionProvider ouProvider = requireOUProvider(transitionProvider);
        final OUProcessModel processModel = ouProvider.getProcessModel();
        final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalSelection =
                processModel.getSelectionMatrixParameterization()
                        instanceof OrthogonalBlockDiagonalSelectionMatrixParameterization
                        ? (OrthogonalBlockDiagonalSelectionMatrixParameterization)
                        processModel.getSelectionMatrixParameterization()
                        : null;
        final double[] orthogonalStationaryMean =
                orthogonalSelection == null ? null : mainWorkspace.orthogonalStationaryMeanScratch;
        if (orthogonalSelection != null) {
            processModel.getInitialMean(orthogonalStationaryMean);
        }

        final int rootIndex = tree.getRoot().getNumber();
        if (branchGradientWorkspaces.length <= 1 || nodeCount <= 2) {
            for (int childIndex = 0; childIndex < nodeCount; ++childIndex) {
                if (childIndex == rootIndex) {
                    continue;
                }
                if (!fillLocalAdjointsForBranch(childIndex, transitionProvider, mainWorkspace)) {
                    continue;
                }
                final double branchLength = transitionProvider.getEffectiveBranchLength(childIndex);
                out.addBranch(
                        childIndex,
                        branchLength,
                        mainWorkspace.adjoints,
                        orthogonalSelection,
                        orthogonalStationaryMean);
            }
        } else {
            branchGradientTaskPool.forkDynamicBalanced(
                    nodeCount,
                    branchGradientPreparationChunkSize(nodeCount),
                    (childIndex, thread) -> {
                if (childIndex == rootIndex) {
                    out.clearStagedBranch(childIndex);
                    return;
                }

                final BranchGradientWorkspace workspace = branchGradientWorkspaces[thread];
                if (!fillLocalAdjointsForBranch(childIndex, transitionProvider, workspace)) {
                    out.clearStagedBranch(childIndex);
                    return;
                }

                out.stageBranch(
                        childIndex,
                        transitionProvider.getEffectiveBranchLength(childIndex),
                        workspace.adjoints,
                        null,
                        null);
                    });
            out.compactStagedBranches(rootIndex, false);
            out.prepareOrthogonalBasisForActiveBranches(orthogonalSelection, orthogonalStationaryMean);
        }

        out.rootDiffusionScale = lastRootDiffusionScale;
        copyState(preOrder[rootIndex], out.rootPreOrder);
        copyState(postOrder[rootIndex], out.rootPostOrder);
    }

    public void computeJointGradients(final CanonicalBranchTransitionProvider transitionProvider,
                                      final BranchGradientInputs inputs,
                                      final double[] gradA,
                                      final double[] gradQ,
                                      final double[] gradMu) {
        Arrays.fill(gradA, 0.0);
        Arrays.fill(gradQ, 0.0);
        Arrays.fill(gradMu, 0.0);
        inputs.checkCompatible(Math.max(0, nodeCount - 1), dim);

        final HomogeneousCanonicalOUBranchTransitionProvider ouProvider = requireOUProvider(transitionProvider);
        final OUProcessModel processModel = ouProvider.getProcessModel();
        final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalSelection =
                processModel.getSelectionMatrixParameterization()
                        instanceof OrthogonalBlockDiagonalSelectionMatrixParameterization
                        ? (OrthogonalBlockDiagonalSelectionMatrixParameterization)
                        processModel.getSelectionMatrixParameterization()
                        : null;
        if (branchGradientWorkspaces.length <= 1 || inputs.activeBranchCount <= 1) {
            computeJointGradientsSequential(inputs, processModel, orthogonalSelection, gradA, gradQ, gradMu);
            return;
        }

        computeJointGradientsParallel(inputs, processModel, orthogonalSelection, gradA, gradQ, gradMu);
    }

    private void computeJointGradientsSequential(final BranchGradientInputs inputs,
                                                 final OUProcessModel processModel,
                                                 final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalSelection,
                                                 final double[] gradA,
                                                 final double[] gradQ,
                                                 final double[] gradMu) {
        final BranchGradientWorkspace workspace = mainWorkspace;
        final OrthogonalGradientLayout orthogonalLayout =
                orthogonalSelection == null ? null : validateOrthogonalGradientLayout(orthogonalSelection, gradA, workspace);

        if (orthogonalLayout != null) {
            Arrays.fill(workspace.orthogonalCompressedGradientScratch, 0, orthogonalLayout.compressedBlockDim, 0.0);
            Arrays.fill(workspace.orthogonalNativeGradientScratch, 0, orthogonalLayout.nativeBlockDim, 0.0);
            clearSquare(workspace.orthogonalRotationGradientScratch);
            workspace.ensureOrthogonalBranchWorkspace(orthogonalSelection);
        }

        for (int activeIndex = 0; activeIndex < inputs.activeBranchCount; ++activeIndex) {
            copyAdjoints(inputs.activeAdjoints[activeIndex], workspace.adjoints);
            accumulateJointGradientForBranch(
                    inputs,
                    processModel,
                    orthogonalSelection,
                    activeIndex,
                    workspace,
                    gradA,
                    gradQ,
                    gradMu);
        }

        finalizeJointGradients(orthogonalLayout, inputs, workspace, gradA, gradQ);
    }

    private void computeJointGradientsParallel(final BranchGradientInputs inputs,
                                               final OUProcessModel processModel,
                                               final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalSelection,
                                               final double[] gradA,
                                               final double[] gradQ,
                                               final double[] gradMu) {
        final BranchGradientWorkspace reductionWorkspace = mainWorkspace;
        final OrthogonalGradientLayout orthogonalLayout =
                orthogonalSelection == null ? null : validateOrthogonalGradientLayout(orthogonalSelection, gradA, reductionWorkspace);

        for (int worker = 0; worker < branchGradientWorkspaces.length; ++worker) {
            final BranchGradientWorkspace workspace = branchGradientWorkspaces[worker];
            workspace.clearLocalGradientBuffers(
                    gradA.length,
                    gradMu.length,
                    dim,
                    orthogonalLayout != null,
                    orthogonalLayout == null ? 0 : orthogonalLayout.compressedBlockDim);
            if (orthogonalLayout != null) {
                workspace.ensureOrthogonalBranchWorkspace(orthogonalSelection);
            }
        }

        if (orthogonalLayout != null) {
            Arrays.fill(reductionWorkspace.orthogonalCompressedGradientScratch, 0, orthogonalLayout.compressedBlockDim, 0.0);
            Arrays.fill(reductionWorkspace.orthogonalNativeGradientScratch, 0, orthogonalLayout.nativeBlockDim, 0.0);
            clearSquare(reductionWorkspace.orthogonalRotationGradientScratch);
        }

        branchGradientTaskPool.forkDynamicBalanced(
                inputs.activeBranchCount,
                branchGradientJointChunkSize(inputs.activeBranchCount),
                (activeIndex, thread) -> {
            final BranchGradientWorkspace workspace = branchGradientWorkspaces[thread];
            copyAdjoints(inputs.activeAdjoints[activeIndex], workspace.adjoints);

            accumulateJointGradientForBranch(
                    inputs,
                    processModel,
                    orthogonalSelection,
                    activeIndex,
                    workspace,
                    workspace.localGradientA,
                    workspace.localGradientQ,
                    workspace.localGradientMu(gradMu.length, dim));
                });

        for (int worker = 0; worker < branchGradientWorkspaces.length; ++worker) {
            final BranchGradientWorkspace workspace = branchGradientWorkspaces[worker];
            accumulateVectorInPlace(gradQ, workspace.localGradientQ, gradQ.length);
            accumulateVectorInPlace(gradMu, workspace.localGradientMu(gradMu.length, dim), gradMu.length);
            if (orthogonalLayout != null) {
                accumulateVectorInPlace(
                        reductionWorkspace.orthogonalCompressedGradientScratch,
                        workspace.orthogonalCompressedGradientScratch,
                        orthogonalLayout.compressedBlockDim);
                accumulateSquareInPlace(
                        reductionWorkspace.orthogonalRotationGradientScratch,
                        workspace.orthogonalRotationGradientScratch);
            } else {
                accumulateVectorInPlace(gradA, workspace.localGradientA, gradA.length);
            }
        }

        finalizeJointGradients(orthogonalLayout, inputs, reductionWorkspace, gradA, gradQ);
    }

    private void accumulateJointGradientForBranch(final BranchGradientInputs inputs,
                                                  final OUProcessModel processModel,
                                                  final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalSelection,
                                                  final int activeIndex,
                                                  final BranchGradientWorkspace workspace,
                                                  final double[] gradA,
                                                  final double[] gradQ,
                                                  final double[] gradMu) {
        if (orthogonalSelection != null) {
            accumulateOrthogonalJointGradientForBranch(
                    inputs.getOrthogonalPreparedBasis(activeIndex),
                    processModel,
                    orthogonalSelection,
                    workspace,
                    gradQ,
                    gradMu);
            return;
        }

        final double branchLength = inputs.activeBranchLengths[activeIndex];
        accumulateDenseJointGradientForBranch(processModel, workspace, branchLength, gradA, gradQ, gradMu);
    }

    private void accumulateOrthogonalJointGradientForBranch(
                                                            final OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis preparedBasis,
                                                            final OUProcessModel processModel,
                                                            final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalSelection,
                                                            final BranchGradientWorkspace workspace,
                                                            final double[] gradQ,
                                                            final double[] gradMu) {
        final OrthogonalBlockDiagonalSelectionMatrixParameterization.BranchGradientWorkspace orthogonalWorkspace =
                workspace.ensureOrthogonalBranchWorkspace(orthogonalSelection);
        orthogonalSelection.accumulateNativeGradientFromAdjointsPrepared(
                preparedBasis,
                processModel.getDiffusionMatrix(),
                workspace.adjoints,
                orthogonalWorkspace,
                workspace.orthogonalCompressedGradientScratch,
                workspace.orthogonalRotationGradientScratch);

        transposeFromFlatInto(workspace.adjoints.dLogL_dOmega, workspace.covarianceAdjoint, dim);
        orthogonalSelection.accumulateDiffusionGradientPrepared(
                preparedBasis,
                workspace.covarianceAdjoint,
                gradQ,
                orthogonalWorkspace);

        orthogonalSelection.accumulateMeanGradientPrepared(
                preparedBasis,
                workspace.adjoints.dLogL_df,
                gradMu,
                orthogonalWorkspace);
    }

    private void accumulateDenseJointGradientForBranch(final OUProcessModel processModel,
                                                       final BranchGradientWorkspace workspace,
                                                       final double branchLength,
                                                       final double[] gradA,
                                                       final double[] gradQ,
                                                       final double[] gradMu) {
        copyFlatToSquare(workspace.adjoints.dLogL_dF, workspace.transitionMatrix, dim);
        processModel.accumulateSelectionGradient(
                branchLength,
                workspace.transitionMatrix,
                workspace.adjoints.dLogL_df,
                gradA);

        transposeFromFlatInto(workspace.adjoints.dLogL_dOmega, workspace.covarianceAdjoint, dim);
        processModel.accumulateSelectionGradientFromCovariance(
                branchLength,
                workspace.covarianceAdjoint,
                gradA);

        copyFlatToSquare(workspace.adjoints.dLogL_dOmega, workspace.covarianceAdjoint, dim);
        processModel.accumulateDiffusionGradient(
                branchLength,
                workspace.covarianceAdjoint,
                gradQ);

        processModel.fillTransitionMatrix(branchLength, workspace.transitionMatrix);
        accumulateStationaryMeanGradient(
                workspace.transitionMatrix,
                workspace.adjoints.dLogL_df,
                gradMu);
    }

    private void finalizeJointGradients(final OrthogonalGradientLayout orthogonalLayout,
                                        final BranchGradientInputs inputs,
                                        final BranchGradientWorkspace workspace,
                                        final double[] gradA,
                                        final double[] gradQ) {
        if (orthogonalLayout != null) {
            orthogonalLayout.blockParameter.chainGradient(
                    workspace.orthogonalCompressedGradientScratch,
                    workspace.orthogonalNativeGradientScratch);
            final double[] angleGradient =
                    orthogonalLayout.orthogonalRotation.pullBackGradient(workspace.orthogonalRotationGradientScratch);
            System.arraycopy(workspace.orthogonalNativeGradientScratch, 0, gradA, 0, orthogonalLayout.nativeBlockDim);
            System.arraycopy(angleGradient, 0, gradA, orthogonalLayout.nativeBlockDim, angleGradient.length);
        } else {
            transposeFlatSquareInPlace(gradA, dim);
        }

        if (inputs.rootDiffusionScale > 0.0) {
            accumulateRootDiffusionGradient(
                    inputs.rootPreOrder,
                    inputs.rootPostOrder,
                    inputs.rootDiffusionScale,
                    gradQ,
                    workspace);
        }
    }

    private OrthogonalGradientLayout validateOrthogonalGradientLayout(
            final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalSelection,
            final double[] gradA,
            final BranchGradientWorkspace workspace) {
        final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter =
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) orthogonalSelection.getMatrixParameter();
        if (!(blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider)) {
            throw new IllegalStateException(
                    "Orthogonal block native gradient requires an OrthogonalMatrixProvider rotation parameter.");
        }

        final OrthogonalMatrixProvider orthogonalRotation =
                (OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter();
        final int nativeBlockDim = blockParameter.getBlockDiagonalNParameters();
        final int angleDim = orthogonalRotation.getOrthogonalParameter().getDimension();
        final int nativeDim = nativeBlockDim + angleDim;
        if (gradA.length != nativeDim) {
            throw new IllegalArgumentException(
                    "Orthogonal block selection gradient expects native parameter length "
                            + nativeDim + ", found " + gradA.length);
        }

        final int compressedBlockDim = blockParameter.getCompressedDDimension();
        if (compressedBlockDim > workspace.orthogonalCompressedGradientScratch.length
                || nativeBlockDim > workspace.orthogonalNativeGradientScratch.length) {
            throw new IllegalStateException(
                    "Orthogonal block scratch is too small for native gradient dimensions "
                            + compressedBlockDim + " and " + nativeBlockDim + ".");
        }
        return new OrthogonalGradientLayout(
                blockParameter,
                orthogonalRotation,
                compressedBlockDim,
                nativeBlockDim);
    }

    private static final class OrthogonalGradientLayout {
        final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
        final OrthogonalMatrixProvider orthogonalRotation;
        final int compressedBlockDim;
        final int nativeBlockDim;

        private OrthogonalGradientLayout(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                         final OrthogonalMatrixProvider orthogonalRotation,
                                         final int compressedBlockDim,
                                         final int nativeBlockDim) {
            this.blockParameter = blockParameter;
            this.orthogonalRotation = orthogonalRotation;
            this.compressedBlockDim = compressedBlockDim;
            this.nativeBlockDim = nativeBlockDim;
        }
    }

    @Override
    public void storeState() {
        storedHasPostOrderState = hasPostOrderState;
        storedHasPreOrderState = hasPreOrderState;
        storedLastRootDiffusionScale = lastRootDiffusionScale;
        storedHasFixedRootValue = hasFixedRootValue;
        System.arraycopy(fixedRootValue, 0, storedFixedRootValue, 0, dim);
        for (int i = 0; i < nodeCount; i++) {
            storedTipObservations[i].copyFrom(tipObservations[i]);
            copyState(postOrder[i], storedPostOrder[i]);
            copyState(preOrder[i], storedPreOrder[i]);
            copyState(branchAboveParent[i], storedBranchAboveParent[i]);
        }
    }

    @Override
    public void restoreState() {
        final CanonicalGaussianState[] tmpPost = postOrder;
        postOrder = storedPostOrder;
        storedPostOrder = tmpPost;

        final CanonicalGaussianState[] tmpPre = preOrder;
        preOrder = storedPreOrder;
        storedPreOrder = tmpPre;

        final CanonicalGaussianState[] tmpAbove = branchAboveParent;
        branchAboveParent = storedBranchAboveParent;
        storedBranchAboveParent = tmpAbove;

        final boolean tmpHasPost = hasPostOrderState;
        hasPostOrderState = storedHasPostOrderState;
        storedHasPostOrderState = tmpHasPost;

        final boolean tmpHasPre = hasPreOrderState;
        hasPreOrderState = storedHasPreOrderState;
        storedHasPreOrderState = tmpHasPre;

        final double tmpRootScale = lastRootDiffusionScale;
        lastRootDiffusionScale = storedLastRootDiffusionScale;
        storedLastRootDiffusionScale = tmpRootScale;

        final boolean tmpHasFixedRoot = hasFixedRootValue;
        hasFixedRootValue = storedHasFixedRootValue;
        storedHasFixedRootValue = tmpHasFixedRoot;

        final double[] tmpRootValue = fixedRootValue;
        fixedRootValue = storedFixedRootValue;
        storedFixedRootValue = tmpRootValue;

        final CanonicalTipObservation[] tmpTips = tipObservations;
        tipObservations = storedTipObservations;
        storedTipObservations = tmpTips;
    }

    @Override
    public void acceptState() { }

    private static CanonicalGaussianState[] allocateStates(final int count, final int dim) {
        final CanonicalGaussianState[] out = new CanonicalGaussianState[count];
        for (int i = 0; i < count; i++) {
            out[i] = new CanonicalGaussianState(dim);
        }
        return out;
    }

    private static CanonicalTipObservation[] allocateTipObservations(final int count, final int dim) {
        final CanonicalTipObservation[] out = new CanonicalTipObservation[count];
        for (int i = 0; i < count; i++) {
            out[i] = new CanonicalTipObservation(dim);
        }
        return out;
    }

    private static void clearState(final CanonicalGaussianState state) {
        CanonicalGaussianMessageOps.clearState(state);
    }

    private static void copyState(final CanonicalGaussianState source, final CanonicalGaussianState target) {
        CanonicalGaussianMessageOps.copyState(source, target);
    }

    private void computePostOrderAtInternalNode(final int nodeIndex,
                                                final CanonicalBranchTransitionProvider transitionProvider,
                                                final BranchGradientWorkspace workspace) {
        final CanonicalGaussianState dest = postOrder[nodeIndex];
        clearState(dest);

        final int childCount = tree.getChildCount(tree.getNode(nodeIndex));
        boolean first = true;
        for (int c = 0; c < childCount; c++) {
            final int childIndex = tree.getChild(tree.getNode(nodeIndex), c).getNumber();
            buildUpwardParentMessage(childIndex, transitionProvider, workspace.state, workspace);
            if (first) {
                copyState(workspace.state, dest);
                first = false;
            } else {
                CanonicalGaussianMessageOps.combineStateInPlace(dest, workspace.state);
            }
        }
    }

    private void computePreOrderRecursive(final int parentIndex,
                                          final CanonicalBranchTransitionProvider transitionProvider,
                                          final BranchGradientWorkspace workspace) {
        if (tree.isExternal(tree.getNode(parentIndex))) {
            return;
        }

        final int rootIndex = tree.getRoot().getNumber();
        final int childCount = tree.getChildCount(tree.getNode(parentIndex));
        for (int c = 0; c < childCount; c++) {
            final int childIndex = tree.getChild(tree.getNode(parentIndex), c).getNumber();
            if (hasFixedRootValue && parentIndex == rootIndex) {
                clearState(branchAboveParent[childIndex]);
                transitionProvider.fillCanonicalTransition(childIndex, workspace.transition);
                CanonicalGaussianMessageOps.conditionOnObservedFirstBlock(
                        workspace.transition,
                        fixedRootValue,
                        preOrder[childIndex]);
                computePreOrderRecursive(childIndex, transitionProvider, workspace);
                continue;
            }

            final boolean hasSiblings = buildSiblingProduct(
                    parentIndex, childIndex, transitionProvider, workspace.siblingProduct, workspace);

            if (hasSiblings) {
                CanonicalGaussianMessageOps.combineStates(
                        preOrder[parentIndex],
                        workspace.siblingProduct,
                        workspace.downwardParentState);
            } else {
                copyState(preOrder[parentIndex], workspace.downwardParentState);
            }

            copyState(workspace.downwardParentState, branchAboveParent[childIndex]);
            transitionProvider.fillCanonicalTransition(childIndex, workspace.transition);
            CanonicalGaussianMessageOps.pushForward(
                    workspace.downwardParentState,
                    workspace.transition,
                    workspace.gaussianWorkspace,
                    preOrder[childIndex]);
            computePreOrderRecursive(childIndex, transitionProvider, workspace);
        }
    }

    private boolean buildSiblingProduct(final int parentIndex,
                                        final int excludedChildIndex,
                                        final CanonicalBranchTransitionProvider transitionProvider,
                                        final CanonicalGaussianState out,
                                        final BranchGradientWorkspace workspace) {
        clearState(out);
        final int childCount = tree.getChildCount(tree.getNode(parentIndex));
        boolean found = false;
        for (int c = 0; c < childCount; c++) {
            final int siblingIndex = tree.getChild(tree.getNode(parentIndex), c).getNumber();
            if (siblingIndex == excludedChildIndex) {
                continue;
            }
            buildUpwardParentMessage(siblingIndex, transitionProvider, workspace.state, workspace);
            if (!found) {
                copyState(workspace.state, out);
                found = true;
            } else {
                CanonicalGaussianMessageOps.combineStateInPlace(out, workspace.state);
            }
        }
        return found;
    }

    private void buildUpwardParentMessage(final int childIndex,
                                          final CanonicalBranchTransitionProvider transitionProvider,
                                          final CanonicalGaussianState out,
                                          final BranchGradientWorkspace workspace) {
        if (tree.isExternal(tree.getNode(childIndex))) {
            final CanonicalTipObservation tipObservation = tipObservations[childIndex];
            if (tipObservation.observedCount < dim) {
                buildTipParentMessage(childIndex, tipObservation, transitionProvider, out, workspace);
                return;
            }
        }
        transitionProvider.fillCanonicalTransition(childIndex, workspace.transition);
        buildUpwardParentMessageForTransition(childIndex, workspace.transition, out, workspace);
    }

    private void buildUpwardParentMessageForTransition(final int childIndex,
                                                       final CanonicalGaussianTransition transition,
                                                       final CanonicalGaussianState out,
                                                       final BranchGradientWorkspace workspace) {
        if (tree.isExternal(tree.getNode(childIndex))) {
            final CanonicalTipObservation tipObservation = tipObservations[childIndex];
            if (tipObservation.observedCount == 0) {
                clearState(out);
            } else if (tipObservation.observedCount == dim) {
                CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transition, tipObservation.values, out);
            } else {
                throw new IllegalStateException(
                        "Partially observed canonical tips must use the missing-mask branch update.");
            }
        } else {
            CanonicalGaussianMessageOps.pushBackward(postOrder[childIndex], transition, workspace.gaussianWorkspace, out);
        }
    }

    private void buildTipParentMessage(final int childIndex,
                                       final CanonicalTipObservation tipObservation,
                                       final CanonicalBranchTransitionProvider transitionProvider,
                                       final CanonicalGaussianState out,
                                       final BranchGradientWorkspace workspace) {
        if (tipObservation.observedCount == 0) {
            clearState(out);
            return;
        }

        if (!(transitionProvider instanceof HomogeneousCanonicalOUBranchTransitionProvider)) {
            throw new UnsupportedOperationException(
                    "Not yet implemented: canonical OU missing-tip propagation currently supports only "
                            + "HomogeneousCanonicalOUBranchTransitionProvider.");
        }

        final HomogeneousCanonicalOUBranchTransitionProvider ouProvider =
                (HomogeneousCanonicalOUBranchTransitionProvider) transitionProvider;
        buildTipParentMessageWithMissingMask(
                tipObservation,
                ouProvider.getProcessModel(),
                transitionProvider.getEffectiveBranchLength(childIndex),
                out,
                workspace);
    }

    private void buildTipParentMessageWithMissingMask(final CanonicalTipObservation tipObservation,
                                                      final OUProcessModel processModel,
                                                      final double branchLength,
                                                      final CanonicalGaussianState out,
                                                      final BranchGradientWorkspace workspace) {
        final int observedCount = collectObservationPartition(tipObservation, workspace);
        if (observedCount == 0) {
            clearState(out);
            return;
        }

        processModel.fillTransitionMatrix(branchLength, workspace.transitionMatrix);
        processModel.fillTransitionOffset(branchLength, workspace.mean);
        processModel.fillTransitionCovariance(branchLength, workspace.covariance);

        for (int observed = 0; observed < observedCount; ++observed) {
            final int observedTrait = workspace.observedIndexScratch[observed];
            workspace.mean2[observed] = tipObservation.values[observedTrait] - workspace.mean[observedTrait];
            final int rowOffset = observed * observedCount;
            for (int otherObserved = 0; otherObserved < observedCount; ++otherObserved) {
                workspace.varianceFlat[rowOffset + otherObserved] =
                        workspace.covariance[observedTrait][workspace.observedIndexScratch[otherObserved]];
            }
        }

        final double logDetVariance = MatrixUtils.invertSymmetricPositiveDefiniteCompact(
                workspace.varianceFlat,
                workspace.precisionFlat,
                observedCount,
                workspace.mean,
                workspace.reducedMeanScratch);

        for (int row = 0; row < observedCount; ++row) {
            double sum = 0.0;
            final int rowOffset = row * observedCount;
            for (int k = 0; k < observedCount; ++k) {
                sum += workspace.precisionFlat[rowOffset + k] * workspace.mean2[k];
            }
            workspace.mean[row] = sum;
        }

        for (int row = 0; row < observedCount; ++row) {
            final int observedTrait = workspace.observedIndexScratch[row];
            final int precisionRowOffset = row * observedCount;
            final int projectedRowOffset = row * dim;
            for (int col = 0; col < dim; ++col) {
                double sum = 0.0;
                for (int k = 0; k < observedCount; ++k) {
                    sum += workspace.precisionFlat[precisionRowOffset + k]
                            * workspace.transitionMatrix[workspace.observedIndexScratch[k]][col];
                }
                workspace.reducedPrecisionFlatScratch[projectedRowOffset + col] = sum;
            }
        }

        for (int i = 0; i < dim; ++i) {
            double information = 0.0;
            for (int observed = 0; observed < observedCount; ++observed) {
                information += workspace.transitionMatrix[workspace.observedIndexScratch[observed]][i]
                        * workspace.mean[observed];
            }
            out.information[i] = information;

            for (int j = 0; j < dim; ++j) {
                double precision = 0.0;
                for (int observed = 0; observed < observedCount; ++observed) {
                    precision += workspace.transitionMatrix[workspace.observedIndexScratch[observed]][i]
                            * workspace.reducedPrecisionFlatScratch[observed * dim + j];
                }
                out.precision[i * dim + j] = precision;
            }
        }
        symmetrizeFlatSquare(out.precision, dim);

        double quadratic = 0.0;
        for (int observed = 0; observed < observedCount; ++observed) {
            quadratic += workspace.mean2[observed] * workspace.mean[observed];
        }
        out.logNormalizer = 0.5 * (
                observedCount * LOG_TWO_PI
                        + logDetVariance
                        + quadratic);
    }

    private int collectObservationPartition(final CanonicalTipObservation tipObservation,
                                            final BranchGradientWorkspace workspace) {
        int observedCount = 0;
        int missingCount = 0;

        for (int i = 0; i < dim; i++) {
            if (tipObservation.observed[i]) {
                workspace.observedIndexScratch[observedCount++] = i;
                workspace.reducedIndexByTraitScratch[i] = -1;
            } else {
                workspace.missingIndexScratch[missingCount++] = i;
                workspace.reducedIndexByTraitScratch[i] = dim + missingCount - 1;
            }
        }

        if (observedCount != tipObservation.observedCount || observedCount + missingCount != dim) {
            throw new UnsupportedOperationException(
                    "Canonical tip observation partition is inconsistent with observedCount.");
        }
        return observedCount;
    }

    private void ensureGradientState() {
        if (!hasPostOrderState || !hasPreOrderState) {
            throw new IllegalStateException(
                    "Canonical gradients require both computePostOrderLogLikelihood and computePreOrder to have been called.");
        }
    }

    private HomogeneousCanonicalOUBranchTransitionProvider requireOUProvider(
            final CanonicalBranchTransitionProvider transitionProvider) {
        if (!(transitionProvider instanceof HomogeneousCanonicalOUBranchTransitionProvider)) {
            throw new UnsupportedOperationException(
                    "Not yet implemented: canonical OU gradients currently support only "
                            + "HomogeneousCanonicalOUBranchTransitionProvider; heterogeneous "
                            + "CanonicalBranchTransitionProvider implementations are not yet supported.");
        }
        return (HomogeneousCanonicalOUBranchTransitionProvider) transitionProvider;
    }

    private boolean fillLocalAdjointsForBranch(final int childIndex,
                                               final CanonicalBranchTransitionProvider transitionProvider,
                                               final BranchGradientWorkspace workspace) {
        transitionProvider.fillCanonicalTransition(childIndex, workspace.transition);
        if (hasFixedRootValue && tree.isRoot(tree.getParent(tree.getNode(childIndex)))) {
            if (tree.isExternal(tree.getNode(childIndex))) {
                final CanonicalTipObservation tipObservation = tipObservations[childIndex];
                if (tipObservation.observedCount == 0) {
                    return false;
                }
                final int observedCount = collectObservationPartition(tipObservation, workspace);
                if (observedCount == dim) {
                    fillContributionForFixedParentObservedTip(workspace.transition, tipObservation, workspace);
                } else {
                    fillContributionForFixedParentPartiallyObservedTip(
                            workspace.transition, tipObservation, observedCount, workspace);
                }
            } else {
                fillContributionForFixedParentInternalNode(workspace.transition, postOrder[childIndex], workspace);
            }

            CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                    workspace.transition,
                    workspace.contribution,
                    workspace.transitionAdjointWorkspace,
                    workspace.adjoints);
            return true;
        }

        if (tree.isExternal(tree.getNode(childIndex))) {
            final CanonicalTipObservation tipObservation = tipObservations[childIndex];
            if (tipObservation.observedCount == 0) {
                return false;
            }
            final int observedCount = collectObservationPartition(tipObservation, workspace);
            if (observedCount == dim) {
                fillContributionForObservedTip(branchAboveParent[childIndex], workspace.transition, tipObservation, workspace);
            } else {
                fillContributionForPartiallyObservedTip(
                        branchAboveParent[childIndex], workspace.transition, tipObservation, observedCount, workspace);
            }
        } else {
            CanonicalGaussianMessageOps.buildPairPosterior(
                    branchAboveParent[childIndex],
                    workspace.transition,
                    postOrder[childIndex],
                    workspace.pairState);
            CanonicalBranchMessageContributionUtils.fillFromPairState(
                    workspace.pairState,
                    workspace.contributionWorkspace,
                    workspace.contribution);
        }

        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                workspace.transition,
                workspace.contribution,
                workspace.transitionAdjointWorkspace,
                workspace.adjoints);
        return true;
    }

    private void fillContributionForObservedTip(final CanonicalGaussianState aboveState,
                                                final CanonicalGaussianTransition transition,
                                                final CanonicalTipObservation tipObservation,
                                                final BranchGradientWorkspace workspace) {
        for (int i = 0; i < dim; ++i) {
            double info = transition.informationX[i] + aboveState.information[i];
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                info -= transition.precisionXY[iOff + j] * tipObservation.values[j];
                workspace.parentPosterior.precision[iOff + j] =
                        transition.precisionXX[iOff + j] + aboveState.precision[iOff + j];
            }
            workspace.parentPosterior.information[i] = info;
        }
        workspace.parentPosterior.logNormalizer = 0.0;

        CanonicalGaussianUtils.fillMomentsFromCanonical(
                workspace.parentPosterior, workspace.mean, workspace.covariance);

        for (int i = 0; i < dim; ++i) {
            final double xi = workspace.mean[i];
            final double yi = tipObservation.values[i];
            workspace.contribution.dLogL_dInformationX[i] = xi;
            workspace.contribution.dLogL_dInformationY[i] = yi;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = workspace.mean[j];
                final double yj = tipObservation.values[j];
                final double exx = workspace.covariance[i][j] + xi * xj;
                workspace.contribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * exx;
                workspace.contribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                workspace.contribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                workspace.contribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * (yi * yj);
            }
        }
        workspace.contribution.dLogL_dLogNormalizer = -1.0;
    }

    private void fillContributionForFixedParentObservedTip(final CanonicalGaussianTransition transition,
                                                           final CanonicalTipObservation tipObservation,
                                                           final BranchGradientWorkspace workspace) {
        clearContribution(workspace);
        for (int i = 0; i < dim; ++i) {
            final double xi = fixedRootValue[i];
            final double yi = tipObservation.values[i];
            workspace.contribution.dLogL_dInformationX[i] = xi;
            workspace.contribution.dLogL_dInformationY[i] = yi;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = fixedRootValue[j];
                final double yj = tipObservation.values[j];
                workspace.contribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * (xi * xj);
                workspace.contribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                workspace.contribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                workspace.contribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * (yi * yj);
            }
        }
        workspace.contribution.dLogL_dLogNormalizer = -1.0;
    }

    private void fillContributionForFixedParentInternalNode(final CanonicalGaussianTransition transition,
                                                            final CanonicalGaussianState childMessage,
                                                            final BranchGradientWorkspace workspace) {
        CanonicalGaussianMessageOps.conditionOnObservedFirstBlock(transition, fixedRootValue, workspace.state);
        CanonicalGaussianMessageOps.combineStates(workspace.state, childMessage, workspace.parentPosterior);
        CanonicalGaussianUtils.fillMomentsFromCanonical(
                workspace.parentPosterior, workspace.mean, workspace.covariance);
        fillContributionFromFixedParentChildMoments(workspace.mean, workspace.covariance, workspace);
    }

    private void fillContributionForFixedParentPartiallyObservedTip(final CanonicalGaussianTransition transition,
                                                                    final CanonicalTipObservation tipObservation,
                                                                    final int observedCount,
                                                                    final BranchGradientWorkspace workspace) {
        final int missingCount = dim - observedCount;
        CanonicalGaussianMessageOps.conditionOnObservedFirstBlock(transition, fixedRootValue, workspace.state);

        for (int missing = 0; missing < missingCount; ++missing) {
            final int missingTrait = workspace.missingIndexScratch[missing];
            double info = workspace.state.information[missingTrait];
            final int missingTraitOff = missingTrait * dim;
            for (int observed = 0; observed < observedCount; ++observed) {
                final int observedTrait = workspace.observedIndexScratch[observed];
                info -= workspace.state.precision[missingTraitOff + observedTrait] * tipObservation.values[observedTrait];
            }
            workspace.reducedInformationScratch[missing] = info;
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                workspace.reducedPrecisionFlatScratch[missing * missingCount + otherMissing] =
                        workspace.state.precision[missingTraitOff + workspace.missingIndexScratch[otherMissing]];
            }
        }

        MatrixUtils.safeInvertPrecision(
                workspace.reducedPrecisionFlatScratch,
                workspace.reducedCovarianceFlatScratch,
                missingCount);

        for (int missing = 0; missing < missingCount; ++missing) {
            double sum = 0.0;
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                sum += workspace.reducedCovarianceFlatScratch[missing * missingCount + otherMissing]
                        * workspace.reducedInformationScratch[otherMissing];
            }
            workspace.mean2[missing] = sum;
        }

        clearContribution(workspace);
        for (int i = 0; i < dim; ++i) {
            workspace.contribution.dLogL_dInformationX[i] = fixedRootValue[i];
            workspace.contribution.dLogL_dInformationY[i] = tipObservation.observed[i]
                    ? tipObservation.values[i]
                    : workspace.mean2[workspace.reducedIndexByTraitScratch[i] - dim];
        }

        for (int i = 0; i < dim; ++i) {
            final double xi = fixedRootValue[i];
            final double yi = workspace.contribution.dLogL_dInformationY[i];
            final int missingI = tipObservation.observed[i] ? -1 : workspace.reducedIndexByTraitScratch[i] - dim;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = fixedRootValue[j];
                final double yj = workspace.contribution.dLogL_dInformationY[j];
                final int missingJ = tipObservation.observed[j] ? -1 : workspace.reducedIndexByTraitScratch[j] - dim;
                final double eyy = (missingI < 0 || missingJ < 0)
                        ? yi * yj
                        : workspace.reducedCovarianceFlatScratch[missingI * missingCount + missingJ] + yi * yj;

                workspace.contribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * (xi * xj);
                workspace.contribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                workspace.contribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                workspace.contribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * eyy;
            }
        }
        workspace.contribution.dLogL_dLogNormalizer = -1.0;
    }

    private void fillContributionFromFixedParentChildMoments(final double[] childMean,
                                                             final double[][] childCovariance,
                                                             final BranchGradientWorkspace workspace) {
        clearContribution(workspace);
        for (int i = 0; i < dim; ++i) {
            final double xi = fixedRootValue[i];
            final double yi = childMean[i];
            workspace.contribution.dLogL_dInformationX[i] = xi;
            workspace.contribution.dLogL_dInformationY[i] = yi;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = fixedRootValue[j];
                final double yj = childMean[j];
                final double eyy = childCovariance[i][j] + yi * yj;
                workspace.contribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * (xi * xj);
                workspace.contribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                workspace.contribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                workspace.contribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * eyy;
            }
        }
        workspace.contribution.dLogL_dLogNormalizer = -1.0;
    }

    private void fillContributionForPartiallyObservedTip(final CanonicalGaussianState aboveState,
                                                         final CanonicalGaussianTransition transition,
                                                         final CanonicalTipObservation tipObservation,
                                                         final int observedCount,
                                                         final BranchGradientWorkspace workspace) {
        final int missingCount = dim - observedCount;
        final int reducedDimension = dim + missingCount;

        for (int i = 0; i < dim; ++i) {
            double info = transition.informationX[i] + aboveState.information[i];
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                workspace.reducedPrecisionFlatScratch[i * reducedDimension + j] =
                        transition.precisionXX[iOff + j] + aboveState.precision[iOff + j];
            }
            for (int observed = 0; observed < observedCount; ++observed) {
                final int observedIndex = workspace.observedIndexScratch[observed];
                info -= transition.precisionXY[iOff + observedIndex] * tipObservation.values[observedIndex];
            }
            workspace.reducedInformationScratch[i] = info;
            for (int missing = 0; missing < missingCount; ++missing) {
                workspace.reducedPrecisionFlatScratch[i * reducedDimension + dim + missing] =
                        transition.precisionXY[iOff + workspace.missingIndexScratch[missing]];
            }
        }

        for (int missing = 0; missing < missingCount; ++missing) {
            final int childIndex = workspace.missingIndexScratch[missing];
            final int row = dim + missing;
            final int childOff = childIndex * dim;
            double info = transition.informationY[childIndex];
            for (int observed = 0; observed < observedCount; ++observed) {
                final int observedIndex = workspace.observedIndexScratch[observed];
                info -= transition.precisionYY[childOff + observedIndex] * tipObservation.values[observedIndex];
            }
            workspace.reducedInformationScratch[row] = info;
            for (int j = 0; j < dim; ++j) {
                workspace.reducedPrecisionFlatScratch[row * reducedDimension + j] =
                        transition.precisionYX[childOff + j];
            }
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                workspace.reducedPrecisionFlatScratch[row * reducedDimension + dim + otherMissing] =
                        transition.precisionYY[childOff + workspace.missingIndexScratch[otherMissing]];
            }
        }

        MatrixUtils.safeInvertPrecision(
                workspace.reducedPrecisionFlatScratch,
                workspace.reducedCovarianceFlatScratch,
                reducedDimension);

        for (int i = 0; i < reducedDimension; ++i) {
            double sum = 0.0;
            for (int j = 0; j < reducedDimension; ++j) {
                sum += workspace.reducedCovarianceFlatScratch[i * reducedDimension + j]
                        * workspace.reducedInformationScratch[j];
            }
            workspace.reducedMeanScratch[i] = Double.isNaN(sum) ? 0.0 : sum;
        }

        clearContribution(workspace);
        for (int i = 0; i < dim; ++i) {
            workspace.contribution.dLogL_dInformationX[i] = workspace.reducedMeanScratch[i];
            workspace.contribution.dLogL_dInformationY[i] = tipObservation.observed[i]
                    ? tipObservation.values[i]
                    : workspace.reducedMeanScratch[workspace.reducedIndexByTraitScratch[i]];
        }

        for (int i = 0; i < dim; ++i) {
            final double xi = workspace.reducedMeanScratch[i];
            final int reducedI = workspace.reducedIndexByTraitScratch[i];
            for (int j = 0; j < dim; ++j) {
                final double xj = workspace.reducedMeanScratch[j];
                final int reducedJ = workspace.reducedIndexByTraitScratch[j];
                final double yi = workspace.contribution.dLogL_dInformationY[i];
                final double yj = workspace.contribution.dLogL_dInformationY[j];

                final double exx = workspace.reducedCovarianceFlatScratch[i * reducedDimension + j] + xi * xj;
                final double exy = tipObservation.observed[j]
                        ? xi * yj
                        : workspace.reducedCovarianceFlatScratch[i * reducedDimension + reducedJ] + xi * yj;
                final double eyx = tipObservation.observed[i]
                        ? yi * xj
                        : workspace.reducedCovarianceFlatScratch[reducedI * reducedDimension + j] + yi * xj;
                final double eyy = (tipObservation.observed[i] || tipObservation.observed[j])
                        ? yi * yj
                        : workspace.reducedCovarianceFlatScratch[reducedI * reducedDimension + reducedJ] + yi * yj;

                final int ij = i * dim + j;
                workspace.contribution.dLogL_dPrecisionXX[ij] = -0.5 * exx;
                workspace.contribution.dLogL_dPrecisionXY[ij] = -0.5 * exy;
                workspace.contribution.dLogL_dPrecisionYX[ij] = -0.5 * eyx;
                workspace.contribution.dLogL_dPrecisionYY[ij] = -0.5 * eyy;
            }
        }
        workspace.contribution.dLogL_dLogNormalizer = -1.0;
    }

    private void clearContribution(final BranchGradientWorkspace workspace) {
        Arrays.fill(workspace.contribution.dLogL_dInformationX, 0.0);
        Arrays.fill(workspace.contribution.dLogL_dInformationY, 0.0);
        Arrays.fill(workspace.contribution.dLogL_dPrecisionXX, 0.0);
        Arrays.fill(workspace.contribution.dLogL_dPrecisionXY, 0.0);
        Arrays.fill(workspace.contribution.dLogL_dPrecisionYX, 0.0);
        Arrays.fill(workspace.contribution.dLogL_dPrecisionYY, 0.0);
        workspace.contribution.dLogL_dLogNormalizer = -1.0;
    }

    private double finiteDifferenceBranchLengthGradient(
            final int childIndex,
            final HomogeneousCanonicalOUBranchTransitionProvider provider,
            final double branchLength) {
        final double step = Math.max(BRANCH_LENGTH_FD_ABSOLUTE_STEP,
                BRANCH_LENGTH_FD_RELATIVE_STEP * Math.max(1.0, Math.abs(branchLength)));

        if (branchLength > step) {
            final double plus = evaluateFrozenLocalLogFactor(childIndex, provider, branchLength + step);
            final double minus = evaluateFrozenLocalLogFactor(childIndex, provider, branchLength - step);
            return (plus - minus) / (2.0 * step);
        }

        final double plus = evaluateFrozenLocalLogFactor(childIndex, provider, branchLength + step);
        final double base = evaluateFrozenLocalLogFactor(childIndex, provider, branchLength);
        return (plus - base) / step;
    }

    private double evaluateFrozenLocalLogFactor(final int childIndex,
                                                final HomogeneousCanonicalOUBranchTransitionProvider provider,
                                                final double branchLength) {
        final BranchGradientWorkspace workspace = mainWorkspace;
        final CanonicalTipObservation tipObservation =
                tree.isExternal(tree.getNode(childIndex)) ? tipObservations[childIndex] : null;
        if (tipObservation != null && tipObservation.observedCount < dim) {
            buildTipParentMessageWithMissingMask(
                    tipObservation,
                    provider.getProcessModel(),
                    branchLength,
                    workspace.state,
                    workspace);
        } else {
            provider.getProcessModel().fillCanonicalTransition(branchLength, workspace.transition);
            buildUpwardParentMessageForTransition(childIndex, workspace.transition, workspace.state, workspace);
        }
        CanonicalGaussianMessageOps.combineStates(
                branchAboveParent[childIndex], workspace.state, workspace.combinedState);
        return CanonicalGaussianMessageOps.normalizationShift(workspace.combinedState, workspace.gaussianWorkspace);
    }

    private void accumulateStationaryMeanGradient(final double[][] transitionMatrix,
                                                  final double[] adjointB,
                                                  final double[] gradient) {
        if (gradient.length == 1) {
            double sum = 0.0;
            for (int i = 0; i < dim; ++i) {
                double ftAdjoint = 0.0;
                for (int j = 0; j < dim; ++j) {
                    ftAdjoint += transitionMatrix[j][i] * adjointB[j];
                }
                sum += adjointB[i] - ftAdjoint;
            }
            gradient[0] += sum;
            return;
        }
        if (gradient.length != dim) {
            throw new IllegalArgumentException(
                    "Stationary-mean gradient length must be 1 or " + dim + ", found " + gradient.length);
        }

        for (int i = 0; i < dim; ++i) {
            double ftAdjoint = 0.0;
            for (int j = 0; j < dim; ++j) {
                ftAdjoint += transitionMatrix[j][i] * adjointB[j];
            }
            gradient[i] += adjointB[i] - ftAdjoint;
        }
    }

    private void accumulateRootDiffusionGradient(final double[] gradQ,
                                                 final BranchGradientWorkspace workspace) {
        final int rootIndex = tree.getRoot().getNumber();
        accumulateRootDiffusionGradient(
                preOrder[rootIndex],
                postOrder[rootIndex],
                lastRootDiffusionScale,
                gradQ,
                workspace);
    }

    private void accumulateRootDiffusionGradient(final CanonicalGaussianState rootPreOrder,
                                                 final CanonicalGaussianState rootPostOrder,
                                                 final double rootDiffusionScale,
                                                 final double[] gradQ,
                                                 final BranchGradientWorkspace workspace) {
        CanonicalGaussianMessageOps.combineStates(rootPreOrder, rootPostOrder, workspace.combinedState);

        CanonicalGaussianUtils.fillMomentsFromCanonical(
                workspace.combinedState, workspace.mean, workspace.covariance);
        CanonicalGaussianUtils.fillMomentsFromCanonical(
                rootPreOrder, workspace.mean2, workspace.covariance2);

        final double[] priorPrecision = rootPreOrder.precision;
        for (int i = 0; i < dim; ++i) {
            final double deltaI = workspace.mean[i] - workspace.mean2[i];
            for (int j = 0; j < dim; ++j) {
                workspace.covarianceAdjoint[i][j] =
                        workspace.covariance[i][j] + deltaI * (workspace.mean[j] - workspace.mean2[j]);
            }
        }

        multiplyFlatBySquare(priorPrecision, workspace.covarianceAdjoint, workspace.transitionMatrix, dim);
        multiplySquareByFlat(workspace.transitionMatrix, priorPrecision, workspace.covarianceAdjoint, dim);
        for (int i = 0; i < dim; ++i) {
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                gradQ[iOff + j] += rootDiffusionScale
                        * (-0.5 * priorPrecision[iOff + j] + 0.5 * workspace.covarianceAdjoint[i][j]);
            }
        }
    }

    private static void copyAdjoints(final CanonicalLocalTransitionAdjoints source,
                                     final CanonicalLocalTransitionAdjoints target) {
        System.arraycopy(source.dLogL_dF, 0, target.dLogL_dF, 0, source.dLogL_dF.length);
        System.arraycopy(source.dLogL_df, 0, target.dLogL_df, 0, source.dLogL_df.length);
        System.arraycopy(source.dLogL_dOmega, 0, target.dLogL_dOmega, 0, source.dLogL_dOmega.length);
    }

    private static void transposeFromFlatInto(final double[] source, final double[][] out, final int dim) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                out[j][i] = source[i * dim + j];
            }
        }
    }

    private static void copyFlatToSquare(final double[] source, final double[][] out, final int dim) {
        for (int i = 0; i < dim; ++i) {
            System.arraycopy(source, i * dim, out[i], 0, dim);
        }
    }

    private static void symmetrizeFlatSquare(final double[] matrix, final int dim) {
        for (int i = 0; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                final double avg = 0.5 * (matrix[i * dim + j] + matrix[j * dim + i]);
                matrix[i * dim + j] = avg;
                matrix[j * dim + i] = avg;
            }
        }
    }

    private static void multiplyFlatBySquare(final double[] left, final double[][] right,
                                              final double[][] out, final int dim) {
        for (int i = 0; i < dim; ++i) {
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dim; ++k) {
                    sum += left[iOff + k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void multiplySquareByFlat(final double[][] left, final double[] right,
                                              final double[][] out, final int dim) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dim; ++k) {
                    sum += left[i][k] * right[k * dim + j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void transposeFlatSquareInPlace(final double[] matrix, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final int ij = i * dimension + j;
                final int ji = j * dimension + i;
                final double tmp = matrix[ij];
                matrix[ij] = matrix[ji];
                matrix[ji] = tmp;
            }
        }
    }

    private static void accumulateVectorInPlace(final double[] target,
                                                final double[] source,
                                                final int length) {
        for (int i = 0; i < length; ++i) {
            target[i] += source[i];
        }
    }

    private static void accumulateSquareInPlace(final double[][] target,
                                                final double[][] source) {
        for (int i = 0; i < target.length; ++i) {
            for (int j = 0; j < target[i].length; ++j) {
                target[i][j] += source[i][j];
            }
        }
    }

    private static void symmetrizeSquare(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = i + 1; j < matrix[i].length; ++j) {
                final double average = 0.5 * (matrix[i][j] + matrix[j][i]);
                matrix[i][j] = average;
                matrix[j][i] = average;
            }
        }
    }

    private static double dot(final double[] left, final double[] right) {
        double sum = 0.0;
        for (int i = 0; i < left.length; ++i) {
            sum += left[i] * right[i];
        }
        return sum;
    }

    private static void clearSquare(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            Arrays.fill(matrix[i], 0.0);
        }
    }
}
