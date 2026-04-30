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
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalOUTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalPreparedBranchBasisProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTransitionCacheDiagnostics;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTransitionMomentProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalRootPrior;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTreeMessagePasser;
import dr.evomodel.treedatalikelihood.continuous.framework.MatrixUtils;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.continuous.ou.OrthogonalBlockDiagonalSelectionMatrixParameterization;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianUtils;
import dr.util.TaskPool;

import java.util.Arrays;

/**
 * <p>This class implements the canonical OU tree traversals and exact/partially
 * observed tip elimination in canonical form.
 */
public final class SequentialCanonicalOUMessagePasser implements CanonicalTreeMessagePasser {

    private static final String PHASE_POSTORDER = "postorder";
    private static final String PHASE_PREORDER = "preorder";
    private static final String PHASE_GRADIENT_PREP = "gradientPrep";
    private static final String PHASE_BRANCH_LENGTH_GRADIENT = "branchLengthGradient";
    private static final double BRANCH_LENGTH_FD_RELATIVE_STEP = 1.0e-6;
    private static final double BRANCH_LENGTH_FD_ABSOLUTE_STEP = 1.0e-8;
    private static final String BRANCH_GRADIENT_THREADS_PROPERTY =
            "beast.experimental.canonicalBranchGradientThreads";

    private final Tree tree;
    private final int dim;
    private final int nodeCount;
    private final int tipCount;
    private final CanonicalTreeStateStore stateStore;
    private final CanonicalTipProjector tipProjector;
    private final CanonicalTreeTraversal treeTraversal;
    private final BranchGradientInputs preparedBranchGradientInputs;
    private final BranchGradientWorkspace mainWorkspace;
    private final TaskPool branchGradientTaskPool;
    private final BranchGradientWorkspace[] branchGradientWorkspaces;
    private final CanonicalBranchAdjointPreparer branchAdjointPreparer;
    private final CanonicalTreeGradientEngine treeGradientEngine;

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

        void clear() {
            Arrays.fill(stagedActiveByChild, false);
            activeBranchCount = 0;
            rootDiffusionScale = 0.0;
            hasOrthogonalPreparedBasis = false;
        }

        void checkCompatible(final int expectedCapacity,
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

        void addBranch(final int childIndex,
                       final double branchLength,
                       final CanonicalLocalTransitionAdjoints source,
                       final OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis orthogonalPreparedBasis) {
            if (activeBranchCount >= capacity) {
                throw new IllegalStateException("BranchGradientInputs capacity exceeded.");
            }
            if (orthogonalPreparedBasis != null) {
                hasOrthogonalPreparedBasis = true;
                activeOrthogonalPreparedBasis[activeBranchCount] = orthogonalPreparedBasis;
            }
            activeChildIndices[activeBranchCount] = childIndex;
            activeBranchLengths[activeBranchCount] = branchLength;
            copyAdjoints(source, activeAdjoints[activeBranchCount]);
            activeBranchCount++;
        }

        void stageBranch(final int childIndex,
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

        void clearStagedBranch(final int childIndex) {
            checkChildIndex(childIndex);
            stagedActiveByChild[childIndex] = false;
        }

        void compactStagedBranches(final int rootIndex,
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

        void prepareOrthogonalBasisForActiveBranches(
                final CanonicalPreparedBranchBasisProvider provider) {
            if (provider == null) {
                hasOrthogonalPreparedBasis = false;
                return;
            }

            for (int activeIndex = 0; activeIndex < activeBranchCount; ++activeIndex) {
                final OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis prepared =
                        provider.getOrthogonalPreparedBranchBasis(activeChildIndices[activeIndex]);
                if (prepared == null) {
                    hasOrthogonalPreparedBasis = false;
                    return;
                }
                activeOrthogonalPreparedBasis[activeIndex] = prepared;
            }
            hasOrthogonalPreparedBasis = activeBranchCount > 0;
        }

        OrthogonalBlockDiagonalSelectionMatrixParameterization.PreparedBranchBasis
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

        void setRoot(final double rootDiffusionScale,
                     final CanonicalGaussianState rootPreOrder,
                     final CanonicalGaussianState rootPostOrder) {
            this.rootDiffusionScale = rootDiffusionScale;
            copyState(rootPreOrder, this.rootPreOrder);
            copyState(rootPostOrder, this.rootPostOrder);
        }
    }

    static final class TraversalWorkspace {
        final double[][] traitCovariance;
        final CanonicalGaussianTransition transition;
        final CanonicalGaussianState state;
        final CanonicalGaussianState siblingProduct;
        final CanonicalGaussianState downwardParentState;
        final CanonicalGaussianMessageOps.Workspace gaussianWorkspace;

        TraversalWorkspace(final int dim) {
            this.traitCovariance = new double[dim][dim];
            this.transition = new CanonicalGaussianTransition(dim);
            this.state = new CanonicalGaussianState(dim);
            this.siblingProduct = new CanonicalGaussianState(dim);
            this.downwardParentState = new CanonicalGaussianState(dim);
            this.gaussianWorkspace = new CanonicalGaussianMessageOps.Workspace(dim);
        }
    }

    static final class BranchAdjointWorkspace {
        final CanonicalGaussianState combinedState;
        final CanonicalGaussianState parentPosterior;
        final CanonicalGaussianState pairState;
        final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
        final CanonicalTransitionAdjointUtils.Workspace transitionAdjointWorkspace;
        final CanonicalBranchMessageContribution contribution;
        final CanonicalLocalTransitionAdjoints adjoints;
        final int[] observedIndexScratch;
        final int[] missingIndexScratch;
        final int[] reducedIndexByTraitScratch;
        final double[] mean;
        final double[] mean2;
        final double[][] covariance;
        final double[] varianceFlat;
        final double[] precisionFlat;
        final double[] reducedPrecisionFlatScratch;
        final double[] reducedCovarianceFlatScratch;
        final double[] reducedInformationScratch;
        final double[] reducedMeanScratch;

        BranchAdjointWorkspace(final int dim) {
            this.combinedState = new CanonicalGaussianState(dim);
            this.parentPosterior = new CanonicalGaussianState(dim);
            this.pairState = new CanonicalGaussianState(2 * dim);
            this.contributionWorkspace = new CanonicalBranchMessageContributionUtils.Workspace(dim);
            this.transitionAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(dim);
            this.contribution = new CanonicalBranchMessageContribution(dim);
            this.adjoints = new CanonicalLocalTransitionAdjoints(dim);
            this.observedIndexScratch = new int[dim];
            this.missingIndexScratch = new int[dim];
            this.reducedIndexByTraitScratch = new int[dim];
            this.mean = new double[dim];
            this.mean2 = new double[dim];
            this.covariance = new double[dim][dim];
            this.varianceFlat = new double[dim * dim];
            this.precisionFlat = new double[dim * dim];
            this.reducedPrecisionFlatScratch = new double[4 * dim * dim];
            this.reducedCovarianceFlatScratch = new double[4 * dim * dim];
            this.reducedInformationScratch = new double[2 * dim];
            this.reducedMeanScratch = new double[2 * dim];
        }
    }

    static final class GradientPullbackWorkspace {
        final double[] orthogonalStationaryMeanScratch;
        final double[] orthogonalCompressedGradientScratch;
        final double[] orthogonalNativeGradientScratch;
        final double[] orthogonalRotationGradientFlatScratch;
        final double[][] orthogonalRotationGradientScratch;
        final double[][] covariance2;
        final double[] transitionMatrixFlat;
        final double[][] transitionMatrix;
        final double[][] covarianceAdjoint;
        final double[] localGradientA;
        final double[] localGradientQ;
        final double[] localGradientMuVector;
        final double[] localGradientMuScalar;
        OrthogonalBlockDiagonalSelectionMatrixParameterization.BranchGradientWorkspace orthogonalBranchWorkspace;

        GradientPullbackWorkspace(final int dim) {
            this.orthogonalStationaryMeanScratch = new double[dim];
            this.orthogonalCompressedGradientScratch = new double[dim + 2 * (dim / 2)];
            this.orthogonalNativeGradientScratch = new double[((dim & 1) == 1 ? 1 : 0) + 3 * (dim / 2)];
            this.orthogonalRotationGradientFlatScratch = new double[dim * dim];
            this.orthogonalRotationGradientScratch = new double[dim][dim];
            this.covariance2 = new double[dim][dim];
            this.transitionMatrixFlat = new double[dim * dim];
            this.transitionMatrix = new double[dim][dim];
            this.covarianceAdjoint = new double[dim][dim];
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
                Arrays.fill(orthogonalRotationGradientFlatScratch, 0.0);
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

    static final class BranchGradientWorkspace {
        final TraversalWorkspace traversal;
        final BranchAdjointWorkspace adjoint;
        final GradientPullbackWorkspace gradient;

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
        final double[] orthogonalRotationGradientFlatScratch;
        final double[][] orthogonalRotationGradientScratch;
        final double[][] covariance;
        final double[][] covariance2;
        final double[] transitionMatrixFlat;
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

        BranchGradientWorkspace(final int dim) {
            this.traversal = new TraversalWorkspace(dim);
            this.adjoint = new BranchAdjointWorkspace(dim);
            this.gradient = new GradientPullbackWorkspace(dim);

            this.traitCovariance = traversal.traitCovariance;
            this.transition = traversal.transition;
            this.state = traversal.state;
            this.siblingProduct = traversal.siblingProduct;
            this.downwardParentState = traversal.downwardParentState;
            this.gaussianWorkspace = traversal.gaussianWorkspace;

            this.combinedState = adjoint.combinedState;
            this.parentPosterior = adjoint.parentPosterior;
            this.pairState = adjoint.pairState;
            this.contributionWorkspace = adjoint.contributionWorkspace;
            this.transitionAdjointWorkspace = adjoint.transitionAdjointWorkspace;
            this.contribution = adjoint.contribution;
            this.adjoints = adjoint.adjoints;
            this.observedIndexScratch = adjoint.observedIndexScratch;
            this.missingIndexScratch = adjoint.missingIndexScratch;
            this.reducedIndexByTraitScratch = adjoint.reducedIndexByTraitScratch;
            this.mean = adjoint.mean;
            this.mean2 = adjoint.mean2;
            this.covariance = adjoint.covariance;
            this.varianceFlat = adjoint.varianceFlat;
            this.precisionFlat = adjoint.precisionFlat;
            this.reducedPrecisionFlatScratch = adjoint.reducedPrecisionFlatScratch;
            this.reducedCovarianceFlatScratch = adjoint.reducedCovarianceFlatScratch;
            this.reducedInformationScratch = adjoint.reducedInformationScratch;
            this.reducedMeanScratch = adjoint.reducedMeanScratch;

            this.orthogonalStationaryMeanScratch = gradient.orthogonalStationaryMeanScratch;
            this.orthogonalCompressedGradientScratch = gradient.orthogonalCompressedGradientScratch;
            this.orthogonalNativeGradientScratch = gradient.orthogonalNativeGradientScratch;
            this.orthogonalRotationGradientFlatScratch = gradient.orthogonalRotationGradientFlatScratch;
            this.orthogonalRotationGradientScratch = gradient.orthogonalRotationGradientScratch;
            this.covariance2 = gradient.covariance2;
            this.transitionMatrixFlat = gradient.transitionMatrixFlat;
            this.transitionMatrix = gradient.transitionMatrix;
            this.covarianceAdjoint = gradient.covarianceAdjoint;
            this.localGradientA = gradient.localGradientA;
            this.localGradientQ = gradient.localGradientQ;
            this.localGradientMuVector = gradient.localGradientMuVector;
            this.localGradientMuScalar = gradient.localGradientMuScalar;
        }

        double[] localGradientMu(final int gradientLength, final int dim) {
            return gradient.localGradientMu(gradientLength, dim);
        }

        void clearLocalGradientBuffers(final int gradALength,
                                       final int gradMuLength,
                                       final int dim,
                                       final boolean orthogonalSelection,
                                       final int compressedGradientLength) {
            gradient.clearLocalGradientBuffers(
                    gradALength, gradMuLength, dim, orthogonalSelection, compressedGradientLength);
        }

        OrthogonalBlockDiagonalSelectionMatrixParameterization.BranchGradientWorkspace
        ensureOrthogonalBranchWorkspace(final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalSelection) {
            return gradient.ensureOrthogonalBranchWorkspace(orthogonalSelection);
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
        this.stateStore = new CanonicalTreeStateStore(nodeCount, dim);
        this.tipProjector = new CanonicalTipProjector(dim);
        this.treeTraversal = new CanonicalTreeTraversal(tree, dim, tipProjector);
        this.preparedBranchGradientInputs = new BranchGradientInputs(Math.max(0, nodeCount - 1), dim);
        this.mainWorkspace = new BranchGradientWorkspace(dim);
        this.branchGradientTaskPool = new TaskPool(nodeCount, Math.max(1, branchGradientParallelism));
        final int branchGradientWorkspaceCount =
                branchGradientTaskPool.getNumThreads() <= 1 ? 1 : branchGradientTaskPool.getNumThreads() + 1;
        this.branchGradientWorkspaces = new BranchGradientWorkspace[branchGradientWorkspaceCount];
        for (int i = 0; i < branchGradientWorkspaces.length; ++i) {
            branchGradientWorkspaces[i] = new BranchGradientWorkspace(dim);
        }
        this.branchAdjointPreparer = new CanonicalBranchAdjointPreparer(
                tree,
                dim,
                branchGradientTaskPool,
                mainWorkspace,
                branchGradientWorkspaces,
                this::fillLocalAdjointsForBranch);
        this.treeGradientEngine = new CanonicalTreeGradientEngine(
                dim,
                Math.max(0, nodeCount - 1),
                branchGradientTaskPool,
                mainWorkspace,
                branchGradientWorkspaces);
    }

    private static int resolveBranchGradientParallelism() {
        final String propertyValue = System.getProperty(BRANCH_GRADIENT_THREADS_PROPERTY);
        if (propertyValue != null) {
            try {
                return Math.max(1, Integer.parseInt(propertyValue));
            } catch (NumberFormatException ignored) {
                // Fall back to the reproducible default below.
            }
        }
        return 1;
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
        stateStore.tipObservations[tipIndex].copyFrom(observation);
    }

    @Override
    public double computePostOrderLogLikelihood(final CanonicalBranchTransitionProvider transitionProvider,
                                                final CanonicalRootPrior rootPrior) {
        final String previousPhase = pushTransitionCachePhase(transitionProvider, PHASE_POSTORDER);
        try {
            return treeTraversal.computePostOrderLogLikelihood(
                    transitionProvider,
                    rootPrior,
                    stateStore,
                    mainWorkspace);
        } finally {
            popTransitionCachePhase(transitionProvider, previousPhase);
        }
    }

    @Override
    public CanonicalGaussianState getPostOrderState(final int nodeIndex) {
        return stateStore.postOrder[nodeIndex];
    }

    @Override
    public void computePreOrder(final CanonicalBranchTransitionProvider transitionProvider,
                                final CanonicalRootPrior rootPrior) {
        final String previousPhase = pushTransitionCachePhase(transitionProvider, PHASE_PREORDER);
        try {
            treeTraversal.computePreOrder(
                    transitionProvider,
                    rootPrior,
                    stateStore,
                    mainWorkspace);
        } finally {
            popTransitionCachePhase(transitionProvider, previousPhase);
        }
    }

    @Override
    public CanonicalGaussianState getPreOrderState(final int nodeIndex) {
        return stateStore.preOrder[nodeIndex];
    }

    @Override
    @Deprecated
    public void computeGradientQ(final CanonicalBranchTransitionProvider transitionProvider, final double[] gradQ) {
        Arrays.fill(gradQ, 0.0);
        ensureGradientState();

        final String previousPhase = pushTransitionCachePhase(transitionProvider, PHASE_GRADIENT_PREP);
        try {
            final CanonicalOUTransitionProvider ouProvider = requireOUProvider(transitionProvider);
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

            if (stateStore.lastRootDiffusionScale > 0.0) {
                accumulateRootDiffusionGradient(gradQ, workspace);
            }
        } finally {
            popTransitionCachePhase(transitionProvider, previousPhase);
        }
    }

    @Override
    public void computeGradientBranchLengths(final CanonicalBranchTransitionProvider transitionProvider,
                                             final double[] gradT) {
        Arrays.fill(gradT, 0.0);
        ensureGradientState();

        final String previousPhase = pushTransitionCachePhase(transitionProvider, PHASE_BRANCH_LENGTH_GRADIENT);
        try {
            final CanonicalOUTransitionProvider ouProvider = requireOUProvider(transitionProvider);
            final int rootIndex = tree.getRoot().getNumber();
            for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
                if (childIndex == rootIndex) {
                    continue;
                }
                final CanonicalTipObservation tipObservation =
                        tree.isExternal(tree.getNode(childIndex)) ? stateStore.tipObservations[childIndex] : null;
                if (tipObservation != null && tipObservation.observedCount == 0) {
                    gradT[childIndex] = 0.0;
                    continue;
                }
                final double branchLength = transitionProvider.getEffectiveBranchLength(childIndex);
                gradT[childIndex] = finiteDifferenceBranchLengthGradient(childIndex, ouProvider, branchLength);
            }
        } finally {
            popTransitionCachePhase(transitionProvider, previousPhase);
        }
    }

    @Override
    @Deprecated
    public void computeGradientA(final CanonicalBranchTransitionProvider transitionProvider, final double[] gradA) {
        Arrays.fill(gradA, 0.0);
        ensureGradientState();

        final CanonicalOUTransitionProvider ouProvider = requireOUProvider(transitionProvider);
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
    @Deprecated
    public void computeGradientMu(final CanonicalBranchTransitionProvider transitionProvider, final double[] gradMu) {
        Arrays.fill(gradMu, 0.0);
        ensureGradientState();

        final CanonicalOUTransitionProvider ouProvider = requireOUProvider(transitionProvider);
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
        final String previousPhase = pushTransitionCachePhase(transitionProvider, PHASE_GRADIENT_PREP);
        try {
            ensureGradientState();
            branchAdjointPreparer.prepare(transitionProvider, stateStore, out);
        } finally {
            popTransitionCachePhase(transitionProvider, previousPhase);
        }
    }

    public void computeJointGradients(final CanonicalBranchTransitionProvider transitionProvider,
                                      final BranchGradientInputs inputs,
                                      final double[] gradA,
                                      final double[] gradQ,
                                      final double[] gradMu) {
        treeGradientEngine.computeJointGradients(transitionProvider, inputs, gradA, gradQ, gradMu);
    }

    @Override
    public void storeState() {
        stateStore.storeState();
    }

    @Override
    public void restoreState() {
        stateStore.restoreState();
    }

    @Override
    public void acceptState() { }

    private static void clearState(final CanonicalGaussianState state) {
        CanonicalGaussianMessageOps.clearState(state);
    }

    private static void copyState(final CanonicalGaussianState source, final CanonicalGaussianState target) {
        CanonicalGaussianMessageOps.copyState(source, target);
    }

    private void buildUpwardParentMessageForTransition(final int childIndex,
                                                       final CanonicalGaussianTransition transition,
                                                       final CanonicalGaussianState out,
                                                       final BranchGradientWorkspace workspace) {
        if (tree.isExternal(tree.getNode(childIndex))) {
            final CanonicalTipObservation tipObservation = stateStore.tipObservations[childIndex];
            if (tipObservation.observedCount == 0) {
                clearState(out);
            } else if (tipObservation.observedCount == dim) {
                CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transition, tipObservation.values, out);
            } else {
                throw new IllegalStateException(
                        "Partially observed canonical tips must use the missing-mask branch update.");
            }
        } else {
            CanonicalGaussianMessageOps.pushBackward(stateStore.postOrder[childIndex], transition, workspace.gaussianWorkspace, out);
        }
    }

    private void buildTipParentMessageWithMissingMask(final CanonicalTipObservation tipObservation,
                                                      final CanonicalTransitionMomentProvider transitionMomentProvider,
                                                      final double branchLength,
                                                      final CanonicalGaussianState out,
                                                      final BranchGradientWorkspace workspace) {
        tipProjector.projectObservedChildToParent(tipObservation, transitionMomentProvider, branchLength, out);
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
        if (!stateStore.hasPostOrderState || !stateStore.hasPreOrderState) {
            throw new IllegalStateException(
                    "Canonical gradients require both computePostOrderLogLikelihood and computePreOrder to have been called.");
        }
    }

    private CanonicalOUTransitionProvider requireOUProvider(
            final CanonicalBranchTransitionProvider transitionProvider) {
        if (!(transitionProvider instanceof CanonicalOUTransitionProvider)) {
            throw new UnsupportedOperationException(
                    "Not yet implemented: canonical OU gradients currently support only "
                            + "CanonicalOUTransitionProvider implementations.");
        }
        return (CanonicalOUTransitionProvider) transitionProvider;
    }

    private String pushTransitionCachePhase(final CanonicalBranchTransitionProvider transitionProvider,
                                            final String phase) {
        if (transitionProvider instanceof CanonicalTransitionCacheDiagnostics) {
            return ((CanonicalTransitionCacheDiagnostics) transitionProvider).pushDiagnosticPhase(phase);
        }
        return null;
    }

    private void popTransitionCachePhase(final CanonicalBranchTransitionProvider transitionProvider,
                                         final String previousPhase) {
        if (transitionProvider instanceof CanonicalTransitionCacheDiagnostics) {
            ((CanonicalTransitionCacheDiagnostics) transitionProvider).popDiagnosticPhase(previousPhase);
        }
    }

    private boolean fillLocalAdjointsForBranch(final int childIndex,
                                               final CanonicalBranchTransitionProvider transitionProvider,
                                               final BranchGradientWorkspace workspace) {
        transitionProvider.fillCanonicalTransition(childIndex, workspace.transition);
        if (stateStore.hasFixedRootValue && tree.isRoot(tree.getParent(tree.getNode(childIndex)))) {
            if (tree.isExternal(tree.getNode(childIndex))) {
                final CanonicalTipObservation tipObservation = stateStore.tipObservations[childIndex];
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
                fillContributionForFixedParentInternalNode(workspace.transition, stateStore.postOrder[childIndex], workspace);
            }

            CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                    workspace.transition,
                    workspace.contribution,
                    workspace.transitionAdjointWorkspace,
                    workspace.adjoints);
            return true;
        }

        if (tree.isExternal(tree.getNode(childIndex))) {
            final CanonicalTipObservation tipObservation = stateStore.tipObservations[childIndex];
            if (tipObservation.observedCount == 0) {
                return false;
            }
            final int observedCount = collectObservationPartition(tipObservation, workspace);
            if (observedCount == dim) {
                fillContributionForObservedTip(stateStore.branchAboveParent[childIndex], workspace.transition, tipObservation, workspace);
            } else {
                fillContributionForPartiallyObservedTip(
                        stateStore.branchAboveParent[childIndex], workspace.transition, tipObservation, observedCount, workspace);
            }
        } else {
            CanonicalGaussianMessageOps.buildPairPosterior(
                    stateStore.branchAboveParent[childIndex],
                    workspace.transition,
                    stateStore.postOrder[childIndex],
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
            final double xi = stateStore.fixedRootValue[i];
            final double yi = tipObservation.values[i];
            workspace.contribution.dLogL_dInformationX[i] = xi;
            workspace.contribution.dLogL_dInformationY[i] = yi;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = stateStore.fixedRootValue[j];
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
        CanonicalGaussianMessageOps.conditionOnObservedFirstBlock(transition, stateStore.fixedRootValue, workspace.state);
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
        CanonicalGaussianMessageOps.conditionOnObservedFirstBlock(transition, stateStore.fixedRootValue, workspace.state);

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
            workspace.contribution.dLogL_dInformationX[i] = stateStore.fixedRootValue[i];
            workspace.contribution.dLogL_dInformationY[i] = tipObservation.observed[i]
                    ? tipObservation.values[i]
                    : workspace.mean2[workspace.reducedIndexByTraitScratch[i] - dim];
        }

        for (int i = 0; i < dim; ++i) {
            final double xi = stateStore.fixedRootValue[i];
            final double yi = workspace.contribution.dLogL_dInformationY[i];
            final int missingI = tipObservation.observed[i] ? -1 : workspace.reducedIndexByTraitScratch[i] - dim;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = stateStore.fixedRootValue[j];
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
            final double xi = stateStore.fixedRootValue[i];
            final double yi = childMean[i];
            workspace.contribution.dLogL_dInformationX[i] = xi;
            workspace.contribution.dLogL_dInformationY[i] = yi;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = stateStore.fixedRootValue[j];
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
            final CanonicalOUTransitionProvider provider,
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
                                                final CanonicalOUTransitionProvider provider,
                                                final double branchLength) {
        final BranchGradientWorkspace workspace = mainWorkspace;
        final CanonicalTipObservation tipObservation =
                tree.isExternal(tree.getNode(childIndex)) ? stateStore.tipObservations[childIndex] : null;
        if (tipObservation != null && tipObservation.observedCount < dim) {
            buildTipParentMessageWithMissingMask(
                    tipObservation,
                    provider,
                    branchLength,
                    workspace.state,
                    workspace);
        } else {
            provider.fillCanonicalTransitionForLength(branchLength, workspace.transition);
            buildUpwardParentMessageForTransition(childIndex, workspace.transition, workspace.state, workspace);
        }
        CanonicalGaussianMessageOps.combineStates(
                stateStore.branchAboveParent[childIndex], workspace.state, workspace.combinedState);
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
                stateStore.preOrder[rootIndex],
                stateStore.postOrder[rootIndex],
                stateStore.lastRootDiffusionScale,
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

    private static void clearSquare(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            Arrays.fill(matrix[i], 0.0);
        }
    }
}
