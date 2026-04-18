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

import java.util.Arrays;

/**
 * <p>This class implements the canonical OU tree traversals and exact/partially
 * observed tip elimination in canonical form.
 */
public final class SequentialCanonicalOUMessagePasser implements CanonicalTreeMessagePasser {

    private static final double BRANCH_LENGTH_FD_RELATIVE_STEP = 1.0e-6;
    private static final double BRANCH_LENGTH_FD_ABSOLUTE_STEP = 1.0e-8;
    private static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);

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

    private final double[][] scratchTraitCovariance;
    private final CanonicalGaussianTransition scratchTransition;
    private final CanonicalGaussianState scratchState;
    private final CanonicalGaussianState scratchSiblingProduct;
    private final CanonicalGaussianState scratchDownwardParentState;
    private final CanonicalGaussianState scratchCombinedState;
    private final CanonicalGaussianState scratchParentPosterior;
    private final CanonicalGaussianState scratchPairState;
    private final CanonicalGaussianMessageOps.Workspace workspace;
    private final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
    private final CanonicalTransitionAdjointUtils.Workspace transitionAdjointWorkspace;
    private final CanonicalBranchMessageContribution scratchContribution;
    private final CanonicalLocalTransitionAdjoints scratchAdjoints;
    private final int[] observedIndexScratch;
    private final int[] missingIndexScratch;
    private final int[] reducedIndexByTraitScratch;
    private final double[] scratchMean;
    private final double[] scratchMean2;
    private final double[] orthogonalStationaryMeanScratch;
    private final double[] orthogonalCompressedGradientScratch;
    private final double[] orthogonalNativeGradientScratch;
    private final double[][] orthogonalRotationGradientScratch;
    private double[] fixedRootValue;
    private double[] storedFixedRootValue;
    private final double[][] scratchCovariance;
    private final double[][] scratchCovariance2;
    private final double[][] scratchTransitionMatrix;
    private final double[][] scratchCovarianceAdjoint;
    private final double[] scratchVarianceFlat;
    private final double[] scratchPrecisionFlat;
    private final double[] reducedPrecisionFlatScratch;
    private final double[] reducedCovarianceFlatScratch;
    private final double[] reducedInformationScratch;
    private final double[] reducedMeanScratch;

    public SequentialCanonicalOUMessagePasser(final Tree tree, final int dim) {
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
        this.scratchTraitCovariance = new double[dim][dim];
        this.scratchTransition = new CanonicalGaussianTransition(dim);
        this.scratchState = new CanonicalGaussianState(dim);
        this.scratchSiblingProduct = new CanonicalGaussianState(dim);
        this.scratchDownwardParentState = new CanonicalGaussianState(dim);
        this.scratchCombinedState = new CanonicalGaussianState(dim);
        this.scratchParentPosterior = new CanonicalGaussianState(dim);
        this.scratchPairState = new CanonicalGaussianState(2 * dim);
        this.workspace = new CanonicalGaussianMessageOps.Workspace(dim);
        this.contributionWorkspace = new CanonicalBranchMessageContributionUtils.Workspace(dim);
        this.transitionAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(dim);
        this.scratchContribution = new CanonicalBranchMessageContribution(dim);
        this.scratchAdjoints = new CanonicalLocalTransitionAdjoints(dim);
        this.observedIndexScratch = new int[dim];
        this.missingIndexScratch = new int[dim];
        this.reducedIndexByTraitScratch = new int[dim];
        this.scratchMean = new double[dim];
        this.scratchMean2 = new double[dim];
        this.orthogonalStationaryMeanScratch = new double[dim];
        this.orthogonalCompressedGradientScratch = new double[dim + 2 * (dim / 2)];
        this.orthogonalNativeGradientScratch = new double[((dim & 1) == 1 ? 1 : 0) + 3 * (dim / 2)];
        this.orthogonalRotationGradientScratch = new double[dim][dim];
        this.fixedRootValue = new double[dim];
        this.storedFixedRootValue = new double[dim];
        this.scratchCovariance = new double[dim][dim];
        this.scratchCovariance2 = new double[dim][dim];
        this.scratchTransitionMatrix = new double[dim][dim];
        this.scratchCovarianceAdjoint = new double[dim][dim];
        this.scratchVarianceFlat = new double[dim * dim];
        this.scratchPrecisionFlat = new double[dim * dim];
        this.reducedPrecisionFlatScratch = new double[4 * dim * dim];
        this.reducedCovarianceFlatScratch = new double[4 * dim * dim];
        this.reducedInformationScratch = new double[2 * dim];
        this.reducedMeanScratch = new double[2 * dim];
        this.hasPostOrderState = false;
        this.storedHasPostOrderState = false;
        this.hasPreOrderState = false;
        this.storedHasPreOrderState = false;
        this.lastRootDiffusionScale = 0.0;
        this.storedLastRootDiffusionScale = 0.0;
        this.hasFixedRootValue = false;
        this.storedHasFixedRootValue = false;
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
                computePostOrderAtInternalNode(i, transitionProvider);
            }
        }

        final int rootIndex = tree.getRoot().getNumber();
        transitionProvider.fillTraitCovariance(scratchTraitCovariance);
        hasPostOrderState = true;
        return rootPrior.computeLogMarginalLikelihood(postOrder[rootIndex], scratchTraitCovariance);
    }

    @Override
    public CanonicalGaussianState getPostOrderState(final int nodeIndex) {
        return postOrder[nodeIndex];
    }

    @Override
    public void computePreOrder(final CanonicalBranchTransitionProvider transitionProvider,
                                final CanonicalRootPrior rootPrior) {
        final int rootIndex = tree.getRoot().getNumber();
        transitionProvider.fillTraitCovariance(scratchTraitCovariance);
        CanonicalGaussianMessageOps.clearState(preOrder[rootIndex]);
        hasFixedRootValue = rootPrior.isFixedRoot();
        if (hasFixedRootValue) {
            rootPrior.fillFixedRootValue(fixedRootValue);
            lastRootDiffusionScale = 0.0;
        } else {
            rootPrior.fillRootPriorState(scratchTraitCovariance, preOrder[rootIndex]);
            lastRootDiffusionScale = rootPrior.getDiffusionScale();
        }
        clearState(branchAboveParent[rootIndex]);
        computePreOrderRecursive(tree.getRoot().getNumber(), transitionProvider);
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

        for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
            if (childIndex == rootIndex) {
                continue;
            }
            if (!fillLocalAdjointsForBranch(childIndex, transitionProvider)) {
                continue;
            }
            if (orthogonalSelection != null) {
                transposeFromFlatInto(scratchAdjoints.dLogL_dOmega, scratchCovarianceAdjoint, dim);
                orthogonalSelection.accumulateDiffusionGradient(
                        processModel.getDiffusionMatrix(),
                        transitionProvider.getEffectiveBranchLength(childIndex),
                        scratchCovarianceAdjoint,
                        gradQ);
            } else {
                copyFlatToSquare(scratchAdjoints.dLogL_dOmega, scratchCovarianceAdjoint, dim);
                processModel.accumulateDiffusionGradient(
                        transitionProvider.getEffectiveBranchLength(childIndex),
                        scratchCovarianceAdjoint,
                        gradQ);
            }
        }

        if (lastRootDiffusionScale > 0.0) {
            accumulateRootDiffusionGradient(gradQ);
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

        for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
            if (childIndex == rootIndex) {
                continue;
            }
            if (!fillLocalAdjointsForBranch(childIndex, transitionProvider)) {
                continue;
            }

            final double branchLength = transitionProvider.getEffectiveBranchLength(childIndex);
            copyFlatToSquare(scratchAdjoints.dLogL_dF, scratchTransitionMatrix, dim);
            processModel.accumulateSelectionGradient(
                    branchLength,
                    scratchTransitionMatrix,
                    scratchAdjoints.dLogL_df,
                    gradA);

            transposeFromFlatInto(scratchAdjoints.dLogL_dOmega, scratchCovarianceAdjoint, dim);
            processModel.accumulateSelectionGradientFromCovariance(
                    branchLength,
                    scratchCovarianceAdjoint,
                    gradA);
        }

        transposeFlatSquareInPlace(gradA, dim);
    }

    private void computeOrthogonalBlockGradientA(final CanonicalBranchTransitionProvider transitionProvider,
                                                 final OUProcessModel processModel,
                                                 final OrthogonalBlockDiagonalSelectionMatrixParameterization parameterization,
                                                 final double[] gradA) {
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
        if (compressedBlockDim > orthogonalCompressedGradientScratch.length
                || nativeBlockDim > orthogonalNativeGradientScratch.length) {
            throw new IllegalStateException(
                    "Orthogonal block scratch is too small for native gradient dimensions "
                            + compressedBlockDim + " and " + nativeBlockDim + ".");
        }
        final double[] stationaryMean = orthogonalStationaryMeanScratch;
        processModel.getInitialMean(stationaryMean);
        final double[] compressedBlockGradient = orthogonalCompressedGradientScratch;
        Arrays.fill(compressedBlockGradient, 0, compressedBlockDim, 0.0);
        final double[] nativeBlockGradient = orthogonalNativeGradientScratch;
        final double[][] rotationGradient = orthogonalRotationGradientScratch;
        clearSquare(rotationGradient);

        for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
            if (childIndex == rootIndex) {
                continue;
            }
            if (!fillLocalAdjointsForBranch(childIndex, transitionProvider)) {
                continue;
            }

            parameterization.accumulateNativeGradientFromAdjoints(
                    processModel.getDiffusionMatrix(),
                    stationaryMean,
                    transitionProvider.getEffectiveBranchLength(childIndex),
                    scratchAdjoints,
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

        for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
            if (childIndex == rootIndex) {
                continue;
            }
            if (!fillLocalAdjointsForBranch(childIndex, transitionProvider)) {
                continue;
            }

            final double branchLength = transitionProvider.getEffectiveBranchLength(childIndex);
            processModel.fillTransitionMatrix(branchLength, scratchTransitionMatrix);
            accumulateStationaryMeanGradient(
                    scratchTransitionMatrix,
                    scratchAdjoints.dLogL_df,
                    gradMu);
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
                                                final CanonicalBranchTransitionProvider transitionProvider) {
        final CanonicalGaussianState dest = postOrder[nodeIndex];
        clearState(dest);

        final int childCount = tree.getChildCount(tree.getNode(nodeIndex));
        boolean first = true;
        for (int c = 0; c < childCount; c++) {
            final int childIndex = tree.getChild(tree.getNode(nodeIndex), c).getNumber();
            buildUpwardParentMessage(childIndex, transitionProvider, scratchState);
            if (first) {
                copyState(scratchState, dest);
                first = false;
            } else {
                CanonicalGaussianMessageOps.combineStateInPlace(dest, scratchState);
            }
        }
    }

    private void computePreOrderRecursive(final int parentIndex,
                                          final CanonicalBranchTransitionProvider transitionProvider) {
        if (tree.isExternal(tree.getNode(parentIndex))) {
            return;
        }

        final int rootIndex = tree.getRoot().getNumber();
        final int childCount = tree.getChildCount(tree.getNode(parentIndex));
        for (int c = 0; c < childCount; c++) {
            final int childIndex = tree.getChild(tree.getNode(parentIndex), c).getNumber();
            if (hasFixedRootValue && parentIndex == rootIndex) {
                clearState(branchAboveParent[childIndex]);
                transitionProvider.fillCanonicalTransition(childIndex, scratchTransition);
                CanonicalGaussianMessageOps.conditionOnObservedFirstBlock(
                        scratchTransition,
                        fixedRootValue,
                        preOrder[childIndex]);
                computePreOrderRecursive(childIndex, transitionProvider);
                continue;
            }

            final boolean hasSiblings = buildSiblingProduct(parentIndex, childIndex, transitionProvider, scratchSiblingProduct);

            if (hasSiblings) {
                CanonicalGaussianMessageOps.combineStates(preOrder[parentIndex], scratchSiblingProduct, scratchDownwardParentState);
            } else {
                copyState(preOrder[parentIndex], scratchDownwardParentState);
            }

            copyState(scratchDownwardParentState, branchAboveParent[childIndex]);
            transitionProvider.fillCanonicalTransition(childIndex, scratchTransition);
            CanonicalGaussianMessageOps.pushForward(
                    scratchDownwardParentState, scratchTransition, workspace, preOrder[childIndex]);
            computePreOrderRecursive(childIndex, transitionProvider);
        }
    }

    private boolean buildSiblingProduct(final int parentIndex,
                                        final int excludedChildIndex,
                                        final CanonicalBranchTransitionProvider transitionProvider,
                                        final CanonicalGaussianState out) {
        clearState(out);
        final int childCount = tree.getChildCount(tree.getNode(parentIndex));
        boolean found = false;
        for (int c = 0; c < childCount; c++) {
            final int siblingIndex = tree.getChild(tree.getNode(parentIndex), c).getNumber();
            if (siblingIndex == excludedChildIndex) {
                continue;
            }
            buildUpwardParentMessage(siblingIndex, transitionProvider, scratchState);
            if (!found) {
                copyState(scratchState, out);
                found = true;
            } else {
                CanonicalGaussianMessageOps.combineStateInPlace(out, scratchState);
            }
        }
        return found;
    }

    private void buildUpwardParentMessage(final int childIndex,
                                          final CanonicalBranchTransitionProvider transitionProvider,
                                          final CanonicalGaussianState out) {
        if (tree.isExternal(tree.getNode(childIndex))) {
            final CanonicalTipObservation tipObservation = tipObservations[childIndex];
            if (tipObservation.observedCount < dim) {
                buildTipParentMessage(childIndex, tipObservation, transitionProvider, out);
                return;
            }
        }
        transitionProvider.fillCanonicalTransition(childIndex, scratchTransition);
        buildUpwardParentMessageForTransition(childIndex, scratchTransition, out);
    }

    private void buildUpwardParentMessageForTransition(final int childIndex,
                                                       final CanonicalGaussianTransition transition,
                                                       final CanonicalGaussianState out) {
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
            CanonicalGaussianMessageOps.pushBackward(postOrder[childIndex], transition, workspace, out);
        }
    }

    private void buildTipParentMessage(final int childIndex,
                                       final CanonicalTipObservation tipObservation,
                                       final CanonicalBranchTransitionProvider transitionProvider,
                                       final CanonicalGaussianState out) {
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
                out);
    }

    private void buildTipParentMessageWithMissingMask(final CanonicalTipObservation tipObservation,
                                                      final OUProcessModel processModel,
                                                      final double branchLength,
                                                      final CanonicalGaussianState out) {
        final int observedCount = collectObservationPartition(tipObservation);
        if (observedCount == 0) {
            clearState(out);
            return;
        }

        processModel.fillTransitionMatrix(branchLength, scratchTransitionMatrix);
        processModel.fillTransitionOffset(branchLength, scratchMean);
        processModel.fillTransitionCovariance(branchLength, scratchCovariance);

        for (int observed = 0; observed < observedCount; ++observed) {
            final int observedTrait = observedIndexScratch[observed];
            scratchMean2[observed] = tipObservation.values[observedTrait] - scratchMean[observedTrait];
            final int rowOffset = observed * observedCount;
            for (int otherObserved = 0; otherObserved < observedCount; ++otherObserved) {
                scratchVarianceFlat[rowOffset + otherObserved] =
                        scratchCovariance[observedTrait][observedIndexScratch[otherObserved]];
            }
        }

        final double logDetVariance = MatrixUtils.invertSymmetricPositiveDefiniteCompact(
                scratchVarianceFlat,
                scratchPrecisionFlat,
                observedCount,
                scratchMean,
                reducedMeanScratch);

        for (int row = 0; row < observedCount; ++row) {
            double sum = 0.0;
            final int rowOffset = row * observedCount;
            for (int k = 0; k < observedCount; ++k) {
                sum += scratchPrecisionFlat[rowOffset + k] * scratchMean2[k];
            }
            scratchMean[row] = sum;
        }

        for (int row = 0; row < observedCount; ++row) {
            final int observedTrait = observedIndexScratch[row];
            final int precisionRowOffset = row * observedCount;
            final int projectedRowOffset = row * dim;
            for (int col = 0; col < dim; ++col) {
                double sum = 0.0;
                for (int k = 0; k < observedCount; ++k) {
                    sum += scratchPrecisionFlat[precisionRowOffset + k]
                            * scratchTransitionMatrix[observedIndexScratch[k]][col];
                }
                reducedPrecisionFlatScratch[projectedRowOffset + col] = sum;
            }
        }

        for (int i = 0; i < dim; ++i) {
            double information = 0.0;
            for (int observed = 0; observed < observedCount; ++observed) {
                information += scratchTransitionMatrix[observedIndexScratch[observed]][i] * scratchMean[observed];
            }
            out.information[i] = information;

            for (int j = 0; j < dim; ++j) {
                double precision = 0.0;
                for (int observed = 0; observed < observedCount; ++observed) {
                    precision += scratchTransitionMatrix[observedIndexScratch[observed]][i]
                            * reducedPrecisionFlatScratch[observed * dim + j];
                }
                out.precision[i * dim + j] = precision;
            }
        }
        symmetrizeFlatSquare(out.precision, dim);

        double quadratic = 0.0;
        for (int observed = 0; observed < observedCount; ++observed) {
            quadratic += scratchMean2[observed] * scratchMean[observed];
        }
        out.logNormalizer = 0.5 * (
                observedCount * LOG_TWO_PI
                        + logDetVariance
                        + quadratic);
    }

    private int collectObservationPartition(final CanonicalTipObservation tipObservation) {
        int observedCount = 0;
        int missingCount = 0;

        for (int i = 0; i < dim; i++) {
            if (tipObservation.observed[i]) {
                observedIndexScratch[observedCount++] = i;
                reducedIndexByTraitScratch[i] = -1;
            } else {
                missingIndexScratch[missingCount++] = i;
                reducedIndexByTraitScratch[i] = dim + missingCount - 1;
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
                                               final CanonicalBranchTransitionProvider transitionProvider) {
        transitionProvider.fillCanonicalTransition(childIndex, scratchTransition);
        if (hasFixedRootValue && tree.isRoot(tree.getParent(tree.getNode(childIndex)))) {
            if (tree.isExternal(tree.getNode(childIndex))) {
                final CanonicalTipObservation tipObservation = tipObservations[childIndex];
                if (tipObservation.observedCount == 0) {
                    return false;
                }
                final int observedCount = collectObservationPartition(tipObservation);
                if (observedCount == dim) {
                    fillContributionForFixedParentObservedTip(scratchTransition, tipObservation);
                } else {
                    fillContributionForFixedParentPartiallyObservedTip(
                            scratchTransition, tipObservation, observedCount);
                }
            } else {
                fillContributionForFixedParentInternalNode(scratchTransition, postOrder[childIndex]);
            }

            CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                    scratchTransition,
                    scratchContribution,
                    transitionAdjointWorkspace,
                    scratchAdjoints);
            return true;
        }

        if (tree.isExternal(tree.getNode(childIndex))) {
            final CanonicalTipObservation tipObservation = tipObservations[childIndex];
            if (tipObservation.observedCount == 0) {
                return false;
            }
            final int observedCount = collectObservationPartition(tipObservation);
            if (observedCount == dim) {
                fillContributionForObservedTip(branchAboveParent[childIndex], scratchTransition, tipObservation);
            } else {
                fillContributionForPartiallyObservedTip(
                        branchAboveParent[childIndex], scratchTransition, tipObservation, observedCount);
            }
        } else {
            CanonicalGaussianMessageOps.buildPairPosterior(
                    branchAboveParent[childIndex],
                    scratchTransition,
                    postOrder[childIndex],
                    scratchPairState);
            CanonicalBranchMessageContributionUtils.fillFromPairState(
                    scratchPairState,
                    contributionWorkspace,
                    scratchContribution);
        }

        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                scratchTransition,
                scratchContribution,
                transitionAdjointWorkspace,
                scratchAdjoints);
        return true;
    }

    private void fillContributionForObservedTip(final CanonicalGaussianState aboveState,
                                                final CanonicalGaussianTransition transition,
                                                final CanonicalTipObservation tipObservation) {
        for (int i = 0; i < dim; ++i) {
            double info = transition.informationX[i] + aboveState.information[i];
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                info -= transition.precisionXY[iOff + j] * tipObservation.values[j];
                scratchParentPosterior.precision[iOff + j] =
                        transition.precisionXX[iOff + j] + aboveState.precision[iOff + j];
            }
            scratchParentPosterior.information[i] = info;
        }
        scratchParentPosterior.logNormalizer = 0.0;

        CanonicalGaussianUtils.fillMomentsFromCanonical(
                scratchParentPosterior, scratchMean, scratchCovariance);

        for (int i = 0; i < dim; ++i) {
            final double xi = scratchMean[i];
            final double yi = tipObservation.values[i];
            scratchContribution.dLogL_dInformationX[i] = xi;
            scratchContribution.dLogL_dInformationY[i] = yi;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = scratchMean[j];
                final double yj = tipObservation.values[j];
                final double exx = scratchCovariance[i][j] + xi * xj;
                scratchContribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * exx;
                scratchContribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                scratchContribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                scratchContribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * (yi * yj);
            }
        }
        scratchContribution.dLogL_dLogNormalizer = -1.0;
    }

    private void fillContributionForFixedParentObservedTip(final CanonicalGaussianTransition transition,
                                                           final CanonicalTipObservation tipObservation) {
        clearContribution();
        for (int i = 0; i < dim; ++i) {
            final double xi = fixedRootValue[i];
            final double yi = tipObservation.values[i];
            scratchContribution.dLogL_dInformationX[i] = xi;
            scratchContribution.dLogL_dInformationY[i] = yi;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = fixedRootValue[j];
                final double yj = tipObservation.values[j];
                scratchContribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * (xi * xj);
                scratchContribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                scratchContribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                scratchContribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * (yi * yj);
            }
        }
        scratchContribution.dLogL_dLogNormalizer = -1.0;
    }

    private void fillContributionForFixedParentInternalNode(final CanonicalGaussianTransition transition,
                                                            final CanonicalGaussianState childMessage) {
        CanonicalGaussianMessageOps.conditionOnObservedFirstBlock(transition, fixedRootValue, scratchState);
        CanonicalGaussianMessageOps.combineStates(scratchState, childMessage, scratchParentPosterior);
        CanonicalGaussianUtils.fillMomentsFromCanonical(
                scratchParentPosterior, scratchMean, scratchCovariance);
        fillContributionFromFixedParentChildMoments(scratchMean, scratchCovariance);
    }

    private void fillContributionForFixedParentPartiallyObservedTip(final CanonicalGaussianTransition transition,
                                                                    final CanonicalTipObservation tipObservation,
                                                                    final int observedCount) {
        final int missingCount = dim - observedCount;
        CanonicalGaussianMessageOps.conditionOnObservedFirstBlock(transition, fixedRootValue, scratchState);

        for (int missing = 0; missing < missingCount; ++missing) {
            final int missingTrait = missingIndexScratch[missing];
            double info = scratchState.information[missingTrait];
            final int missingTraitOff = missingTrait * dim;
            for (int observed = 0; observed < observedCount; ++observed) {
                final int observedTrait = observedIndexScratch[observed];
                info -= scratchState.precision[missingTraitOff + observedTrait] * tipObservation.values[observedTrait];
            }
            reducedInformationScratch[missing] = info;
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                reducedPrecisionFlatScratch[missing * missingCount + otherMissing] =
                        scratchState.precision[missingTraitOff + missingIndexScratch[otherMissing]];
            }
        }

        MatrixUtils.safeInvertPrecision(reducedPrecisionFlatScratch, reducedCovarianceFlatScratch, missingCount);

        for (int missing = 0; missing < missingCount; ++missing) {
            double sum = 0.0;
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                sum += reducedCovarianceFlatScratch[missing * missingCount + otherMissing]
                        * reducedInformationScratch[otherMissing];
            }
            scratchMean2[missing] = sum;
        }

        clearContribution();
        for (int i = 0; i < dim; ++i) {
            scratchContribution.dLogL_dInformationX[i] = fixedRootValue[i];
            scratchContribution.dLogL_dInformationY[i] = tipObservation.observed[i]
                    ? tipObservation.values[i]
                    : scratchMean2[reducedIndexByTraitScratch[i] - dim];
        }

        for (int i = 0; i < dim; ++i) {
            final double xi = fixedRootValue[i];
            final double yi = scratchContribution.dLogL_dInformationY[i];
            final int missingI = tipObservation.observed[i] ? -1 : reducedIndexByTraitScratch[i] - dim;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = fixedRootValue[j];
                final double yj = scratchContribution.dLogL_dInformationY[j];
                final int missingJ = tipObservation.observed[j] ? -1 : reducedIndexByTraitScratch[j] - dim;
                final double eyy = (missingI < 0 || missingJ < 0)
                        ? yi * yj
                        : reducedCovarianceFlatScratch[missingI * missingCount + missingJ] + yi * yj;

                scratchContribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * (xi * xj);
                scratchContribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                scratchContribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                scratchContribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * eyy;
            }
        }
        scratchContribution.dLogL_dLogNormalizer = -1.0;
    }

    private void fillContributionFromFixedParentChildMoments(final double[] childMean,
                                                             final double[][] childCovariance) {
        clearContribution();
        for (int i = 0; i < dim; ++i) {
            final double xi = fixedRootValue[i];
            final double yi = childMean[i];
            scratchContribution.dLogL_dInformationX[i] = xi;
            scratchContribution.dLogL_dInformationY[i] = yi;
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                final double xj = fixedRootValue[j];
                final double yj = childMean[j];
                final double eyy = childCovariance[i][j] + yi * yj;
                scratchContribution.dLogL_dPrecisionXX[iOff + j] = -0.5 * (xi * xj);
                scratchContribution.dLogL_dPrecisionXY[iOff + j] = -0.5 * (xi * yj);
                scratchContribution.dLogL_dPrecisionYX[iOff + j] = -0.5 * (yi * xj);
                scratchContribution.dLogL_dPrecisionYY[iOff + j] = -0.5 * eyy;
            }
        }
        scratchContribution.dLogL_dLogNormalizer = -1.0;
    }

    private void fillContributionForPartiallyObservedTip(final CanonicalGaussianState aboveState,
                                                         final CanonicalGaussianTransition transition,
                                                         final CanonicalTipObservation tipObservation,
                                                         final int observedCount) {
        final int missingCount = dim - observedCount;
        final int reducedDimension = dim + missingCount;

        for (int i = 0; i < dim; ++i) {
            double info = transition.informationX[i] + aboveState.information[i];
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                reducedPrecisionFlatScratch[i * reducedDimension + j] =
                        transition.precisionXX[iOff + j] + aboveState.precision[iOff + j];
            }
            for (int observed = 0; observed < observedCount; ++observed) {
                final int observedIndex = observedIndexScratch[observed];
                info -= transition.precisionXY[iOff + observedIndex] * tipObservation.values[observedIndex];
            }
            reducedInformationScratch[i] = info;
            for (int missing = 0; missing < missingCount; ++missing) {
                reducedPrecisionFlatScratch[i * reducedDimension + dim + missing] =
                        transition.precisionXY[iOff + missingIndexScratch[missing]];
            }
        }

        for (int missing = 0; missing < missingCount; ++missing) {
            final int childIndex = missingIndexScratch[missing];
            final int row = dim + missing;
            final int childOff = childIndex * dim;
            double info = transition.informationY[childIndex];
            for (int observed = 0; observed < observedCount; ++observed) {
                final int observedIndex = observedIndexScratch[observed];
                info -= transition.precisionYY[childOff + observedIndex] * tipObservation.values[observedIndex];
            }
            reducedInformationScratch[row] = info;
            for (int j = 0; j < dim; ++j) {
                reducedPrecisionFlatScratch[row * reducedDimension + j] =
                        transition.precisionYX[childOff + j];
            }
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                reducedPrecisionFlatScratch[row * reducedDimension + dim + otherMissing] =
                        transition.precisionYY[childOff + missingIndexScratch[otherMissing]];
            }
        }

        MatrixUtils.safeInvertPrecision(reducedPrecisionFlatScratch, reducedCovarianceFlatScratch, reducedDimension);

        for (int i = 0; i < reducedDimension; ++i) {
            double sum = 0.0;
            for (int j = 0; j < reducedDimension; ++j) {
                sum += reducedCovarianceFlatScratch[i * reducedDimension + j] * reducedInformationScratch[j];
            }
            reducedMeanScratch[i] = Double.isNaN(sum) ? 0.0 : sum;
        }

        clearContribution();
        for (int i = 0; i < dim; ++i) {
            scratchContribution.dLogL_dInformationX[i] = reducedMeanScratch[i];
            scratchContribution.dLogL_dInformationY[i] = tipObservation.observed[i]
                    ? tipObservation.values[i]
                    : reducedMeanScratch[reducedIndexByTraitScratch[i]];
        }

        for (int i = 0; i < dim; ++i) {
            final double xi = reducedMeanScratch[i];
            final int reducedI = reducedIndexByTraitScratch[i];
            for (int j = 0; j < dim; ++j) {
                final double xj = reducedMeanScratch[j];
                final int reducedJ = reducedIndexByTraitScratch[j];
                final double yi = scratchContribution.dLogL_dInformationY[i];
                final double yj = scratchContribution.dLogL_dInformationY[j];

                final double exx = reducedCovarianceFlatScratch[i * reducedDimension + j] + xi * xj;
                final double exy = tipObservation.observed[j]
                        ? xi * yj
                        : reducedCovarianceFlatScratch[i * reducedDimension + reducedJ] + xi * yj;
                final double eyx = tipObservation.observed[i]
                        ? yi * xj
                        : reducedCovarianceFlatScratch[reducedI * reducedDimension + j] + yi * xj;
                final double eyy = (tipObservation.observed[i] || tipObservation.observed[j])
                        ? yi * yj
                        : reducedCovarianceFlatScratch[reducedI * reducedDimension + reducedJ] + yi * yj;

                final int ij = i * dim + j;
                scratchContribution.dLogL_dPrecisionXX[ij] = -0.5 * exx;
                scratchContribution.dLogL_dPrecisionXY[ij] = -0.5 * exy;
                scratchContribution.dLogL_dPrecisionYX[ij] = -0.5 * eyx;
                scratchContribution.dLogL_dPrecisionYY[ij] = -0.5 * eyy;
            }
        }
        scratchContribution.dLogL_dLogNormalizer = -1.0;
    }

    private void clearContribution() {
        Arrays.fill(scratchContribution.dLogL_dInformationX, 0.0);
        Arrays.fill(scratchContribution.dLogL_dInformationY, 0.0);
        Arrays.fill(scratchContribution.dLogL_dPrecisionXX, 0.0);
        Arrays.fill(scratchContribution.dLogL_dPrecisionXY, 0.0);
        Arrays.fill(scratchContribution.dLogL_dPrecisionYX, 0.0);
        Arrays.fill(scratchContribution.dLogL_dPrecisionYY, 0.0);
        scratchContribution.dLogL_dLogNormalizer = -1.0;
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
        final CanonicalTipObservation tipObservation =
                tree.isExternal(tree.getNode(childIndex)) ? tipObservations[childIndex] : null;
        if (tipObservation != null && tipObservation.observedCount < dim) {
            buildTipParentMessageWithMissingMask(
                    tipObservation,
                    provider.getProcessModel(),
                    branchLength,
                    scratchState);
        } else {
            provider.getProcessModel().fillCanonicalTransition(branchLength, scratchTransition);
            buildUpwardParentMessageForTransition(childIndex, scratchTransition, scratchState);
        }
        CanonicalGaussianMessageOps.combineStates(
                branchAboveParent[childIndex], scratchState, scratchCombinedState);
        return CanonicalGaussianMessageOps.normalizationShift(scratchCombinedState, workspace);
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

    private void accumulateRootDiffusionGradient(final double[] gradQ) {
        final int rootIndex = tree.getRoot().getNumber();
        CanonicalGaussianMessageOps.combineStates(
                preOrder[rootIndex], postOrder[rootIndex], scratchCombinedState);

        CanonicalGaussianUtils.fillMomentsFromCanonical(
                scratchCombinedState, scratchMean, scratchCovariance);
        CanonicalGaussianUtils.fillMomentsFromCanonical(
                preOrder[rootIndex], scratchMean2, scratchCovariance2);

        final double[] priorPrecision = preOrder[rootIndex].precision;
        for (int i = 0; i < dim; ++i) {
            final double deltaI = scratchMean[i] - scratchMean2[i];
            for (int j = 0; j < dim; ++j) {
                scratchCovarianceAdjoint[i][j] = scratchCovariance[i][j] + deltaI * (scratchMean[j] - scratchMean2[j]);
            }
        }

        multiplyFlatBySquare(priorPrecision, scratchCovarianceAdjoint, scratchTransitionMatrix, dim);
        multiplySquareByFlat(scratchTransitionMatrix, priorPrecision, scratchCovarianceAdjoint, dim);
        for (int i = 0; i < dim; ++i) {
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                gradQ[iOff + j] += lastRootDiffusionScale
                        * (-0.5 * priorPrecision[iOff + j] + 0.5 * scratchCovarianceAdjoint[i][j]);
            }
        }
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
