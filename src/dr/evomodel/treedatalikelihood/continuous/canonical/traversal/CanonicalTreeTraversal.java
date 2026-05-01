/*
 * CanonicalTreeTraversal.java
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

package dr.evomodel.treedatalikelihood.continuous.canonical.traversal;

import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.*;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalRootPrior;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTransitionMomentProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservationModel;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.TipObservationMode;

/**
 * Canonical post-order and pre-order tree message propagation.
 *
 * <p>Tip observation models are dispatched polymorphically through
 * {@link CanonicalTipObservationModel#fillParentMessage}. This class contains no
 * per-mode logic for how individual tip types produce their parent messages.
 */
public final class CanonicalTreeTraversal {

    private final Tree tree;
    private final int dimension;

    public CanonicalTreeTraversal(final Tree tree, final int dimension) {
        this.tree      = tree;
        this.dimension = dimension;
    }

    public double computePostOrderLogLikelihood(final CanonicalBranchTransitionProvider transitionProvider,
                                         final CanonicalRootPrior rootPrior,
                                         final CanonicalTreeStateStore stateStore,
                                         final BranchGradientWorkspace workspace) {
        if (tree.isExternal(tree.getRoot())) {
            throw new UnsupportedOperationException(
                    "Single-node trees are not yet supported by the canonical OU passer.");
        }

        final int nodeCount = tree.getNodeCount();
        for (int i = 0; i < nodeCount; i++) {
            if (!tree.isExternal(tree.getNode(i))) {
                computePostOrderAtInternalNode(i, transitionProvider, stateStore, workspace);
            }
        }

        final int rootIndex = tree.getRoot().getNumber();
        transitionProvider.fillTraitCovariance(workspace.traitCovariance);
        stateStore.hasPostOrderState = true;
        return rootPrior.computeLogMarginalLikelihood(
                stateStore.postOrder[rootIndex],
                workspace.traitCovariance);
    }

    public void computePreOrder(final CanonicalBranchTransitionProvider transitionProvider,
                         final CanonicalRootPrior rootPrior,
                         final CanonicalTreeStateStore stateStore,
                         final BranchGradientWorkspace workspace) {
        final int rootIndex = tree.getRoot().getNumber();
        transitionProvider.fillTraitCovariance(workspace.traitCovariance);
        CanonicalGaussianMessageOps.clearState(stateStore.preOrder[rootIndex]);
        stateStore.hasFixedRootValue = rootPrior.isFixedRoot();
        if (stateStore.hasFixedRootValue) {
            rootPrior.fillFixedRootValue(stateStore.fixedRootValue);
            stateStore.lastRootDiffusionScale = 0.0;
        } else {
            rootPrior.fillRootPriorState(workspace.traitCovariance, stateStore.preOrder[rootIndex]);
            stateStore.lastRootDiffusionScale = rootPrior.getDiffusionScale();
        }
        CanonicalGaussianMessageOps.clearState(stateStore.branchAboveParent[rootIndex]);
        computePreOrderRecursive(rootIndex, transitionProvider, stateStore, workspace);
        stateStore.hasPreOrderState = true;
    }

    private void computePostOrderAtInternalNode(
            final int nodeIndex,
            final CanonicalBranchTransitionProvider transitionProvider,
            final CanonicalTreeStateStore stateStore,
            final BranchGradientWorkspace workspace) {
        final CanonicalGaussianState dest = stateStore.postOrder[nodeIndex];
        CanonicalGaussianMessageOps.clearState(dest);

        final int childCount = tree.getChildCount(tree.getNode(nodeIndex));
        boolean first = true;
        for (int c = 0; c < childCount; c++) {
            final int childIndex = tree.getChild(tree.getNode(nodeIndex), c).getNumber();
            buildUpwardParentMessage(childIndex, transitionProvider, stateStore, workspace.state, workspace);
            if (first) {
                CanonicalGaussianMessageOps.copyState(workspace.state, dest);
                first = false;
            } else {
                CanonicalGaussianMessageOps.combineStateInPlace(dest, workspace.state);
            }
        }
    }

    private void computePreOrderRecursive(
            final int parentIndex,
            final CanonicalBranchTransitionProvider transitionProvider,
            final CanonicalTreeStateStore stateStore,
            final BranchGradientWorkspace workspace) {
        if (tree.isExternal(tree.getNode(parentIndex))) {
            return;
        }

        final int rootIndex = tree.getRoot().getNumber();
        final int childCount = tree.getChildCount(tree.getNode(parentIndex));
        for (int c = 0; c < childCount; c++) {
            final int childIndex = tree.getChild(tree.getNode(parentIndex), c).getNumber();
            if (stateStore.hasFixedRootValue && parentIndex == rootIndex) {
                CanonicalGaussianMessageOps.clearState(stateStore.branchAboveParent[childIndex]);
                final CanonicalGaussianTransition transition =
                        transitionFor(childIndex, transitionProvider, workspace);
                CanonicalGaussianMessageOps.conditionOnObservedFirstBlock(
                        transition,
                        stateStore.fixedRootValue,
                        stateStore.preOrder[childIndex]);
                computePreOrderRecursive(childIndex, transitionProvider, stateStore, workspace);
                continue;
            }

            final boolean hasSiblings = buildSiblingProduct(
                    parentIndex, childIndex, transitionProvider, stateStore, workspace.siblingProduct, workspace);

            if (hasSiblings) {
                CanonicalGaussianMessageOps.combineStates(
                        stateStore.preOrder[parentIndex],
                        workspace.siblingProduct,
                        workspace.downwardParentState);
            } else {
                CanonicalGaussianMessageOps.copyState(
                        stateStore.preOrder[parentIndex],
                        workspace.downwardParentState);
            }

            CanonicalGaussianMessageOps.copyState(
                    workspace.downwardParentState,
                    stateStore.branchAboveParent[childIndex]);
            final CanonicalGaussianTransition transition =
                    transitionFor(childIndex, transitionProvider, workspace);
            CanonicalGaussianMessageOps.pushForward(
                    workspace.downwardParentState,
                    transition,
                    workspace.gaussianWorkspace,
                    stateStore.preOrder[childIndex]);
            computePreOrderRecursive(childIndex, transitionProvider, stateStore, workspace);
        }
    }

    private boolean buildSiblingProduct(
            final int parentIndex,
            final int excludedChildIndex,
            final CanonicalBranchTransitionProvider transitionProvider,
            final CanonicalTreeStateStore stateStore,
            final CanonicalGaussianState out,
            final BranchGradientWorkspace workspace) {
        CanonicalGaussianMessageOps.clearState(out);
        final int childCount = tree.getChildCount(tree.getNode(parentIndex));
        boolean found = false;
        for (int c = 0; c < childCount; c++) {
            final int siblingIndex = tree.getChild(tree.getNode(parentIndex), c).getNumber();
            if (siblingIndex == excludedChildIndex) {
                continue;
            }
            buildUpwardParentMessage(siblingIndex, transitionProvider, stateStore, workspace.state, workspace);
            if (!found) {
                CanonicalGaussianMessageOps.copyState(workspace.state, out);
                found = true;
            } else {
                CanonicalGaussianMessageOps.combineStateInPlace(out, workspace.state);
            }
        }
        return found;
    }

    /**
     * Builds the upward parent message from child {@code childIndex} into {@code out}.
     *
     * <p>For tip nodes, dispatches to the tip's observation model via
     * {@link CanonicalTipObservationModel#fillParentMessage}, which owns all mode-specific
     * logic (exact identity, partial identity, missing, Gaussian link). For internal
     * nodes, applies the standard Gaussian backward push.
     */
    private void buildUpwardParentMessage(
            final int childIndex,
            final CanonicalBranchTransitionProvider transitionProvider,
            final CanonicalTreeStateStore stateStore,
            final CanonicalGaussianState out,
            final BranchGradientWorkspace workspace) {
        final CanonicalGaussianTransition transition = transitionFor(childIndex, transitionProvider, workspace);
        if (tree.isExternal(tree.getNode(childIndex))) {
            final CanonicalTipObservationModel observationModel =
                    stateStore.tipObservationModels[childIndex];
            final CanonicalTransitionMomentProvider momentProvider =
                    momentProviderIfNeeded(observationModel, transitionProvider);
            final double branchLength = transitionProvider.getEffectiveBranchLength(childIndex);
            observationModel.fillParentMessage(
                    transition,
                    momentProvider,
                    branchLength,
                    workspace.tipParentMessageWorkspace,
                    workspace.gaussianWorkspace,
                    out);
        } else {
            CanonicalGaussianMessageOps.pushBackward(
                    stateStore.postOrder[childIndex], transition, workspace.gaussianWorkspace, out);
        }
    }

    /**
     * Returns the {@link CanonicalTransitionMomentProvider} required for partial-identity tip
     * observations, or {@code null} for all other modes.
     *
     * @throws UnsupportedOperationException if the mode is {@code PARTIAL_EXACT_IDENTITY} but the
     *                                        transition provider does not implement
     *                                        {@link CanonicalTransitionMomentProvider}
     */
    private CanonicalTransitionMomentProvider momentProviderIfNeeded(
            final CanonicalTipObservationModel observationModel,
            final CanonicalBranchTransitionProvider transitionProvider) {
        if (observationModel.getMode() != TipObservationMode.PARTIAL_EXACT_IDENTITY) {
            return null;
        }
        if (transitionProvider instanceof CanonicalTransitionMomentProvider) {
            return (CanonicalTransitionMomentProvider) transitionProvider;
        }
        throw new UnsupportedOperationException(
                "Canonical OU partial exact identity observations require a "
                + "CanonicalTransitionMomentProvider, but the current provider ("
                + transitionProvider.getClass().getSimpleName() + ") does not implement it.");
    }

    private CanonicalGaussianTransition transitionFor(final int childIndex,
                                                      final CanonicalBranchTransitionProvider transitionProvider,
                                                      final BranchGradientWorkspace workspace) {
        final CanonicalGaussianTransition transitionView =
                transitionProvider.getCanonicalTransitionView(childIndex);
        if (transitionView != null) {
            return transitionView;
        }
        transitionProvider.fillCanonicalTransition(childIndex, workspace.transition);
        return workspace.transition;
    }
}
