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

package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTransitionMomentProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalRootPrior;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservationModel;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.TipObservationMode;

/**
 * Canonical post-order and pre-order tree message propagation.
 */
final class CanonicalTreeTraversal {

    private final Tree tree;
    private final int dimension;
    private final CanonicalTipProjector tipProjector;

    CanonicalTreeTraversal(final Tree tree,
                           final int dimension,
                           final CanonicalTipProjector tipProjector) {
        this.tree = tree;
        this.dimension = dimension;
        this.tipProjector = tipProjector;
    }

    double computePostOrderLogLikelihood(final CanonicalBranchTransitionProvider transitionProvider,
                                         final CanonicalRootPrior rootPrior,
                                         final CanonicalTreeStateStore stateStore,
                                         final BranchGradientWorkspace workspace) {
        if (tree.isExternal(tree.getRoot())) {
            throw new UnsupportedOperationException("Single-node trees are not yet supported by the canonical OU passer.");
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

    void computePreOrder(final CanonicalBranchTransitionProvider transitionProvider,
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

    private void buildUpwardParentMessage(
            final int childIndex,
            final CanonicalBranchTransitionProvider transitionProvider,
            final CanonicalTreeStateStore stateStore,
            final CanonicalGaussianState out,
            final BranchGradientWorkspace workspace) {
        if (tree.isExternal(tree.getNode(childIndex))) {
            final CanonicalTipObservationModel observationModel =
                    stateStore.tipObservationModels[childIndex];
            if (observationModel.getMode() == TipObservationMode.GAUSSIAN_LINK) {
                final CanonicalGaussianTransition transition =
                        transitionFor(childIndex, transitionProvider, workspace);
                observationModel.fillChildCanonicalState(workspace.siblingProduct, workspace.observationWorkspace);
                CanonicalGaussianMessageOps.pushBackward(
                        workspace.siblingProduct, transition, workspace.gaussianWorkspace, out);
                return;
            }
            final CanonicalTipObservation tipObservation = stateStore.tipObservations[childIndex];
            if (observationModel.getMode() != TipObservationMode.EXACT_IDENTITY) {
                buildTipParentMessage(childIndex, tipObservation, transitionProvider, out);
                return;
            }
        }
        final CanonicalGaussianTransition transition =
                transitionFor(childIndex, transitionProvider, workspace);
        buildUpwardParentMessageForTransition(childIndex, transition, stateStore, out, workspace);
    }

    private void buildUpwardParentMessageForTransition(
            final int childIndex,
            final CanonicalGaussianTransition transition,
            final CanonicalTreeStateStore stateStore,
            final CanonicalGaussianState out,
            final BranchGradientWorkspace workspace) {
        if (tree.isExternal(tree.getNode(childIndex))) {
            final CanonicalTipObservationModel observationModel =
                    stateStore.tipObservationModels[childIndex];
            if (observationModel.getMode() == TipObservationMode.GAUSSIAN_LINK) {
                observationModel.fillChildCanonicalState(workspace.siblingProduct, workspace.observationWorkspace);
                CanonicalGaussianMessageOps.pushBackward(
                        workspace.siblingProduct, transition, workspace.gaussianWorkspace, out);
                return;
            }
            final CanonicalTipObservation tipObservation = stateStore.tipObservations[childIndex];
            if (tipObservation.observedCount == 0) {
                CanonicalGaussianMessageOps.clearState(out);
            } else if (tipObservation.observedCount == dimension) {
                CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transition, tipObservation.values, out);
            } else {
                throw new IllegalStateException(
                        "Partially observed canonical tips must use the missing-mask branch update.");
            }
        } else {
            CanonicalGaussianMessageOps.pushBackward(
                    stateStore.postOrder[childIndex],
                    transition,
                    workspace.gaussianWorkspace,
                    out);
        }
    }

    private void buildTipParentMessage(final int childIndex,
                                       final CanonicalTipObservation tipObservation,
                                       final CanonicalBranchTransitionProvider transitionProvider,
                                       final CanonicalGaussianState out) {
        if (tipObservation.observedCount == 0) {
            CanonicalGaussianMessageOps.clearState(out);
            return;
        }

        if (!(transitionProvider instanceof CanonicalTransitionMomentProvider)) {
            throw new UnsupportedOperationException(
                    "Not yet implemented: canonical OU missing-tip propagation currently supports only "
                            + "CanonicalTransitionMomentProvider.");
        }

        transitionProvider.getCanonicalTransitionView(childIndex);
        tipProjector.projectObservedChildToParent(
                tipObservation,
                (CanonicalTransitionMomentProvider) transitionProvider,
                transitionProvider.getEffectiveBranchLength(childIndex),
                out);
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
