/*
 * CanonicalBranchAdjointPreparer.java
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
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalOUTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalPreparedBranchSnapshot;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalPreparedBranchSnapshotProvider;
import dr.util.TaskPool;

/**
 * Prepares frozen branch-local canonical adjoints for tree OU gradients.
 */
final class CanonicalBranchAdjointPreparer {

    interface Source {
        boolean fillLocalAdjointsForBranch(int childIndex,
                                           CanonicalBranchTransitionProvider transitionProvider,
                                           BranchGradientWorkspace workspace);
    }

    private final Tree tree;
    private final int dimension;
    private final int nodeCount;
    private final TaskPool taskPool;
    private final BranchGradientWorkspace mainWorkspace;
    private final BranchGradientWorkspace[] workspaces;
    private final Source source;

    CanonicalBranchAdjointPreparer(final Tree tree,
                                   final int dimension,
                                   final TaskPool taskPool,
                                   final BranchGradientWorkspace mainWorkspace,
                                   final BranchGradientWorkspace[] workspaces,
                                   final Source source) {
        this.tree = tree;
        this.dimension = dimension;
        this.nodeCount = tree.getNodeCount();
        this.taskPool = taskPool;
        this.mainWorkspace = mainWorkspace;
        this.workspaces = workspaces;
        this.source = source;
    }

    void prepare(final CanonicalBranchTransitionProvider transitionProvider,
                 final CanonicalTreeStateStore stateStore,
                 final BranchGradientInputs out) {
        out.checkCompatible(Math.max(0, nodeCount - 1), dimension);
        out.clear();

        final CanonicalOUTransitionProvider ouProvider = requireOUProvider(transitionProvider);
        final OUProcessModel processModel = ouProvider.getProcessModel();
        final OrthogonalBlockCanonicalParameterization orthogonalSelection =
                processModel.getSelectionMatrixParameterization()
                        instanceof OrthogonalBlockCanonicalParameterization
                        ? (OrthogonalBlockCanonicalParameterization)
                        processModel.getSelectionMatrixParameterization()
                        : null;
        final CanonicalPreparedBranchSnapshotProvider snapshotProvider =
                requirePreparedBranchSnapshotProvider(transitionProvider);
        final int rootIndex = tree.getRoot().getNumber();
        if (workspaces.length <= 1 || nodeCount <= 2) {
            for (int childIndex = 0; childIndex < nodeCount; ++childIndex) {
                if (childIndex == rootIndex) {
                    continue;
                }
                if (!source.fillLocalAdjointsForBranch(childIndex, transitionProvider, mainWorkspace)) {
                    continue;
                }
                out.addBranch(requireSnapshot(snapshotProvider, childIndex), mainWorkspace.adjoints);
            }
        } else {
            taskPool.forkDynamicBalanced(
                    nodeCount,
                    branchGradientPreparationChunkSize(nodeCount),
                    (childIndex, thread) -> {
                        if (childIndex == rootIndex) {
                            out.clearStagedBranch(childIndex);
                            return;
                        }

                        final BranchGradientWorkspace workspace = workspaces[thread];
                        if (!source.fillLocalAdjointsForBranch(childIndex, transitionProvider, workspace)) {
                            out.clearStagedBranch(childIndex);
                            return;
                        }

                        out.stageBranch(
                                requireSnapshot(snapshotProvider, childIndex),
                                workspace.adjoints);
                    });
            out.compactStagedBranches(rootIndex, orthogonalSelection != null);
        }

        out.setRoot(
                stateStore.lastRootDiffusionScale,
                stateStore.preOrder[rootIndex],
                stateStore.postOrder[rootIndex]);
    }

    private int branchGradientPreparationChunkSize(final int taskLimit) {
        return branchGradientChunkSize(taskLimit, 2, 16);
    }

    private int branchGradientChunkSize(final int taskLimit,
                                        final int targetChunksPerWorker,
                                        final int maxChunkSize) {
        final int workerCount = Math.max(1, workspaces.length);
        final int suggested =
                (taskLimit + workerCount * targetChunksPerWorker - 1) / (workerCount * targetChunksPerWorker);
        return Math.max(1, Math.min(maxChunkSize, suggested));
    }

    private static CanonicalOUTransitionProvider requireOUProvider(
            final CanonicalBranchTransitionProvider transitionProvider) {
        if (!(transitionProvider instanceof CanonicalOUTransitionProvider)) {
            throw new UnsupportedOperationException(
                    "Not yet implemented: canonical OU gradients currently support only "
                            + "CanonicalOUTransitionProvider implementations.");
        }
        return (CanonicalOUTransitionProvider) transitionProvider;
    }

    private static CanonicalPreparedBranchSnapshotProvider requirePreparedBranchSnapshotProvider(
            final CanonicalBranchTransitionProvider transitionProvider) {
        if (!(transitionProvider instanceof CanonicalPreparedBranchSnapshotProvider)) {
            throw new UnsupportedOperationException(
                    "Canonical OU gradients require CanonicalPreparedBranchSnapshotProvider.");
        }
        return (CanonicalPreparedBranchSnapshotProvider) transitionProvider;
    }

    private static CanonicalPreparedBranchSnapshot requireSnapshot(
            final CanonicalPreparedBranchSnapshotProvider snapshotProvider,
            final int childIndex) {
        final CanonicalPreparedBranchSnapshot snapshot = snapshotProvider.getPreparedBranchSnapshot(childIndex);
        if (snapshot == null) {
            throw new IllegalStateException("Missing prepared canonical branch snapshot for child " + childIndex);
        }
        return snapshot;
    }
}
