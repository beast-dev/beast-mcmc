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

package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

import dr.evolution.tree.Tree;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedTransitionCapability;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.SelectionMatrixParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalOUProviderSupport;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalOUTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalPreparedBranchSnapshot;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalPreparedBranchSnapshotProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.scheduling.ChunkSizeStrategy;
import dr.evomodel.treedatalikelihood.continuous.canonical.scheduling.DimensionWeightedChunkSizeStrategy;
import dr.evomodel.treedatalikelihood.continuous.canonical.traversal.CanonicalTreeStateStore;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.util.TaskPool;

/**
 * Prepares frozen branch-local canonical adjoints for tree OU gradients.
 */
public final class CanonicalBranchAdjointPreparer {

    public interface Source {
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
    private final ChunkSizeStrategy chunkSizeStrategy;

    public CanonicalBranchAdjointPreparer(final Tree tree,
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
        this.chunkSizeStrategy =
                new DimensionWeightedChunkSizeStrategy(dimension, Math.max(1, workspaces.length));
    }

    public void prepare(final CanonicalBranchTransitionProvider transitionProvider,
                 final CanonicalTreeStateStore stateStore,
                 final BranchGradientInputs out) {
        out.checkCompatible(Math.max(0, nodeCount - 1), dimension);
        out.clear();

        final CanonicalOUTransitionProvider ouProvider =
                CanonicalOUProviderSupport.requireCapability(transitionProvider, CanonicalOUTransitionProvider.class);
        final OUProcessModel processModel = ouProvider.getProcessModel();
        final SelectionMatrixParameterization parameterization =
                processModel.getSelectionMatrixParameterization();
        final CanonicalPreparedTransitionCapability preparedTransition =
                parameterization instanceof CanonicalPreparedTransitionCapability
                        ? (CanonicalPreparedTransitionCapability) parameterization
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
            out.compactStagedBranches(rootIndex, preparedTransition != null);
        }
        out.sortActiveBranchesByLength();

        out.setRoot(
                stateStore.lastRootDiffusionScale,
                stateStore.preOrder[rootIndex],
                stateStore.postOrder[rootIndex]);
    }

    private int branchGradientPreparationChunkSize(final int taskLimit) {
        return chunkSizeStrategy.chunkSize(taskLimit);
    }

    private static CanonicalPreparedBranchSnapshotProvider requirePreparedBranchSnapshotProvider(
            final CanonicalBranchTransitionProvider transitionProvider) {
        return CanonicalOUProviderSupport.requireCapability(
                transitionProvider, CanonicalPreparedBranchSnapshotProvider.class);
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
