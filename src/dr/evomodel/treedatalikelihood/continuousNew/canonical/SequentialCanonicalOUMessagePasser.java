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

package dr.evomodel.treedatalikelihood.continuousNew.canonical;

import dr.evomodel.continuous.canonical.contribution.CanonicalBranchContributionAssembler;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientInputs;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientMessageProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.CanonicalBranchAdjointPreparer;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.CanonicalBranchLengthGradientEngine;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.CanonicalTreeGradientEngine;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.StateStoreBranchGradientMessageProvider;
import dr.evomodel.continuous.canonical.scheduling.DimensionWeightedChunkSizeStrategy;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.traversal.CanonicalTreeStateStore;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.traversal.CanonicalTreeTraversal;
import dr.evomodel.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.continuous.canonical.workspace.WorkspaceFactory;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuousNew.CanonicalDebugOptions;
import dr.evomodel.treedatalikelihood.continuousNew.CanonicalGradientFallbackPolicy;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.CanonicalTipObservationModel;
import dr.util.TaskPool;

/**
 * <p>This class implements the canonical OU tree traversals and exact/partially
 * observed tip elimination in canonical form.
 */
public final class SequentialCanonicalOUMessagePasser implements CanonicalTreeMessagePasser {

    private final int dim;
    private final int nodeCount;
    private final int tipCount;
    private final CanonicalTreeStateStore stateStore;
    private final CanonicalTreeTraversal treeTraversal;
    private final BranchGradientWorkspace messageWorkspace;
    private final CanonicalBranchContributionAssembler branchContributionAssembler;
    private final BranchGradientInputs preparedBranchGradientInputs;
    private final CanonicalTraversalRunner traversalRunner;
    private final CanonicalGradientPreparationRunner gradientPreparationRunner;
    private final CanonicalTreeGradientEngine treeGradientEngine;
    private final CanonicalBranchLengthGradientRunner branchLengthGradientRunner;

    public SequentialCanonicalOUMessagePasser(final Tree tree, final int dim) {
        this(tree, dim, CanonicalGradientFallbackPolicy.branchGradientParallelismFromSystemProperties());
    }

    public static SequentialCanonicalOUMessagePasser likelihoodOnly(final Tree tree, final int dim) {
        return new SequentialCanonicalOUMessagePasser(
                tree,
                dim,
                1,
                CanonicalDebugOptions.fromSystemProperties(),
                CanonicalGradientFallbackPolicy.fromSystemProperties(),
                false);
    }

    public SequentialCanonicalOUMessagePasser(final Tree tree,
                                              final int dim,
                                              final int branchGradientParallelism) {
        this(tree,
                dim,
                branchGradientParallelism,
                CanonicalDebugOptions.fromSystemProperties(),
                CanonicalGradientFallbackPolicy.fromSystemProperties(),
                true);
    }

    SequentialCanonicalOUMessagePasser(final Tree tree,
                                       final int dim,
                                       final int branchGradientParallelism,
                                       final CanonicalDebugOptions debugOptions,
                                       final CanonicalGradientFallbackPolicy fallbackPolicy) {
        this(tree, dim, branchGradientParallelism, debugOptions, fallbackPolicy, true);
    }

    private SequentialCanonicalOUMessagePasser(final Tree tree,
                                               final int dim,
                                               final int branchGradientParallelism,
                                               final CanonicalDebugOptions debugOptions,
                                               final CanonicalGradientFallbackPolicy fallbackPolicy,
                                               final boolean gradientSupport) {
        if (debugOptions == null) {
            throw new IllegalArgumentException("debugOptions must not be null");
        }
        if (fallbackPolicy == null) {
            throw new IllegalArgumentException("fallbackPolicy must not be null");
        }
        this.dim = dim;
        this.nodeCount = tree.getNodeCount();
        this.tipCount = tree.getExternalNodeCount();
        this.stateStore = new CanonicalTreeStateStore(nodeCount, dim);
        this.branchContributionAssembler = new CanonicalBranchContributionAssembler(dim);
        final BranchGradientWorkspace mainWorkspace = WorkspaceFactory.branchGradientWorkspace(dim);
        this.messageWorkspace = mainWorkspace;
        final BranchGradientMessageProvider branchGradientMessageProvider =
                new StateStoreBranchGradientMessageProvider(tree, stateStore);
        final TaskPool transitionPreloadTaskPool = new TaskPool(
                nodeCount,
                Math.max(1, Math.min(branchGradientParallelism, fallbackPolicy.getTransitionPreloadParallelism())));
        final CanonicalTreeTraversal treeTraversal = new CanonicalTreeTraversal(tree, dim);
        this.treeTraversal = treeTraversal;
        this.traversalRunner = new CanonicalTraversalRunner(
                treeTraversal,
                stateStore,
                mainWorkspace,
                tree.getRoot().getNumber(),
                transitionPreloadTaskPool,
                new DimensionWeightedChunkSizeStrategy(dim, Math.max(1, transitionPreloadTaskPool.getNumThreads()))
                        .chunkSize(nodeCount));
        if (gradientSupport) {
            this.preparedBranchGradientInputs = new BranchGradientInputs(Math.max(0, nodeCount - 1), dim);
            final TaskPool branchGradientTaskPool = new TaskPool(nodeCount, Math.max(1, branchGradientParallelism));
            final int branchGradientWorkspaceCount =
                    branchGradientTaskPool.getNumThreads() <= 1 ? 1 : branchGradientTaskPool.getNumThreads() + 1;
            final BranchGradientWorkspace[] branchGradientWorkspaces =
                    new BranchGradientWorkspace[branchGradientWorkspaceCount];
            for (int i = 0; i < branchGradientWorkspaces.length; ++i) {
                branchGradientWorkspaces[i] = WorkspaceFactory.branchGradientWorkspace(dim);
            }
            final CanonicalBranchAdjointPreparer branchAdjointPreparer = new CanonicalBranchAdjointPreparer(
                    tree,
                    dim,
                    branchGradientTaskPool,
                    mainWorkspace,
                    branchGradientWorkspaces,
                    branchGradientMessageProvider,
                    branchContributionAssembler::fillLocalAdjointsForBranch);
            this.gradientPreparationRunner = new CanonicalGradientPreparationRunner(
                    debugOptions,
                    traversalRunner,
                    branchAdjointPreparer);
            this.treeGradientEngine = new CanonicalTreeGradientEngine(
                    dim,
                    Math.max(0, nodeCount - 1),
                    branchGradientTaskPool,
                    mainWorkspace,
                    branchGradientWorkspaces);
            this.branchLengthGradientRunner = new CanonicalBranchLengthGradientRunner(
                    tree,
                    fallbackPolicy,
                    stateStore,
                    mainWorkspace,
                    traversalRunner,
                    gradientPreparationRunner,
                    preparedBranchGradientInputs,
                    new CanonicalBranchLengthGradientEngine());
        } else {
            this.preparedBranchGradientInputs = null;
            this.gradientPreparationRunner = null;
            this.treeGradientEngine = null;
            this.branchLengthGradientRunner = null;
        }
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
        stateStore.setTipObservation(tipIndex, observation);
    }

    public void setTipObservationModel(final int tipIndex,
                                       final CanonicalTipObservationModel observationModel) {
        stateStore.setTipObservationModel(tipIndex, observationModel);
    }

    @Override
    public double computePostOrderLogLikelihood(final CanonicalBranchTransitionProvider transitionProvider,
                                                final CanonicalRootPrior rootPrior) {
        return traversalRunner.computePostOrderLogLikelihood(transitionProvider, rootPrior);
    }

    @Override
    public CanonicalGaussianState getPostOrderState(final int nodeIndex) {
        return stateStore.postOrder[nodeIndex];
    }

    @Override
    public void fillPostOrderNodeState(final int nodeIndex, final CanonicalGaussianState out) {
        requirePostOrderState();
        checkOutputState(out);
        treeTraversal.fillPostOrderNodeState(nodeIndex, stateStore, messageWorkspace, out);
    }

    @Override
    public void fillPostOrderBranchBottomState(final int childNodeIndex,
                                               final CanonicalGaussianState out) {
        requirePostOrderState();
        checkOutputState(out);
        treeTraversal.fillPostOrderBranchBottomState(childNodeIndex, stateStore, messageWorkspace, out);
    }

    @Override
    public void fillPostOrderBranchTopState(final int childNodeIndex,
                                            final CanonicalBranchTransitionProvider transitionProvider,
                                            final CanonicalGaussianState out) {
        requirePostOrderState();
        checkOutputState(out);
        treeTraversal.fillPostOrderBranchTopState(
                childNodeIndex,
                transitionProvider,
                stateStore,
                messageWorkspace,
                out);
    }

    @Override
    public void computePreOrder(final CanonicalBranchTransitionProvider transitionProvider,
                                final CanonicalRootPrior rootPrior) {
        traversalRunner.computePreOrder(transitionProvider, rootPrior);
    }

    @Override
    public void fillPreOrderBranchTopState(final int childNodeIndex,
                                           final CanonicalGaussianState out) {
        requirePreOrderState();
        checkOutputState(out);
        treeTraversal.fillPreOrderBranchTopState(childNodeIndex, stateStore, out);
    }

    @Override
    public void fillPreOrderBranchBottomState(final int childNodeIndex,
                                              final CanonicalBranchTransitionProvider transitionProvider,
                                              final CanonicalGaussianState out) {
        requirePreOrderState();
        checkOutputState(out);
        treeTraversal.fillPreOrderBranchBottomState(
                childNodeIndex,
                transitionProvider,
                stateStore,
                messageWorkspace,
                out);
    }

    @Override
    public void fillRootPreOrderState(final CanonicalGaussianState out) {
        requirePreOrderState();
        checkOutputState(out);
        treeTraversal.fillRootPreOrderState(stateStore, out);
    }

    @Override
    public void computeGradientBranchLengths(final CanonicalBranchTransitionProvider transitionProvider,
                                             final double[] gradT) {
        requireGradientSupport();
        branchLengthGradientRunner.compute(transitionProvider, gradT);
    }

    @Override
    public void computeGradientBranchLengths(final CanonicalBranchTransitionProvider transitionProvider,
                                             final BranchGradientInputs inputs,
                                             final double[] gradT) {
        requireGradientSupport();
        branchLengthGradientRunner.compute(transitionProvider, inputs, gradT);
    }

    @Override
    public void computeJointGradients(final CanonicalBranchTransitionProvider transitionProvider,
                                      final double[] gradA,
                                      final double[] gradQ,
                                      final double[] gradMu) {
        requireGradientSupport();
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
        prepareBranchGradientInputs(transitionProvider, CanonicalTransitionCachePhases.GRADIENT_PREP, out);
    }

    @Override
    public void prepareBranchGradientInputs(final CanonicalBranchTransitionProvider transitionProvider,
                                            final String phase,
                                            final BranchGradientInputs out) {
        requireGradientSupport();
        gradientPreparationRunner.prepare(transitionProvider, phase, out);
    }

    public void computeJointGradients(final CanonicalBranchTransitionProvider transitionProvider,
                                      final BranchGradientInputs inputs,
                                      final double[] gradA,
                                      final double[] gradQ,
                                      final double[] gradMu) {
        requireGradientSupport();
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

    private void requirePostOrderState() {
        if (!stateStore.hasPostOrderState) {
            throw new IllegalStateException(
                    "Canonical post-order messages have not been computed.");
        }
    }

    private void requirePreOrderState() {
        if (!stateStore.hasPreOrderState) {
            throw new IllegalStateException(
                    "Canonical pre-order messages have not been computed.");
        }
    }

    private void checkOutputState(final CanonicalGaussianState out) {
        if (out == null) {
            throw new IllegalArgumentException("Output state must not be null.");
        }
        if (out.getDimension() != dim) {
            throw new IllegalArgumentException(
                    "Output state dimension mismatch: " + out.getDimension() + " vs " + dim + ".");
        }
    }

    private void requireGradientSupport() {
        if (gradientPreparationRunner == null || treeGradientEngine == null || branchLengthGradientRunner == null) {
            throw new UnsupportedOperationException("This canonical message passer was created for likelihood-only use.");
        }
    }

}
