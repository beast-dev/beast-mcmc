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

package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.contribution.CanonicalBranchContributionAssembler;
import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.BranchGradientInputs;
import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.CanonicalBranchAdjointPreparer;
import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.CanonicalBranchLengthGradientEngine;
import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.CanonicalTreeGradientEngine;
import dr.evomodel.treedatalikelihood.continuous.canonical.traversal.CanonicalTreeStateStore;
import dr.evomodel.treedatalikelihood.continuous.canonical.traversal.CanonicalTreeTraversal;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.WorkspaceFactory;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.CanonicalDebugOptions;
import dr.evomodel.treedatalikelihood.continuous.CanonicalGradientFallbackPolicy;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservationModel;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.TipObservationMode;
import dr.util.TaskPool;

import java.util.Arrays;

/**
 * <p>This class implements the canonical OU tree traversals and exact/partially
 * observed tip elimination in canonical form.
 */
public final class SequentialCanonicalOUMessagePasser implements CanonicalTreeMessagePasser {

    private static final double BRANCH_LENGTH_FD_RELATIVE_STEP = 1.0e-6;
    private static final double BRANCH_LENGTH_FD_ABSOLUTE_STEP = 1.0e-8;

    private final Tree tree;
    private final int dim;
    private final int nodeCount;
    private final int tipCount;
    private final CanonicalDebugOptions debugOptions;
    private final CanonicalGradientFallbackPolicy fallbackPolicy;
    private final CanonicalTreeStateStore stateStore;
    private final CanonicalTreeTraversal treeTraversal;
    private final CanonicalBranchContributionAssembler branchContributionAssembler;
    private final BranchGradientInputs preparedBranchGradientInputs;
    private final BranchGradientWorkspace mainWorkspace;
    private final TaskPool branchGradientTaskPool;
    private final BranchGradientWorkspace[] branchGradientWorkspaces;
    private final CanonicalBranchAdjointPreparer branchAdjointPreparer;
    private final CanonicalTreeGradientEngine treeGradientEngine;
    private final CanonicalBranchLengthGradientEngine branchLengthGradientEngine;

    public SequentialCanonicalOUMessagePasser(final Tree tree, final int dim) {
        this(tree, dim, CanonicalGradientFallbackPolicy.branchGradientParallelismFromSystemProperties());
    }

    public SequentialCanonicalOUMessagePasser(final Tree tree,
                                              final int dim,
                                              final int branchGradientParallelism) {
        this(tree,
                dim,
                branchGradientParallelism,
                CanonicalDebugOptions.fromSystemProperties(),
                CanonicalGradientFallbackPolicy.fromSystemProperties());
    }

    SequentialCanonicalOUMessagePasser(final Tree tree,
                                       final int dim,
                                       final int branchGradientParallelism,
                                       final CanonicalDebugOptions debugOptions,
                                       final CanonicalGradientFallbackPolicy fallbackPolicy) {
        if (debugOptions == null) {
            throw new IllegalArgumentException("debugOptions must not be null");
        }
        if (fallbackPolicy == null) {
            throw new IllegalArgumentException("fallbackPolicy must not be null");
        }
        this.tree = tree;
        this.dim = dim;
        this.nodeCount = tree.getNodeCount();
        this.tipCount = tree.getExternalNodeCount();
        this.debugOptions = debugOptions;
        this.fallbackPolicy = fallbackPolicy;
        this.stateStore = new CanonicalTreeStateStore(nodeCount, dim);
        this.treeTraversal = new CanonicalTreeTraversal(tree, dim);
        this.branchContributionAssembler = new CanonicalBranchContributionAssembler(tree, dim, stateStore);
        this.preparedBranchGradientInputs = new BranchGradientInputs(Math.max(0, nodeCount - 1), dim);
        this.mainWorkspace = WorkspaceFactory.branchGradientWorkspace(dim);
        this.branchGradientTaskPool = new TaskPool(nodeCount, Math.max(1, branchGradientParallelism));
        final int branchGradientWorkspaceCount =
                branchGradientTaskPool.getNumThreads() <= 1 ? 1 : branchGradientTaskPool.getNumThreads() + 1;
        this.branchGradientWorkspaces = new BranchGradientWorkspace[branchGradientWorkspaceCount];
        for (int i = 0; i < branchGradientWorkspaces.length; ++i) {
            branchGradientWorkspaces[i] = WorkspaceFactory.branchGradientWorkspace(dim);
        }
        this.branchAdjointPreparer = new CanonicalBranchAdjointPreparer(
                tree,
                dim,
                branchGradientTaskPool,
                mainWorkspace,
                branchGradientWorkspaces,
                branchContributionAssembler::fillLocalAdjointsForBranch);
        this.treeGradientEngine = new CanonicalTreeGradientEngine(
                dim,
                Math.max(0, nodeCount - 1),
                branchGradientTaskPool,
                mainWorkspace,
                branchGradientWorkspaces);
        this.branchLengthGradientEngine = new CanonicalBranchLengthGradientEngine();
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
        final String previousPhase = pushTransitionCachePhase(
                transitionProvider, CanonicalTransitionCachePhases.POSTORDER);
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
        final String previousPhase = pushTransitionCachePhase(
                transitionProvider, CanonicalTransitionCachePhases.PREORDER);
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
    public void computeGradientBranchLengths(final CanonicalBranchTransitionProvider transitionProvider,
                                             final double[] gradT) {
        ensureGradientState();

        final String previousPhase = pushTransitionCachePhase(
                transitionProvider, CanonicalTransitionCachePhases.BRANCH_LENGTH_GRADIENT);
        try {
            final CanonicalOUTransitionProvider ouProvider =
                    CanonicalOUProviderSupport.requireOUProvider(transitionProvider);
            if (!fallbackPolicy.useBranchLengthFiniteDifference()) {
                final long missesBefore =
                        transitionCacheMisses(
                                transitionProvider,
                                CanonicalTransitionCachePhases.BRANCH_LENGTH_GRADIENT);
                branchAdjointPreparer.prepare(transitionProvider, stateStore, preparedBranchGradientInputs);
                assertNoGradientTransitionMisses(
                        transitionProvider,
                        CanonicalTransitionCachePhases.BRANCH_LENGTH_GRADIENT,
                        missesBefore);
                branchLengthGradientEngine.compute(
                        ouProvider.getProcessModel(),
                        preparedBranchGradientInputs,
                        gradT);
                return;
            }

            Arrays.fill(gradT, 0.0);
            final int rootIndex = tree.getRoot().getNumber();
            for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
                if (childIndex == rootIndex) {
                    continue;
                }
                final CanonicalTipObservationModel observationModel =
                        tree.isExternal(tree.getNode(childIndex)) ? stateStore.tipObservationModels[childIndex] : null;
                if (observationModel != null && observationModel.isEmpty()) {
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
        final String previousPhase = pushTransitionCachePhase(
                transitionProvider, CanonicalTransitionCachePhases.GRADIENT_PREP);
        try {
            ensureGradientState();
            final long missesBefore = transitionCacheMisses(
                    transitionProvider, CanonicalTransitionCachePhases.GRADIENT_PREP);
            branchAdjointPreparer.prepare(transitionProvider, stateStore, out);
            assertNoGradientTransitionMisses(
                    transitionProvider, CanonicalTransitionCachePhases.GRADIENT_PREP, missesBefore);
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

    private void ensureGradientState() {
        if (!stateStore.hasPostOrderState || !stateStore.hasPreOrderState) {
            throw new IllegalStateException(
                    "Canonical gradients require both computePostOrderLogLikelihood and computePreOrder to have been called.");
        }
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

    private long transitionCacheMisses(final CanonicalBranchTransitionProvider transitionProvider,
                                       final String phase) {
        if (transitionProvider instanceof CanonicalTransitionCacheDiagnostics) {
            return ((CanonicalTransitionCacheDiagnostics) transitionProvider)
                    .getTransitionCacheMissCount(phase);
        }
        return 0L;
    }

    private void assertNoGradientTransitionMisses(final CanonicalBranchTransitionProvider transitionProvider,
                                                  final String phase,
                                                  final long missesBefore) {
        if (!debugOptions.isAssertNoGradientCacheMissesEnabled()) {
            return;
        }
        final long missesAfter = transitionCacheMisses(transitionProvider, phase);
        if (missesAfter != missesBefore) {
            throw new IllegalStateException(
                    "Canonical transition cache rebuilt " + (missesAfter - missesBefore)
                            + " branch transitions during " + phase + ".");
        }
    }

    private boolean fillLocalAdjointsForBranch(final int childIndex,
                                               final CanonicalBranchTransitionProvider transitionProvider,
                                               final BranchGradientWorkspace workspace) {
        return branchContributionAssembler.fillLocalAdjointsForBranch(childIndex, transitionProvider, workspace);
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
        provider.fillCanonicalTransitionForLength(branchLength, workspace.transition);
        if (tree.isExternal(tree.getNode(childIndex))) {
            final CanonicalTipObservationModel observationModel = stateStore.tipObservationModels[childIndex];
            final CanonicalTransitionMomentProvider momentProvider =
                    observationModel.getMode() == TipObservationMode.PARTIAL_EXACT_IDENTITY
                    ? provider
                    : null;
            observationModel.fillParentMessage(
                    workspace.transition,
                    momentProvider,
                    branchLength,
                    workspace.tipParentMessageWorkspace,
                    workspace.gaussianWorkspace,
                    workspace.state);
        } else {
            CanonicalGaussianMessageOps.pushBackward(
                    stateStore.postOrder[childIndex],
                    workspace.transition,
                    workspace.gaussianWorkspace,
                    workspace.state);
        }
        CanonicalGaussianMessageOps.combineStates(
                stateStore.branchAboveParent[childIndex], workspace.state, workspace.combinedState);
        return CanonicalGaussianMessageOps.normalizationShift(workspace.combinedState, workspace.gaussianWorkspace);
    }

}
