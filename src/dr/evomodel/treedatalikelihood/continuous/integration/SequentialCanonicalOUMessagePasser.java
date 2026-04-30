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
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTransitionCacheDiagnostics;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTransitionMomentProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalRootPrior;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTreeMessagePasser;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
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
    private static final String BRANCH_LENGTH_FD_FALLBACK_PROPERTY =
            "beast.experimental.canonicalBranchLengthGradientFiniteDifference";
    private static final String ASSERT_NO_GRADIENT_CACHE_MISSES_PROPERTY =
            "beast.debug.canonicalAssertNoGradientCacheMisses";

    private final Tree tree;
    private final int dim;
    private final int nodeCount;
    private final int tipCount;
    private final CanonicalTreeStateStore stateStore;
    private final CanonicalTipProjector tipProjector;
    private final CanonicalTreeTraversal treeTraversal;
    private final CanonicalBranchContributionAssembler branchContributionAssembler;
    private final BranchGradientInputs preparedBranchGradientInputs;
    private final BranchGradientWorkspace mainWorkspace;
    private final TaskPool branchGradientTaskPool;
    private final BranchGradientWorkspace[] branchGradientWorkspaces;
    private final CanonicalBranchAdjointPreparer branchAdjointPreparer;
    private final CanonicalTreeGradientEngine treeGradientEngine;
    private final CanonicalBranchLengthGradientEngine branchLengthGradientEngine;
    private final CanonicalLegacyGradientCompatibility legacyGradientCompatibility;

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
        this.branchContributionAssembler = new CanonicalBranchContributionAssembler(tree, dim, stateStore);
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
                branchContributionAssembler::fillLocalAdjointsForBranch);
        this.treeGradientEngine = new CanonicalTreeGradientEngine(
                dim,
                Math.max(0, nodeCount - 1),
                branchGradientTaskPool,
                mainWorkspace,
                branchGradientWorkspaces);
        this.branchLengthGradientEngine = new CanonicalBranchLengthGradientEngine();
        this.legacyGradientCompatibility = new CanonicalLegacyGradientCompatibility(
                tree,
                dim,
                stateStore,
                mainWorkspace,
                branchContributionAssembler);
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
        ensureGradientState();

        final String previousPhase = pushTransitionCachePhase(transitionProvider, PHASE_GRADIENT_PREP);
        try {
            legacyGradientCompatibility.computeGradientQ(transitionProvider, gradQ);
        } finally {
            popTransitionCachePhase(transitionProvider, previousPhase);
        }
    }

    @Override
    public void computeGradientBranchLengths(final CanonicalBranchTransitionProvider transitionProvider,
                                             final double[] gradT) {
        ensureGradientState();

        final String previousPhase = pushTransitionCachePhase(transitionProvider, PHASE_BRANCH_LENGTH_GRADIENT);
        try {
            final CanonicalOUTransitionProvider ouProvider = requireOUProvider(transitionProvider);
            if (!Boolean.getBoolean(BRANCH_LENGTH_FD_FALLBACK_PROPERTY)) {
                final long missesBefore =
                        transitionCacheMisses(transitionProvider, PHASE_BRANCH_LENGTH_GRADIENT);
                branchAdjointPreparer.prepare(transitionProvider, stateStore, preparedBranchGradientInputs);
                assertNoGradientTransitionMisses(
                        transitionProvider, PHASE_BRANCH_LENGTH_GRADIENT, missesBefore);
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
        ensureGradientState();
        legacyGradientCompatibility.computeGradientA(transitionProvider, gradA);
    }

    @Override
    @Deprecated
    public void computeGradientMu(final CanonicalBranchTransitionProvider transitionProvider, final double[] gradMu) {
        ensureGradientState();
        legacyGradientCompatibility.computeGradientMu(transitionProvider, gradMu);
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
            final long missesBefore = transitionCacheMisses(transitionProvider, PHASE_GRADIENT_PREP);
            branchAdjointPreparer.prepare(transitionProvider, stateStore, out);
            assertNoGradientTransitionMisses(transitionProvider, PHASE_GRADIENT_PREP, missesBefore);
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
        if (!Boolean.getBoolean(ASSERT_NO_GRADIENT_CACHE_MISSES_PROPERTY)) {
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

}
