/*
 * CanonicalTreeGradientEngine.java
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

package dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.SelectionMatrixParameterization;
import dr.evomodel.continuous.ou.canonical.SpecializedCanonicalSelectionParameterization;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalOUProviderSupport;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalOUTransitionProvider;
import dr.evomodel.continuous.canonical.scheduling.ChunkSizeStrategy;
import dr.evomodel.continuous.canonical.scheduling.DimensionWeightedChunkSizeStrategy;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.continuous.canonical.math.MatrixOps;
import dr.evomodel.continuous.canonical.workspace.BranchAdjointWorkspace;
import dr.evomodel.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.continuous.canonical.workspace.GradientPullbackWorkspace;
import dr.util.CanonicalPhaseTimer;
import dr.util.TaskPool;

import java.util.Arrays;

/**
 * Formula-style gradient engine for canonical tree OU branch adjoints.
 *
 * <p>The message passer prepares branch-local adjoints once. This engine owns
 * the target-specific pullbacks from those adjoints to selection, diffusion,
 * and stationary-mean gradients.</p>
 */
public final class CanonicalTreeGradientEngine {

    private final int dimension;
    private final int branchCapacity;
    private final TaskPool taskPool;
    private final BranchGradientWorkspace mainWorkspace;
    private final BranchGradientWorkspace[] workspaces;
    private final ChunkSizeStrategy chunkSizeStrategy;

    public CanonicalTreeGradientEngine(
            final int dimension,
            final int branchCapacity,
            final TaskPool taskPool,
            final BranchGradientWorkspace mainWorkspace,
            final BranchGradientWorkspace[] workspaces) {
        this.dimension = dimension;
        this.branchCapacity = branchCapacity;
        this.taskPool = taskPool;
        this.mainWorkspace = mainWorkspace;
        this.workspaces = workspaces;
        this.chunkSizeStrategy =
                new DimensionWeightedChunkSizeStrategy(dimension, Math.max(1, workspaces.length));
    }

    public void computeJointGradients(final CanonicalBranchTransitionProvider transitionProvider,
                               final BranchGradientInputs inputs,
                               final double[] gradA,
                               final double[] gradQ,
                               final double[] gradMu) {
        Arrays.fill(gradA, 0.0);
        Arrays.fill(gradQ, 0.0);
        Arrays.fill(gradMu, 0.0);
        inputs.checkCompatible(branchCapacity, dimension);

        final CanonicalOUTransitionProvider ouProvider =
                CanonicalOUProviderSupport.requireOUProvider(transitionProvider);
        final OUProcessModel processModel =
                selectGradientProcessModel(ouProvider, inputs);
        final OUProcessModel homogeneousProcessModel =
                ouProvider.getProcessModelCount() == 1 ? processModel : null;
        final CanonicalSelectionGradientPullback selectionPullback =
                CanonicalSelectionGradientPullbacks.create(processModel, dimension, gradA, mainWorkspace);
        if (workspaces.length <= 1 || inputs.getActiveBranchCount() <= 1) {
            CanonicalPhaseTimer.recordGradientShape(inputs.getActiveBranchCount(), 1, false);
            computeSequential(inputs, ouProvider, homogeneousProcessModel, selectionPullback, gradA, gradQ, gradMu);
            return;
        }

        CanonicalPhaseTimer.recordGradientShape(inputs.getActiveBranchCount(), taskPool.getNumThreads(), true);
        computeParallel(inputs, ouProvider, homogeneousProcessModel, selectionPullback, gradA, gradQ, gradMu);
    }

    private void computeSequential(final BranchGradientInputs inputs,
                                   final CanonicalOUTransitionProvider ouProvider,
                                   final OUProcessModel homogeneousProcessModel,
                                   final CanonicalSelectionGradientPullback selectionPullback,
                                   final double[] gradA,
                                   final double[] gradQ,
                                   final double[] gradMu) {
        final BranchGradientWorkspace workspace = mainWorkspace;
        selectionPullback.initialize(workspace, gradA, gradMu);

        long timingStart = CanonicalPhaseTimer.start();
        for (int activeIndex = 0; activeIndex < inputs.getActiveBranchCount(); ++activeIndex) {
            final OUProcessModel processModel =
                    processModelFor(ouProvider, homogeneousProcessModel, inputs, activeIndex);
            selectionPullback.accumulateForBranch(
                    processModel,
                    inputs,
                    activeIndex,
                    inputs.getLocalAdjoints(activeIndex),
                    workspace,
                    gradA,
                    gradQ,
                    gradMu);
        }
        CanonicalPhaseTimer.recordBranchGradient(timingStart);

        timingStart = CanonicalPhaseTimer.start();
        finalizeJointGradients(selectionPullback, inputs, workspace, gradA, gradQ, gradMu);
        CanonicalPhaseTimer.recordFinalize(timingStart);
    }

    private void computeParallel(final BranchGradientInputs inputs,
                                 final CanonicalOUTransitionProvider ouProvider,
                                 final OUProcessModel homogeneousProcessModel,
                                 final CanonicalSelectionGradientPullback selectionPullback,
                                 final double[] gradA,
                                 final double[] gradQ,
                                 final double[] gradMu) {
        final BranchGradientWorkspace reductionWorkspace = mainWorkspace;

        for (int worker = 0; worker < workspaces.length; ++worker) {
            final BranchGradientWorkspace workspace = workspaces[worker];
            selectionPullback.clearWorkerBuffers(workspace, gradA.length, gradMu.length);
            selectionPullback.prepareWorkspace(workspace);
        }

        selectionPullback.initialize(reductionWorkspace, gradA, gradMu);

        long timingStart = CanonicalPhaseTimer.start();
        taskPool.forkDynamicBalanced(
                inputs.getActiveBranchCount(),
                branchGradientJointChunkSize(inputs.getActiveBranchCount()),
                (activeIndex, thread) -> {
                    final BranchGradientWorkspace workspace = workspaces[thread];
                    final OUProcessModel processModel =
                            processModelFor(ouProvider, homogeneousProcessModel, inputs, activeIndex);

                    selectionPullback.accumulateForBranch(
                            processModel,
                            inputs,
                            activeIndex,
                            inputs.getLocalAdjoints(activeIndex),
                            workspace,
                            workspace.localGradientA,
                            workspace.localGradientQ,
                            workspace.localGradientMu(gradMu.length, dimension));
                });
        CanonicalPhaseTimer.recordBranchGradient(timingStart);

        timingStart = CanonicalPhaseTimer.start();
        for (int worker = 0; worker < workspaces.length; ++worker) {
            final BranchGradientWorkspace workspace = workspaces[worker];
            accumulateVectorInPlace(gradQ, workspace.localGradientQ, gradQ.length);
            accumulateVectorInPlace(gradMu, workspace.localGradientMu(gradMu.length, dimension), gradMu.length);
            selectionPullback.reduceWorker(workspace, reductionWorkspace, gradA);
        }
        CanonicalPhaseTimer.recordReduction(timingStart);

        timingStart = CanonicalPhaseTimer.start();
        finalizeJointGradients(selectionPullback, inputs, reductionWorkspace, gradA, gradQ, gradMu);
        CanonicalPhaseTimer.recordFinalize(timingStart);
    }

    private void finalizeJointGradients(final CanonicalSelectionGradientPullback selectionPullback,
                                        final BranchGradientInputs inputs,
                                        final BranchGradientWorkspace workspace,
                                        final double[] gradA,
                                        final double[] gradQ,
                                        final double[] gradMu) {
        selectionPullback.finish(inputs, workspace, gradA, gradQ, gradMu);

        if (inputs.getRootDiffusionScale() > 0.0) {
            accumulateRootDiffusionGradient(
                    inputs.getRootPreOrderState(),
                    inputs.getRootPostOrderState(),
                    inputs.getRootDiffusionScale(),
                    gradQ,
                    workspace);
        }
    }

    private int branchGradientJointChunkSize(final int taskLimit) {
        return chunkSizeStrategy.chunkSize(taskLimit);
    }

    private OUProcessModel selectGradientProcessModel(final CanonicalOUTransitionProvider ouProvider,
                                                      final BranchGradientInputs inputs) {
        final int activeBranchCount = inputs.getActiveBranchCount();
        if (activeBranchCount == 0) {
            return ouProvider.getProcessModel(0);
        }
        final int firstProcessIndex = inputs.getProcessIndex(0);
        final OUProcessModel baseProcessModel = ouProvider.getProcessModel(firstProcessIndex);
        for (int activeIndex = 1; activeIndex < activeBranchCount; ++activeIndex) {
            final int processIndex = inputs.getProcessIndex(activeIndex);
            validateSharedGradientTarget(
                    baseProcessModel,
                    ouProvider.getProcessModel(processIndex),
                    processIndex);
        }
        return baseProcessModel;
    }

    private OUProcessModel processModelFor(final CanonicalOUTransitionProvider ouProvider,
                                           final OUProcessModel homogeneousProcessModel,
                                           final BranchGradientInputs inputs,
                                           final int activeIndex) {
        return homogeneousProcessModel == null
                ? ouProvider.getProcessModel(inputs.getProcessIndex(activeIndex))
                : homogeneousProcessModel;
    }

    private static void validateSharedGradientTarget(final OUProcessModel baseProcessModel,
                                                     final OUProcessModel candidateProcessModel,
                                                     final int processIndex) {
        if (candidateProcessModel == baseProcessModel) {
            return;
        }
        final SelectionMatrixParameterization baseSelection =
                baseProcessModel.getSelectionMatrixParameterization();
        final SelectionMatrixParameterization candidateSelection =
                candidateProcessModel.getSelectionMatrixParameterization();
        if (baseSelection.getClass() != candidateSelection.getClass()
                || baseSelection.getMatrixParameter() != candidateSelection.getMatrixParameter()
                || baseProcessModel.getDiffusionMatrix() != candidateProcessModel.getDiffusionMatrix()
                || baseProcessModel.getStationaryMeanParameter()
                != candidateProcessModel.getStationaryMeanParameter()) {
            throw new UnsupportedOperationException(
                    "Canonical OU gradients for heterogeneous processes require active branch processes "
                            + "to share selection, diffusion, and stationary-mean parameters; process "
                            + processIndex + " does not match the first active process.");
        }
        if (baseSelection instanceof SpecializedCanonicalSelectionParameterization
                && baseSelection != candidateSelection) {
            throw new UnsupportedOperationException(
                    "Canonical OU specialized selection gradients for heterogeneous processes require "
                            + "the same specialized selection parameterization instance for all active branches.");
        }
    }

    private void accumulateRootDiffusionGradient(final CanonicalGaussianState rootPreOrder,
                                                 final CanonicalGaussianState rootPostOrder,
                                                 final double rootDiffusionScale,
                                                 final double[] gradQ,
                                                 final BranchGradientWorkspace workspace) {
        final BranchAdjointWorkspace adjoint = workspace.adjoint;
        final GradientPullbackWorkspace gradient = workspace.gradient;

        CanonicalGaussianMessageOps.combineStates(rootPreOrder, rootPostOrder, adjoint.combinedState);

        adjoint.fillMomentsFromCanonical(
                adjoint.combinedState, adjoint.mean, adjoint.covarianceFlat, dimension);
        adjoint.fillMomentsFromCanonical(
                rootPreOrder, adjoint.mean2, gradient.covariance2, dimension);

        final double[] priorPrecision = rootPreOrder.precision;
        for (int i = 0; i < dimension; ++i) {
            final double deltaI = adjoint.mean[i] - adjoint.mean2[i];
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                gradient.covarianceAdjointFlat[iOff + j] =
                        adjoint.covarianceFlat[iOff + j] + deltaI * (adjoint.mean[j] - adjoint.mean2[j]);
            }
        }

        MatrixOps.matMul(
                priorPrecision, gradient.covarianceAdjointFlat, gradient.matrixProductFlat, dimension);
        MatrixOps.matMul(
                gradient.matrixProductFlat, priorPrecision, gradient.covarianceAdjointFlat, dimension);
        for (int i = 0; i < dimension; ++i) {
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                final int ij = iOff + j;
                gradQ[iOff + j] += rootDiffusionScale
                        * (-0.5 * priorPrecision[ij] + 0.5 * gradient.covarianceAdjointFlat[ij]);
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

}
