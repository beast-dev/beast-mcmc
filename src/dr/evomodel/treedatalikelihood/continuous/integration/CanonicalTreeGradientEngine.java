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

package dr.evomodel.treedatalikelihood.continuous.integration;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalOUTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianUtils;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.GaussianMatrixOps;
import dr.util.TaskPool;

import java.util.Arrays;

/**
 * Formula-style gradient engine for canonical tree OU branch adjoints.
 *
 * <p>The message passer prepares branch-local adjoints once. This engine owns
 * the target-specific pullbacks from those adjoints to selection, diffusion,
 * and stationary-mean gradients.</p>
 */
final class CanonicalTreeGradientEngine {

    private final int dimension;
    private final int branchCapacity;
    private final TaskPool taskPool;
    private final BranchGradientWorkspace mainWorkspace;
    private final BranchGradientWorkspace[] workspaces;

    CanonicalTreeGradientEngine(
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
    }

    void computeJointGradients(final CanonicalBranchTransitionProvider transitionProvider,
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
        final OUProcessModel processModel = ouProvider.getProcessModel();
        final CanonicalSelectionGradientPullback selectionPullback =
                CanonicalSelectionGradientPullbacks.create(processModel, dimension, gradA, mainWorkspace);
        if (workspaces.length <= 1 || inputs.getActiveBranchCount() <= 1) {
            computeSequential(inputs, processModel, selectionPullback, gradA, gradQ, gradMu);
            return;
        }

        computeParallel(inputs, processModel, selectionPullback, gradA, gradQ, gradMu);
    }

    private void computeSequential(final BranchGradientInputs inputs,
                                   final OUProcessModel processModel,
                                   final CanonicalSelectionGradientPullback selectionPullback,
                                   final double[] gradA,
                                   final double[] gradQ,
                                   final double[] gradMu) {
        final BranchGradientWorkspace workspace = mainWorkspace;
        selectionPullback.initialize(workspace, gradA, gradMu);

        for (int activeIndex = 0; activeIndex < inputs.getActiveBranchCount(); ++activeIndex) {
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

        finalizeJointGradients(selectionPullback, inputs, workspace, gradA, gradQ);
    }

    private void computeParallel(final BranchGradientInputs inputs,
                                 final OUProcessModel processModel,
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

        taskPool.forkDynamicBalanced(
                inputs.getActiveBranchCount(),
                branchGradientJointChunkSize(inputs.getActiveBranchCount()),
                (activeIndex, thread) -> {
                    final BranchGradientWorkspace workspace = workspaces[thread];

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

        for (int worker = 0; worker < workspaces.length; ++worker) {
            final BranchGradientWorkspace workspace = workspaces[worker];
            accumulateVectorInPlace(gradQ, workspace.localGradientQ, gradQ.length);
            accumulateVectorInPlace(gradMu, workspace.localGradientMu(gradMu.length, dimension), gradMu.length);
            selectionPullback.reduceWorker(workspace, reductionWorkspace, gradA);
        }

        finalizeJointGradients(selectionPullback, inputs, reductionWorkspace, gradA, gradQ);
    }

    private void finalizeJointGradients(final CanonicalSelectionGradientPullback selectionPullback,
                                        final BranchGradientInputs inputs,
                                        final BranchGradientWorkspace workspace,
                                        final double[] gradA,
                                        final double[] gradQ) {
        selectionPullback.finish(inputs, workspace, gradA);

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
        return branchGradientChunkSize(taskLimit, 3, 8);
    }

    private int branchGradientChunkSize(final int taskLimit,
                                        final int targetChunksPerWorker,
                                        final int maxChunkSize) {
        final int workerCount = Math.max(1, workspaces.length);
        final int suggested =
                (taskLimit + workerCount * targetChunksPerWorker - 1) / (workerCount * targetChunksPerWorker);
        return Math.max(1, Math.min(maxChunkSize, suggested));
    }

    private void accumulateRootDiffusionGradient(final CanonicalGaussianState rootPreOrder,
                                                 final CanonicalGaussianState rootPostOrder,
                                                 final double rootDiffusionScale,
                                                 final double[] gradQ,
                                                 final BranchGradientWorkspace workspace) {
        final BranchAdjointWorkspace adjoint = workspace.adjoint;
        final GradientPullbackWorkspace gradient = workspace.gradient;

        CanonicalGaussianMessageOps.combineStates(rootPreOrder, rootPostOrder, adjoint.combinedState);

        CanonicalGaussianUtils.fillMomentsFromCanonical(
                adjoint.combinedState, adjoint.mean, adjoint.covariance);
        CanonicalGaussianUtils.fillMomentsFromCanonical(
                rootPreOrder, adjoint.mean2, gradient.covariance2);

        final double[] priorPrecision = rootPreOrder.precision;
        for (int i = 0; i < dimension; ++i) {
            final double deltaI = adjoint.mean[i] - adjoint.mean2[i];
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                gradient.covarianceAdjointFlat[iOff + j] =
                        adjoint.covariance[i][j] + deltaI * (adjoint.mean[j] - adjoint.mean2[j]);
            }
        }

        GaussianMatrixOps.multiplyMatrixMatrixFlat(
                priorPrecision, gradient.covarianceAdjointFlat, gradient.matrixProductFlat, dimension);
        GaussianMatrixOps.multiplyMatrixMatrixFlat(
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
