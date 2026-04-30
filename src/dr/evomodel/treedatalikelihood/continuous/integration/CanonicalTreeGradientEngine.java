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
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockBranchGradientWorkspace;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockPreparedBranchBasis;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalOUTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianUtils;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.GaussianMatrixOps;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
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

        final CanonicalOUTransitionProvider ouProvider = requireOUProvider(transitionProvider);
        final OUProcessModel processModel = ouProvider.getProcessModel();
        final OrthogonalBlockCanonicalParameterization orthogonalSelection =
                processModel.getSelectionMatrixParameterization()
                        instanceof OrthogonalBlockCanonicalParameterization
                        ? (OrthogonalBlockCanonicalParameterization)
                        processModel.getSelectionMatrixParameterization()
                        : null;
        if (workspaces.length <= 1 || inputs.getActiveBranchCount() <= 1) {
            computeSequential(inputs, processModel, orthogonalSelection, gradA, gradQ, gradMu);
            return;
        }

        computeParallel(inputs, processModel, orthogonalSelection, gradA, gradQ, gradMu);
    }

    private void computeSequential(final BranchGradientInputs inputs,
                                   final OUProcessModel processModel,
                                   final OrthogonalBlockCanonicalParameterization orthogonalSelection,
                                   final double[] gradA,
                                   final double[] gradQ,
                                   final double[] gradMu) {
        final BranchGradientWorkspace workspace = mainWorkspace;
        final OrthogonalGradientLayout orthogonalLayout =
                orthogonalSelection == null ? null : validateOrthogonalGradientLayout(orthogonalSelection, gradA, workspace);
        final GradientPullbackWorkspace gradient = workspace.gradient;

        if (orthogonalLayout != null) {
            Arrays.fill(gradient.orthogonalCompressedGradientScratch, 0, orthogonalLayout.compressedBlockDim, 0.0);
            Arrays.fill(gradient.orthogonalNativeGradientScratch, 0, orthogonalLayout.nativeBlockDim, 0.0);
            Arrays.fill(gradient.orthogonalRotationGradientFlatScratch, 0.0);
            workspace.ensureOrthogonalBranchWorkspace(orthogonalSelection);
        }

        for (int activeIndex = 0; activeIndex < inputs.getActiveBranchCount(); ++activeIndex) {
            copyAdjoints(inputs.getLocalAdjoints(activeIndex), workspace.adjoints);
            accumulateForBranch(
                    inputs,
                    processModel,
                    orthogonalSelection,
                    activeIndex,
                    workspace,
                    gradA,
                    gradQ,
                    gradMu);
        }

        finalizeJointGradients(orthogonalLayout, inputs, workspace, gradA, gradQ);
    }

    private void computeParallel(final BranchGradientInputs inputs,
                                 final OUProcessModel processModel,
                                 final OrthogonalBlockCanonicalParameterization orthogonalSelection,
                                 final double[] gradA,
                                 final double[] gradQ,
                                 final double[] gradMu) {
        final BranchGradientWorkspace reductionWorkspace = mainWorkspace;
        final OrthogonalGradientLayout orthogonalLayout =
                orthogonalSelection == null ? null : validateOrthogonalGradientLayout(orthogonalSelection, gradA, reductionWorkspace);

        for (int worker = 0; worker < workspaces.length; ++worker) {
            final BranchGradientWorkspace workspace = workspaces[worker];
            workspace.clearLocalGradientBuffers(
                    gradA.length,
                    gradMu.length,
                    dimension,
                    orthogonalLayout != null,
                    orthogonalLayout == null ? 0 : orthogonalLayout.compressedBlockDim);
            if (orthogonalLayout != null) {
                workspace.ensureOrthogonalBranchWorkspace(orthogonalSelection);
            }
        }

        if (orthogonalLayout != null) {
            final GradientPullbackWorkspace gradient = reductionWorkspace.gradient;
            Arrays.fill(gradient.orthogonalCompressedGradientScratch, 0, orthogonalLayout.compressedBlockDim, 0.0);
            Arrays.fill(gradient.orthogonalNativeGradientScratch, 0, orthogonalLayout.nativeBlockDim, 0.0);
            Arrays.fill(gradient.orthogonalRotationGradientFlatScratch, 0.0);
        }

        taskPool.forkDynamicBalanced(
                inputs.getActiveBranchCount(),
                branchGradientJointChunkSize(inputs.getActiveBranchCount()),
                (activeIndex, thread) -> {
                    final BranchGradientWorkspace workspace = workspaces[thread];
                    copyAdjoints(inputs.getLocalAdjoints(activeIndex), workspace.adjoints);

                    accumulateForBranch(
                            inputs,
                            processModel,
                            orthogonalSelection,
                            activeIndex,
                            workspace,
                            workspace.localGradientA,
                            workspace.localGradientQ,
                            workspace.localGradientMu(gradMu.length, dimension));
                });

        for (int worker = 0; worker < workspaces.length; ++worker) {
            final BranchGradientWorkspace workspace = workspaces[worker];
            accumulateVectorInPlace(gradQ, workspace.localGradientQ, gradQ.length);
            accumulateVectorInPlace(gradMu, workspace.localGradientMu(gradMu.length, dimension), gradMu.length);
            if (orthogonalLayout != null) {
                accumulateVectorInPlace(
                        reductionWorkspace.gradient.orthogonalCompressedGradientScratch,
                        workspace.gradient.orthogonalCompressedGradientScratch,
                        orthogonalLayout.compressedBlockDim);
                accumulateVectorInPlace(
                        reductionWorkspace.gradient.orthogonalRotationGradientFlatScratch,
                        workspace.gradient.orthogonalRotationGradientFlatScratch,
                        reductionWorkspace.gradient.orthogonalRotationGradientFlatScratch.length);
            } else {
                accumulateVectorInPlace(gradA, workspace.localGradientA, gradA.length);
            }
        }

        finalizeJointGradients(orthogonalLayout, inputs, reductionWorkspace, gradA, gradQ);
    }

    private void accumulateForBranch(
            final BranchGradientInputs inputs,
            final OUProcessModel processModel,
            final OrthogonalBlockCanonicalParameterization orthogonalSelection,
            final int activeIndex,
            final BranchGradientWorkspace workspace,
            final double[] gradA,
            final double[] gradQ,
            final double[] gradMu) {
        if (orthogonalSelection != null) {
            accumulateOrthogonalForBranch(
                    inputs.getOrthogonalPreparedBasis(activeIndex),
                    processModel,
                    orthogonalSelection,
                    workspace,
                    gradQ,
                    gradMu);
            return;
        }

        accumulateDenseForBranch(
                processModel,
                workspace,
                inputs.getBranchLength(activeIndex),
                gradA,
                gradQ,
                gradMu);
    }

    private void accumulateOrthogonalForBranch(
            final OrthogonalBlockPreparedBranchBasis preparedBasis,
            final OUProcessModel processModel,
            final OrthogonalBlockCanonicalParameterization orthogonalSelection,
            final BranchGradientWorkspace workspace,
            final double[] gradQ,
            final double[] gradMu) {
        final OrthogonalBlockBranchGradientWorkspace orthogonalWorkspace =
                workspace.ensureOrthogonalBranchWorkspace(orthogonalSelection);
        final GradientPullbackWorkspace gradient = workspace.gradient;
        orthogonalSelection.accumulateNativeGradientFromAdjointsPreparedFlat(
                preparedBasis,
                processModel.getDiffusionMatrix(),
                workspace.adjoints,
                orthogonalWorkspace,
                gradient.orthogonalCompressedGradientScratch,
                gradient.orthogonalRotationGradientFlatScratch);

        orthogonalSelection.accumulateDiffusionGradientPreparedFlat(
                preparedBasis,
                workspace.adjoints.dLogL_dOmega,
                true,
                gradQ,
                orthogonalWorkspace);

        orthogonalSelection.accumulateMeanGradientPrepared(
                preparedBasis,
                workspace.adjoints.dLogL_df,
                gradMu,
                orthogonalWorkspace);
    }

    private void accumulateDenseForBranch(final OUProcessModel processModel,
                                          final BranchGradientWorkspace workspace,
                                          final double branchLength,
                                          final double[] gradA,
                                          final double[] gradQ,
                                          final double[] gradMu) {
        final GradientPullbackWorkspace gradient = workspace.gradient;

        processModel.accumulateSelectionGradientFlat(
                branchLength,
                workspace.adjoints.dLogL_dF,
                workspace.adjoints.dLogL_df,
                gradA);

        processModel.accumulateSelectionGradientFromCovarianceFlat(
                branchLength,
                workspace.adjoints.dLogL_dOmega,
                true,
                gradA);

        processModel.accumulateDiffusionGradientFlat(
                branchLength,
                workspace.adjoints.dLogL_dOmega,
                false,
                gradQ);

        processModel.fillTransitionMatrixFlat(branchLength, gradient.transitionMatrixFlat);
        accumulateStationaryMeanGradientFlat(
                gradient.transitionMatrixFlat,
                workspace.adjoints.dLogL_df,
                gradMu);
    }

    private void finalizeJointGradients(final OrthogonalGradientLayout orthogonalLayout,
                                        final BranchGradientInputs inputs,
                                        final BranchGradientWorkspace workspace,
                                        final double[] gradA,
                                        final double[] gradQ) {
        if (orthogonalLayout != null) {
            final GradientPullbackWorkspace gradient = workspace.gradient;
            orthogonalLayout.blockParameter.chainGradient(
                    gradient.orthogonalCompressedGradientScratch,
                    gradient.orthogonalNativeGradientScratch);
            final double[] angleGradient =
                    orthogonalLayout.orthogonalRotation.pullBackGradientFlat(
                            gradient.orthogonalRotationGradientFlatScratch,
                            dimension);
            System.arraycopy(gradient.orthogonalNativeGradientScratch, 0, gradA, 0, orthogonalLayout.nativeBlockDim);
            System.arraycopy(angleGradient, 0, gradA, orthogonalLayout.nativeBlockDim, angleGradient.length);
        } else {
            GaussianMatrixOps.transposeFlatSquareInPlace(gradA, dimension);
        }

        if (inputs.getRootDiffusionScale() > 0.0) {
            accumulateRootDiffusionGradient(
                    inputs.getRootPreOrderState(),
                    inputs.getRootPostOrderState(),
                    inputs.getRootDiffusionScale(),
                    gradQ,
                    workspace);
        }
    }

    private OrthogonalGradientLayout validateOrthogonalGradientLayout(
            final OrthogonalBlockCanonicalParameterization orthogonalSelection,
            final double[] gradA,
            final BranchGradientWorkspace workspace) {
        final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter =
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) orthogonalSelection.getMatrixParameter();
        if (!(blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider)) {
            throw new IllegalStateException(
                    "Orthogonal block native gradient requires an OrthogonalMatrixProvider rotation parameter.");
        }

        final OrthogonalMatrixProvider orthogonalRotation =
                (OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter();
        final int nativeBlockDim = blockParameter.getBlockDiagonalNParameters();
        final int angleDim = orthogonalRotation.getOrthogonalParameter().getDimension();
        final int nativeDim = nativeBlockDim + angleDim;
        if (gradA.length != nativeDim) {
            throw new IllegalArgumentException(
                    "Orthogonal block selection gradient expects native parameter length "
                            + nativeDim + ", found " + gradA.length);
        }

        final int compressedBlockDim = blockParameter.getCompressedDDimension();
        final GradientPullbackWorkspace gradient = workspace.gradient;
        if (compressedBlockDim > gradient.orthogonalCompressedGradientScratch.length
                || nativeBlockDim > gradient.orthogonalNativeGradientScratch.length) {
            throw new IllegalStateException(
                    "Orthogonal block scratch is too small for native gradient dimensions "
                            + compressedBlockDim + " and " + nativeBlockDim + ".");
        }
        return new OrthogonalGradientLayout(
                blockParameter,
                orthogonalRotation,
                compressedBlockDim,
                nativeBlockDim);
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

    private void accumulateStationaryMeanGradientFlat(final double[] transitionMatrix,
                                                      final double[] adjointB,
                                                      final double[] gradient) {
        if (gradient.length == 1) {
            double sum = 0.0;
            for (int i = 0; i < dimension; ++i) {
                double ftAdjoint = 0.0;
                for (int j = 0; j < dimension; ++j) {
                    ftAdjoint += transitionMatrix[j * dimension + i] * adjointB[j];
                }
                sum += adjointB[i] - ftAdjoint;
            }
            gradient[0] += sum;
            return;
        }
        if (gradient.length != dimension) {
            throw new IllegalArgumentException(
                    "Stationary-mean gradient length must be 1 or " + dimension + ", found " + gradient.length);
        }

        for (int i = 0; i < dimension; ++i) {
            double ftAdjoint = 0.0;
            for (int j = 0; j < dimension; ++j) {
                ftAdjoint += transitionMatrix[j * dimension + i] * adjointB[j];
            }
            gradient[i] += adjointB[i] - ftAdjoint;
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

        CanonicalGaussianUtils.fillMomentsFromCanonical(
                adjoint.combinedState, adjoint.mean, adjoint.covariance);
        CanonicalGaussianUtils.fillMomentsFromCanonical(
                rootPreOrder, adjoint.mean2, gradient.covariance2);

        final double[] priorPrecision = rootPreOrder.precision;
        for (int i = 0; i < dimension; ++i) {
            final double deltaI = adjoint.mean[i] - adjoint.mean2[i];
            for (int j = 0; j < dimension; ++j) {
                gradient.covarianceAdjoint[i][j] =
                        adjoint.covariance[i][j] + deltaI * (adjoint.mean[j] - adjoint.mean2[j]);
            }
        }

        GaussianMatrixOps.multiplyFlatByMatrix(
                priorPrecision, gradient.covarianceAdjoint, gradient.transitionMatrix, dimension);
        GaussianMatrixOps.multiplyMatrixByFlat(
                gradient.transitionMatrix, priorPrecision, gradient.covarianceAdjoint, dimension);
        for (int i = 0; i < dimension; ++i) {
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                gradQ[iOff + j] += rootDiffusionScale
                        * (-0.5 * priorPrecision[iOff + j] + 0.5 * gradient.covarianceAdjoint[i][j]);
            }
        }
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

    private static void copyAdjoints(final CanonicalLocalTransitionAdjoints source,
                                     final CanonicalLocalTransitionAdjoints target) {
        System.arraycopy(source.dLogL_dF, 0, target.dLogL_dF, 0, source.dLogL_dF.length);
        System.arraycopy(source.dLogL_df, 0, target.dLogL_df, 0, source.dLogL_df.length);
        System.arraycopy(source.dLogL_dOmega, 0, target.dLogL_dOmega, 0, source.dLogL_dOmega.length);
    }

    private static void accumulateVectorInPlace(final double[] target,
                                                final double[] source,
                                                final int length) {
        for (int i = 0; i < length; ++i) {
            target[i] += source[i];
        }
    }

    private static final class OrthogonalGradientLayout {
        final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
        final OrthogonalMatrixProvider orthogonalRotation;
        final int compressedBlockDim;
        final int nativeBlockDim;

        private OrthogonalGradientLayout(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                         final OrthogonalMatrixProvider orthogonalRotation,
                                         final int compressedBlockDim,
                                         final int nativeBlockDim) {
            this.blockParameter = blockParameter;
            this.orthogonalRotation = orthogonalRotation;
            this.compressedBlockDim = compressedBlockDim;
            this.nativeBlockDim = nativeBlockDim;
        }
    }
}
