package dr.evomodel.treedatalikelihood.continuous.integration;

import dr.evolution.tree.Tree;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalOUTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianUtils;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.GaussianMatrixOps;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;

import java.util.Arrays;

/**
 * Deprecated single-target gradient entry points kept for older tests/callers.
 */
final class CanonicalLegacyGradientCompatibility {

    private final Tree tree;
    private final int dim;
    private final int nodeCount;
    private final CanonicalTreeStateStore stateStore;
    private final BranchGradientWorkspace workspace;
    private final CanonicalBranchContributionAssembler branchContributionAssembler;

    CanonicalLegacyGradientCompatibility(final Tree tree,
                                         final int dim,
                                         final CanonicalTreeStateStore stateStore,
                                         final BranchGradientWorkspace workspace,
                                         final CanonicalBranchContributionAssembler branchContributionAssembler) {
        this.tree = tree;
        this.dim = dim;
        this.nodeCount = tree.getNodeCount();
        this.stateStore = stateStore;
        this.workspace = workspace;
        this.branchContributionAssembler = branchContributionAssembler;
    }

    void computeGradientQ(final CanonicalBranchTransitionProvider transitionProvider, final double[] gradQ) {
        Arrays.fill(gradQ, 0.0);

        final CanonicalOUTransitionProvider ouProvider =
                CanonicalOUProviderSupport.requireOUProvider(transitionProvider);
        final OUProcessModel processModel = ouProvider.getProcessModel();
        final OrthogonalBlockCanonicalParameterization orthogonalSelection =
                processModel.getSelectionMatrixParameterization()
                        instanceof OrthogonalBlockCanonicalParameterization
                        ? (OrthogonalBlockCanonicalParameterization)
                        processModel.getSelectionMatrixParameterization()
                        : null;
        final int rootIndex = tree.getRoot().getNumber();

        for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
            if (childIndex == rootIndex) {
                continue;
            }
            if (!branchContributionAssembler.fillLocalAdjointsForBranch(childIndex, transitionProvider, workspace)) {
                continue;
            }
            if (orthogonalSelection != null) {
                GaussianMatrixOps.transposeFlatToMatrix(
                        workspace.adjoints.dLogL_dOmega, workspace.covarianceAdjoint, dim);
                orthogonalSelection.accumulateDiffusionGradient(
                        processModel.getDiffusionMatrix(),
                        transitionProvider.getEffectiveBranchLength(childIndex),
                        workspace.covarianceAdjoint,
                        gradQ);
            } else {
                GaussianMatrixOps.copyFlatToMatrix(workspace.adjoints.dLogL_dOmega, workspace.covarianceAdjoint, dim);
                processModel.accumulateDiffusionGradient(
                        transitionProvider.getEffectiveBranchLength(childIndex),
                        workspace.covarianceAdjoint,
                        gradQ);
            }
        }

        if (stateStore.lastRootDiffusionScale > 0.0) {
            accumulateRootDiffusionGradient(gradQ);
        }
    }

    void computeGradientA(final CanonicalBranchTransitionProvider transitionProvider, final double[] gradA) {
        Arrays.fill(gradA, 0.0);

        final CanonicalOUTransitionProvider ouProvider =
                CanonicalOUProviderSupport.requireOUProvider(transitionProvider);
        final OUProcessModel processModel = ouProvider.getProcessModel();
        if (processModel.getSelectionMatrixParameterization()
                instanceof OrthogonalBlockCanonicalParameterization) {
            computeOrthogonalBlockGradientA(
                    transitionProvider,
                    processModel,
                    (OrthogonalBlockCanonicalParameterization)
                            processModel.getSelectionMatrixParameterization(),
                    gradA);
            return;
        }

        final int rootIndex = tree.getRoot().getNumber();
        for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
            if (childIndex == rootIndex) {
                continue;
            }
            if (!branchContributionAssembler.fillLocalAdjointsForBranch(childIndex, transitionProvider, workspace)) {
                continue;
            }

            final double branchLength = transitionProvider.getEffectiveBranchLength(childIndex);
            GaussianMatrixOps.copyFlatToMatrix(workspace.adjoints.dLogL_dF, workspace.transitionMatrix, dim);
            processModel.accumulateSelectionGradient(
                    branchLength,
                    workspace.transitionMatrix,
                    workspace.adjoints.dLogL_df,
                    gradA);

            GaussianMatrixOps.transposeFlatToMatrix(
                    workspace.adjoints.dLogL_dOmega, workspace.covarianceAdjoint, dim);
            processModel.accumulateSelectionGradientFromCovariance(
                    branchLength,
                    workspace.covarianceAdjoint,
                    gradA);
        }

        GaussianMatrixOps.transposeFlatSquareInPlace(gradA, dim);
    }

    void computeGradientMu(final CanonicalBranchTransitionProvider transitionProvider, final double[] gradMu) {
        Arrays.fill(gradMu, 0.0);

        final CanonicalOUTransitionProvider ouProvider =
                CanonicalOUProviderSupport.requireOUProvider(transitionProvider);
        final OUProcessModel processModel = ouProvider.getProcessModel();
        final int rootIndex = tree.getRoot().getNumber();

        for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
            if (childIndex == rootIndex) {
                continue;
            }
            if (!branchContributionAssembler.fillLocalAdjointsForBranch(childIndex, transitionProvider, workspace)) {
                continue;
            }

            final double branchLength = transitionProvider.getEffectiveBranchLength(childIndex);
            processModel.fillTransitionMatrix(branchLength, workspace.transitionMatrix);
            accumulateStationaryMeanGradient(
                    workspace.transitionMatrix,
                    workspace.adjoints.dLogL_df,
                    gradMu);
        }
    }

    private void computeOrthogonalBlockGradientA(final CanonicalBranchTransitionProvider transitionProvider,
                                                 final OUProcessModel processModel,
                                                 final OrthogonalBlockCanonicalParameterization parameterization,
                                                 final double[] gradA) {
        final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter =
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) parameterization.getMatrixParameter();
        if (!(blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider)) {
            throw new IllegalStateException(
                    "Orthogonal block native gradient requires an OrthogonalMatrixProvider rotation parameter.");
        }

        final int compressedBlockDim = blockParameter.getCompressedDDimension();
        final int nativeBlockDim = blockParameter.getBlockDiagonalNParameters();
        final OrthogonalMatrixProvider orthogonalRotation =
                (OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter();
        final int angleDim = orthogonalRotation.getOrthogonalParameter().getDimension();
        final int nativeDim = nativeBlockDim + angleDim;
        if (gradA.length != nativeDim) {
            throw new IllegalArgumentException(
                    "Orthogonal block selection gradient expects native parameter length "
                            + nativeDim + ", found " + gradA.length);
        }

        if (compressedBlockDim > workspace.orthogonalCompressedGradientScratch.length
                || nativeBlockDim > workspace.orthogonalNativeGradientScratch.length) {
            throw new IllegalStateException(
                    "Orthogonal block scratch is too small for native gradient dimensions "
                            + compressedBlockDim + " and " + nativeBlockDim + ".");
        }
        final double[] stationaryMean = workspace.orthogonalStationaryMeanScratch;
        processModel.getInitialMean(stationaryMean);
        final double[] compressedBlockGradient = workspace.orthogonalCompressedGradientScratch;
        Arrays.fill(compressedBlockGradient, 0, compressedBlockDim, 0.0);
        final double[] nativeBlockGradient = workspace.orthogonalNativeGradientScratch;
        Arrays.fill(nativeBlockGradient, 0, nativeBlockDim, 0.0);
        final double[][] rotationGradient = workspace.orthogonalRotationGradientScratch;
        clearSquare(rotationGradient);

        final int rootIndex = tree.getRoot().getNumber();
        for (int childIndex = 0; childIndex < nodeCount; childIndex++) {
            if (childIndex == rootIndex) {
                continue;
            }
            if (!branchContributionAssembler.fillLocalAdjointsForBranch(childIndex, transitionProvider, workspace)) {
                continue;
            }

            parameterization.accumulateNativeGradientFromAdjoints(
                    processModel.getDiffusionMatrix(),
                    stationaryMean,
                    transitionProvider.getEffectiveBranchLength(childIndex),
                    workspace.adjoints,
                    compressedBlockGradient,
                    rotationGradient);
        }

        blockParameter.chainGradient(compressedBlockGradient, nativeBlockGradient);
        final double[] angleGradient = orthogonalRotation.pullBackGradient(rotationGradient);
        System.arraycopy(nativeBlockGradient, 0, gradA, 0, nativeBlockDim);
        System.arraycopy(angleGradient, 0, gradA, nativeBlockDim, angleGradient.length);
    }

    private void accumulateStationaryMeanGradient(final double[][] transitionMatrix,
                                                  final double[] adjointB,
                                                  final double[] gradient) {
        if (gradient.length == 1) {
            double sum = 0.0;
            for (int i = 0; i < dim; ++i) {
                double ftAdjoint = 0.0;
                for (int j = 0; j < dim; ++j) {
                    ftAdjoint += transitionMatrix[j][i] * adjointB[j];
                }
                sum += adjointB[i] - ftAdjoint;
            }
            gradient[0] += sum;
            return;
        }
        if (gradient.length != dim) {
            throw new IllegalArgumentException(
                    "Stationary-mean gradient length must be 1 or " + dim + ", found " + gradient.length);
        }

        for (int i = 0; i < dim; ++i) {
            double ftAdjoint = 0.0;
            for (int j = 0; j < dim; ++j) {
                ftAdjoint += transitionMatrix[j][i] * adjointB[j];
            }
            gradient[i] += adjointB[i] - ftAdjoint;
        }
    }

    private void accumulateRootDiffusionGradient(final double[] gradQ) {
        final int rootIndex = tree.getRoot().getNumber();
        CanonicalGaussianMessageOps.combineStates(
                stateStore.preOrder[rootIndex], stateStore.postOrder[rootIndex], workspace.combinedState);

        CanonicalGaussianUtils.fillMomentsFromCanonical(
                workspace.combinedState, workspace.mean, workspace.covariance);
        CanonicalGaussianUtils.fillMomentsFromCanonical(
                stateStore.preOrder[rootIndex], workspace.mean2, workspace.covariance2);

        final double[] priorPrecision = stateStore.preOrder[rootIndex].precision;
        for (int i = 0; i < dim; ++i) {
            final double deltaI = workspace.mean[i] - workspace.mean2[i];
            for (int j = 0; j < dim; ++j) {
                workspace.covarianceAdjoint[i][j] =
                        workspace.covariance[i][j] + deltaI * (workspace.mean[j] - workspace.mean2[j]);
            }
        }

        GaussianMatrixOps.multiplyFlatByMatrix(
                priorPrecision, workspace.covarianceAdjoint, workspace.transitionMatrix, dim);
        GaussianMatrixOps.multiplyMatrixByFlat(
                workspace.transitionMatrix, priorPrecision, workspace.covarianceAdjoint, dim);
        for (int i = 0; i < dim; ++i) {
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                gradQ[iOff + j] += stateStore.lastRootDiffusionScale
                        * (-0.5 * priorPrecision[iOff + j] + 0.5 * workspace.covarianceAdjoint[i][j]);
            }
        }
    }

    private static void clearSquare(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            Arrays.fill(matrix[i], 0.0);
        }
    }
}
