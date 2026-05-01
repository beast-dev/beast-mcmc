package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

import dr.evomodel.continuous.ou.canonical.CanonicalOUKernel;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedTransitionCapability;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.SelectionMatrixParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.DenseGradientWorkspace;

final class DenseCanonicalSelectionGradientPullback implements CanonicalSelectionGradientPullback {

    private final int dimension;

    DenseCanonicalSelectionGradientPullback(final int dimension) {
        this.dimension = dimension;
    }

    @Override
    public void initialize(final BranchGradientWorkspace workspace,
                           final double[] gradA,
                           final double[] gradMu) {
        workspace.gradient.invalidateTransitionMatrixCache();
    }

    @Override
    public void clearWorkerBuffers(final BranchGradientWorkspace workspace,
                                   final int gradALength,
                                   final int gradMuLength) {
        workspace.clearLocalGradientBuffers(gradALength, gradMuLength, dimension, false, 0);
    }

    @Override
    public void prepareWorkspace(final BranchGradientWorkspace workspace) { }

    @Override
    public void accumulateForBranch(final OUProcessModel processModel,
                                    final BranchGradientInputs inputs,
                                    final int activeIndex,
                                    final CanonicalLocalTransitionAdjoints localAdjoints,
                                    final BranchGradientWorkspace workspace,
                                    final double[] gradA,
                                    final double[] gradQ,
                                    final double[] gradMu) {
        final DenseGradientWorkspace gradient = workspace.denseGradient();
        final CanonicalOUKernel kernel = processModel.getCanonicalKernel();
        final double branchLength = inputs.getBranchLength(activeIndex);
        final SelectionMatrixParameterization parameterization =
                processModel.getSelectionMatrixParameterization();
        final CanonicalPreparedTransitionCapability preparedTransition =
                parameterization instanceof CanonicalPreparedTransitionCapability
                        ? (CanonicalPreparedTransitionCapability) parameterization
                        : null;

        kernel.accumulateSelectionGradientFlat(
                branchLength,
                localAdjoints.dLogL_dF,
                localAdjoints.dLogL_df,
                gradA);

        kernel.accumulateSelectionGradientFromCovarianceFlat(
                branchLength,
                localAdjoints.dLogL_dOmega,
                true,
                gradA);

        kernel.accumulateDiffusionGradientFlat(
                branchLength,
                localAdjoints.dLogL_dOmega,
                false,
                gradQ);

        fillTransitionMatrixFlat(
                kernel,
                preparedTransition,
                inputs,
                activeIndex,
                branchLength,
                workspace,
                gradient);
        accumulateStationaryMeanGradientFlat(
                gradient.transitionMatrixFlat,
                localAdjoints.dLogL_df,
                gradMu);
    }

    private void fillTransitionMatrixFlat(final CanonicalOUKernel kernel,
                                          final CanonicalPreparedTransitionCapability preparedTransition,
                                          final BranchGradientInputs inputs,
                                          final int activeIndex,
                                          final double branchLength,
                                          final BranchGradientWorkspace workspace,
                                          final DenseGradientWorkspace gradient) {
        if (gradient.hasTransitionMatrix(branchLength)) {
            return;
        }
        if (preparedTransition != null && inputs.hasPreparedBranchHandles()) {
            preparedTransition.fillTransitionMatrixPreparedFlat(
                    inputs.getPreparedBranchHandle(activeIndex),
                    workspace.ensureSpecializedBranchWorkspace(preparedTransition),
                    gradient.transitionMatrixFlat);
        } else {
            kernel.fillTransitionMatrixFlat(branchLength, gradient.transitionMatrixFlat);
        }
        gradient.cacheTransitionMatrix(branchLength);
    }

    @Override
    public void reduceWorker(final BranchGradientWorkspace worker,
                             final BranchGradientWorkspace reductionWorkspace,
                             final double[] gradA) {
        accumulateVectorInPlace(gradA, worker.localGradientA, gradA.length);
    }

    @Override
    public void finish(final BranchGradientInputs inputs,
                       final BranchGradientWorkspace workspace,
                       final double[] gradA) {
        MatrixOps.transposeInPlace(gradA, dimension);
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

    private static void accumulateVectorInPlace(final double[] target,
                                                final double[] source,
                                                final int length) {
        for (int i = 0; i < length; ++i) {
            target[i] += source[i];
        }
    }
}
