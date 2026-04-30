package dr.evomodel.treedatalikelihood.continuous.integration;

import dr.evomodel.continuous.ou.CanonicalOUKernel;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.GaussianMatrixOps;

final class DenseCanonicalSelectionGradientPullback implements CanonicalSelectionGradientPullback {

    private final int dimension;

    DenseCanonicalSelectionGradientPullback(final int dimension) {
        this.dimension = dimension;
    }

    @Override
    public void initialize(final BranchGradientWorkspace workspace,
                           final double[] gradA,
                           final double[] gradMu) { }

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
        final GradientPullbackWorkspace gradient = workspace.gradient;
        final CanonicalOUKernel kernel = processModel.getCanonicalKernel();
        final double branchLength = inputs.getBranchLength(activeIndex);

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

        kernel.fillTransitionMatrixFlat(branchLength, gradient.transitionMatrixFlat);
        accumulateStationaryMeanGradientFlat(
                gradient.transitionMatrixFlat,
                localAdjoints.dLogL_df,
                gradMu);
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
        GaussianMatrixOps.transposeFlatSquareInPlace(gradA, dimension);
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
