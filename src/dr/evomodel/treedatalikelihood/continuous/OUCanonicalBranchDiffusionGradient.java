package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.canonical.CanonicalNativeBranchGradientCapability;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;

final class OUCanonicalBranchDiffusionGradient {

    private final OUGaussianBranchTransitionProvider branchTransitionProvider;
    private final OUCanonicalBranchWiring branchWiring;
    private final OUProcessModel processModel;
    private final int dimension;
    private final OUCanonicalBranchGradientScratch scratch;

    OUCanonicalBranchDiffusionGradient(final OUGaussianBranchTransitionProvider branchTransitionProvider,
                                       final OUCanonicalBranchWiring branchWiring,
                                       final OUProcessModel processModel,
                                       final int dimension,
                                       final OUCanonicalBranchGradientScratch scratch) {
        this.branchTransitionProvider = branchTransitionProvider;
        this.branchWiring = branchWiring;
        this.processModel = processModel;
        this.dimension = dimension;
        this.scratch = scratch;
    }

    double[] getGradientWrtVariance(final double branchLength,
                                    final double[] optimum,
                                    final BranchSufficientStatistics statistics) {
        refreshProcessSnapshots();
        final CanonicalLocalTransitionAdjoints localAdjoints = scratch.localAdjoints;
        branchWiring.fillLocalAdjoints(branchLength, optimum, statistics, localAdjoints);
        OUCanonicalBranchGradientUtils.zero(scratch.matrixGradient);
        if (processModel.getSelectionMatrixParameterization()
                instanceof CanonicalNativeBranchGradientCapability) {
            OUCanonicalBranchGradientUtils.transposeFromFlatInto(
                    localAdjoints.dLogL_dOmega, scratch.covarianceAdjoint, dimension);
            ((CanonicalNativeBranchGradientCapability)
                    processModel.getSelectionMatrixParameterization())
                    .accumulateDiffusionGradient(
                            processModel.getDiffusionMatrix(),
                            branchLength,
                            scratch.covarianceAdjoint,
                            scratch.matrixGradient);
        } else {
            OUCanonicalBranchGradientUtils.copyFromFlatInto(
                    localAdjoints.dLogL_dOmega, scratch.covarianceAdjoint, dimension);
            processModel.accumulateDiffusionGradient(branchLength, scratch.covarianceAdjoint, scratch.matrixGradient);
        }
        return scratch.matrixGradient.clone();
    }

    private void refreshProcessSnapshots() {
        branchTransitionProvider.getProcessModel();
    }
}
