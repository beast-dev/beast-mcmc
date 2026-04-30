package dr.evomodel.treedatalikelihood.continuous.integration;

import dr.evomodel.continuous.ou.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.CanonicalPreparedTransitionCapability;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockBranchGradientWorkspace;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalTransitionAdjointUtils;

final class BranchGradientWorkspace {
    final TraversalWorkspace traversal;
    final BranchAdjointWorkspace adjoint;
    final GradientPullbackWorkspace gradient;

    final double[][] traitCovariance;
    final CanonicalGaussianTransition transition;
    final CanonicalGaussianState state;
    final CanonicalGaussianState siblingProduct;
    final CanonicalGaussianState downwardParentState;
    final CanonicalGaussianState combinedState;
    final CanonicalGaussianState parentPosterior;
    final CanonicalGaussianState pairState;
    final CanonicalGaussianMessageOps.Workspace gaussianWorkspace;
    final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
    final CanonicalTransitionAdjointUtils.Workspace transitionAdjointWorkspace;
    final CanonicalBranchMessageContribution contribution;
    final CanonicalLocalTransitionAdjoints adjoints;
    final int[] observedIndexScratch;
    final int[] missingIndexScratch;
    final int[] reducedIndexByTraitScratch;
    final double[] mean;
    final double[] mean2;
    final double[] orthogonalStationaryMeanScratch;
    final double[] orthogonalCompressedGradientScratch;
    final double[] orthogonalNativeGradientScratch;
    final double[] orthogonalRotationGradientFlatScratch;
    final double[][] orthogonalRotationGradientScratch;
    final double[][] covariance;
    final double[][] covariance2;
    final double[] transitionMatrixFlat;
    final double[] covarianceAdjointFlat;
    final double[] matrixProductFlat;
    final double[] varianceFlat;
    final double[] precisionFlat;
    final double[] reducedPrecisionFlatScratch;
    final double[] reducedCovarianceFlatScratch;
    final double[] reducedInformationScratch;
    final double[] reducedMeanScratch;
    final double[] localGradientA;
    final double[] localGradientQ;
    final double[] localGradientMuVector;
    final double[] localGradientMuScalar;

    BranchGradientWorkspace(final int dim) {
        this.traversal = new TraversalWorkspace(dim);
        this.adjoint = new BranchAdjointWorkspace(dim);
        this.gradient = new GradientPullbackWorkspace(dim);

        this.traitCovariance = traversal.traitCovariance;
        this.transition = traversal.transition;
        this.state = traversal.state;
        this.siblingProduct = traversal.siblingProduct;
        this.downwardParentState = traversal.downwardParentState;
        this.gaussianWorkspace = traversal.gaussianWorkspace;

        this.combinedState = adjoint.combinedState;
        this.parentPosterior = adjoint.parentPosterior;
        this.pairState = adjoint.pairState;
        this.contributionWorkspace = adjoint.contributionWorkspace;
        this.transitionAdjointWorkspace = adjoint.transitionAdjointWorkspace;
        this.contribution = adjoint.contribution;
        this.adjoints = adjoint.adjoints;
        this.observedIndexScratch = adjoint.observedIndexScratch;
        this.missingIndexScratch = adjoint.missingIndexScratch;
        this.reducedIndexByTraitScratch = adjoint.reducedIndexByTraitScratch;
        this.mean = adjoint.mean;
        this.mean2 = adjoint.mean2;
        this.covariance = adjoint.covariance;
        this.varianceFlat = adjoint.varianceFlat;
        this.precisionFlat = adjoint.precisionFlat;
        this.reducedPrecisionFlatScratch = adjoint.reducedPrecisionFlatScratch;
        this.reducedCovarianceFlatScratch = adjoint.reducedCovarianceFlatScratch;
        this.reducedInformationScratch = adjoint.reducedInformationScratch;
        this.reducedMeanScratch = adjoint.reducedMeanScratch;

        this.orthogonalStationaryMeanScratch = gradient.orthogonalStationaryMeanScratch;
        this.orthogonalCompressedGradientScratch = gradient.orthogonalCompressedGradientScratch;
        this.orthogonalNativeGradientScratch = gradient.orthogonalNativeGradientScratch;
        this.orthogonalRotationGradientFlatScratch = gradient.orthogonalRotationGradientFlatScratch;
        this.orthogonalRotationGradientScratch = gradient.orthogonalRotationGradientScratch;
        this.covariance2 = gradient.covariance2;
        this.transitionMatrixFlat = gradient.transitionMatrixFlat;
        this.covarianceAdjointFlat = gradient.covarianceAdjointFlat;
        this.matrixProductFlat = gradient.matrixProductFlat;
        this.localGradientA = gradient.localGradientA;
        this.localGradientQ = gradient.localGradientQ;
        this.localGradientMuVector = gradient.localGradientMuVector;
        this.localGradientMuScalar = gradient.localGradientMuScalar;
    }

    double[] localGradientMu(final int gradientLength, final int dim) {
        return gradient.localGradientMu(gradientLength, dim);
    }

    void clearLocalGradientBuffers(final int gradALength,
                                   final int gradMuLength,
                                   final int dim,
                                   final boolean orthogonalSelection,
                                   final int compressedGradientLength) {
        gradient.clearLocalGradientBuffers(
                gradALength, gradMuLength, dim, orthogonalSelection, compressedGradientLength);
    }

    OrthogonalBlockBranchGradientWorkspace
    ensureOrthogonalBranchWorkspace(final OrthogonalBlockCanonicalParameterization orthogonalSelection) {
        return gradient.ensureOrthogonalBranchWorkspace(orthogonalSelection);
    }

    CanonicalBranchWorkspace
    ensureSpecializedBranchWorkspace(final CanonicalPreparedTransitionCapability selection) {
        return gradient.ensureSpecializedBranchWorkspace(selection);
    }
}
