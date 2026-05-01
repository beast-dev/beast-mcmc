package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedTransitionCapability;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockBranchGradientWorkspace;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.TipObservationModelWorkspace;

import java.util.EnumSet;

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
    final TipObservationModelWorkspace observationWorkspace;
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
    private final EnumSet<WorkspaceCapability> capabilities;

    BranchGradientWorkspace(final int dim) {
        this(dim, EnumSet.allOf(WorkspaceCapability.class));
    }

    BranchGradientWorkspace(final int dim, final EnumSet<WorkspaceCapability> capabilities) {
        this.capabilities = EnumSet.copyOf(capabilities);
        this.traversal = capabilities.contains(WorkspaceCapability.TRAVERSAL)
                ? new TraversalWorkspace(dim)
                : null;
        this.adjoint = capabilities.contains(WorkspaceCapability.ADJOINTS)
                ? new BranchAdjointWorkspace(dim)
                : null;
        this.gradient = hasGradientCapability(capabilities)
                ? new GradientPullbackWorkspace(dim, capabilities)
                : null;

        this.traitCovariance = traversal == null ? null : traversal.traitCovariance;
        this.transition = traversal == null ? null : traversal.transition;
        this.state = traversal == null ? null : traversal.state;
        this.siblingProduct = traversal == null ? null : traversal.siblingProduct;
        this.downwardParentState = traversal == null ? null : traversal.downwardParentState;
        this.gaussianWorkspace = traversal == null ? null : traversal.gaussianWorkspace;
        this.observationWorkspace = traversal == null ? null : traversal.observationWorkspace;

        this.combinedState = adjoint == null ? null : adjoint.combinedState;
        this.parentPosterior = adjoint == null ? null : adjoint.parentPosterior;
        this.pairState = adjoint == null ? null : adjoint.pairState;
        this.contributionWorkspace = adjoint == null ? null : adjoint.contributionWorkspace;
        this.transitionAdjointWorkspace = adjoint == null ? null : adjoint.transitionAdjointWorkspace;
        this.contribution = adjoint == null ? null : adjoint.contribution;
        this.adjoints = adjoint == null ? null : adjoint.adjoints;
        this.observedIndexScratch = adjoint == null ? null : adjoint.observedIndexScratch;
        this.missingIndexScratch = adjoint == null ? null : adjoint.missingIndexScratch;
        this.reducedIndexByTraitScratch = adjoint == null ? null : adjoint.reducedIndexByTraitScratch;
        this.mean = adjoint == null ? null : adjoint.mean;
        this.mean2 = adjoint == null ? null : adjoint.mean2;
        this.covariance = adjoint == null ? null : adjoint.covariance;
        this.varianceFlat = adjoint == null ? null : adjoint.varianceFlat;
        this.precisionFlat = adjoint == null ? null : adjoint.precisionFlat;
        this.reducedPrecisionFlatScratch = adjoint == null ? null : adjoint.reducedPrecisionFlatScratch;
        this.reducedCovarianceFlatScratch = adjoint == null ? null : adjoint.reducedCovarianceFlatScratch;
        this.reducedInformationScratch = adjoint == null ? null : adjoint.reducedInformationScratch;
        this.reducedMeanScratch = adjoint == null ? null : adjoint.reducedMeanScratch;

        this.orthogonalStationaryMeanScratch =
                gradient == null ? null : gradient.orthogonalStationaryMeanScratch;
        this.orthogonalCompressedGradientScratch =
                gradient == null ? null : gradient.orthogonalCompressedGradientScratch;
        this.orthogonalNativeGradientScratch =
                gradient == null ? null : gradient.orthogonalNativeGradientScratch;
        this.orthogonalRotationGradientFlatScratch =
                gradient == null ? null : gradient.orthogonalRotationGradientFlatScratch;
        this.orthogonalRotationGradientScratch =
                gradient == null ? null : gradient.orthogonalRotationGradientScratch;
        this.covariance2 = gradient == null ? null : gradient.covariance2;
        this.transitionMatrixFlat = gradient == null ? null : gradient.transitionMatrixFlat;
        this.covarianceAdjointFlat = gradient == null ? null : gradient.covarianceAdjointFlat;
        this.matrixProductFlat = gradient == null ? null : gradient.matrixProductFlat;
        this.localGradientA = gradient == null ? null : gradient.localGradientA;
        this.localGradientQ = gradient == null ? null : gradient.localGradientQ;
        this.localGradientMuVector = gradient == null ? null : gradient.localGradientMuVector;
        this.localGradientMuScalar = gradient == null ? null : gradient.localGradientMuScalar;

        validate(dim);
    }

    void validate(final int dim) {
        if (hasCapability(WorkspaceCapability.TRAVERSAL)) {
            requireLength(traitCovariance, dim, dim, "traitCovariance");
        }
        if (hasCapability(WorkspaceCapability.ADJOINTS)) {
            requireLength(covariance, dim, dim, "covariance");
            requireLength(varianceFlat, dim * dim, "varianceFlat");
            requireLength(precisionFlat, dim * dim, "precisionFlat");
            requireLength(observedIndexScratch, dim, "observedIndexScratch");
            requireLength(missingIndexScratch, dim, "missingIndexScratch");
            requireLength(reducedIndexByTraitScratch, dim, "reducedIndexByTraitScratch");
        }
        if (hasCapability(WorkspaceCapability.DENSE_GRADIENT)) {
            requireLength(covariance2, dim, dim, "covariance2");
            requireLength(transitionMatrixFlat, dim * dim, "transitionMatrixFlat");
            requireLength(covarianceAdjointFlat, dim * dim, "covarianceAdjointFlat");
            requireLength(matrixProductFlat, dim * dim, "matrixProductFlat");
            requireLength(localGradientA, dim * dim, "localGradientA");
            requireLength(localGradientQ, dim * dim, "localGradientQ");
            requireLength(localGradientMuVector, dim, "localGradientMuVector");
            requireLength(localGradientMuScalar, 1, "localGradientMuScalar");
        }
        if (hasCapability(WorkspaceCapability.ORTHOGONAL_BLOCK_GRADIENT)) {
            requireLength(orthogonalCompressedGradientScratch, dim + 2 * (dim / 2),
                    "orthogonalCompressedGradientScratch");
            requireLength(orthogonalNativeGradientScratch, ((dim & 1) == 1 ? 1 : 0) + 3 * (dim / 2),
                    "orthogonalNativeGradientScratch");
            requireLength(orthogonalRotationGradientFlatScratch, dim * dim,
                    "orthogonalRotationGradientFlatScratch");
        }
    }

    private static void requireLength(final double[] array, final int expected, final String name) {
        if (array == null || array.length < expected) {
            throw new IllegalStateException(name + " length must be at least " + expected);
        }
    }

    private static void requireLength(final int[] array, final int expected, final String name) {
        if (array == null || array.length < expected) {
            throw new IllegalStateException(name + " length must be at least " + expected);
        }
    }

    private static void requireLength(final double[][] matrix,
                                      final int rows,
                                      final int cols,
                                      final String name) {
        if (matrix == null || matrix.length < rows) {
            throw new IllegalStateException(name + " row count must be at least " + rows);
        }
        for (int i = 0; i < rows; ++i) {
            if (matrix[i] == null || matrix[i].length < cols) {
                throw new IllegalStateException(name + " column count must be at least " + cols);
            }
        }
    }

    double[] localGradientMu(final int gradientLength, final int dim) {
        return denseGradient().localGradientMu(gradientLength, dim);
    }

    void clearLocalGradientBuffers(final int gradALength,
                                   final int gradMuLength,
                                   final int dim,
                                   final boolean orthogonalSelection,
                                   final int compressedGradientLength) {
        gradient.clearLocalGradientBuffers(
                gradALength, gradMuLength, dim, orthogonalSelection, compressedGradientLength);
    }

    DenseGradientWorkspace denseGradient() {
        if (gradient == null || gradient.dense == null) {
            throw new IllegalStateException("Dense gradient workspace capability was not requested");
        }
        return gradient.dense;
    }

    OrthogonalBlockGradientWorkspace orthogonalGradient() {
        if (gradient == null || gradient.orthogonal == null) {
            throw new IllegalStateException("Orthogonal-block gradient workspace capability was not requested");
        }
        return gradient.orthogonal;
    }

    OrthogonalBlockBranchGradientWorkspace
    ensureOrthogonalBranchWorkspace(final OrthogonalBlockCanonicalParameterization orthogonalSelection) {
        return orthogonalGradient().ensureOrthogonalBranchWorkspace(orthogonalSelection);
    }

    CanonicalBranchWorkspace
    ensureSpecializedBranchWorkspace(final CanonicalPreparedTransitionCapability selection) {
        return orthogonalGradient().ensureSpecializedBranchWorkspace(selection);
    }

    boolean hasCapability(final WorkspaceCapability capability) {
        return capabilities.contains(capability);
    }

    private static boolean hasGradientCapability(final EnumSet<WorkspaceCapability> capabilities) {
        return capabilities.contains(WorkspaceCapability.DENSE_GRADIENT)
                || capabilities.contains(WorkspaceCapability.ORTHOGONAL_BLOCK_GRADIENT);
    }
}
