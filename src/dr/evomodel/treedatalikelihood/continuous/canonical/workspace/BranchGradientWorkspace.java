package dr.evomodel.treedatalikelihood.continuous.canonical.workspace;

import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedTransitionCapability;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.PartialIdentityTipObservationWorkspace;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.TipObservationModelWorkspace;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.TipParentMessageWorkspace;

import java.util.EnumSet;

public final class BranchGradientWorkspace {
    public final TraversalWorkspace traversal;
    public final BranchAdjointWorkspace adjoint;
    public final GradientPullbackWorkspace gradient;

    public final double[][] traitCovariance;
    public final CanonicalGaussianTransition transition;
    public final CanonicalGaussianState state;
    public final CanonicalGaussianState siblingProduct;
    public final CanonicalGaussianState downwardParentState;
    public final CanonicalGaussianState combinedState;
    public final CanonicalGaussianState parentPosterior;
    public final CanonicalGaussianState pairState;
    public final CanonicalGaussianMessageOps.Workspace gaussianWorkspace;
    public final TipObservationModelWorkspace observationWorkspace;
    public final TipParentMessageWorkspace tipParentMessageWorkspace;
    public final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
    public final CanonicalTransitionAdjointUtils.Workspace transitionAdjointWorkspace;
    public final CanonicalBranchMessageContribution contribution;
    public final CanonicalLocalTransitionAdjoints adjoints;
    public final PartialIdentityTipObservationWorkspace partialObservationWorkspace;
    public final int[] observedIndexScratch;
    public final int[] missingIndexScratch;
    public final int[] reducedIndexByTraitScratch;
    public final double[] mean;
    public final double[] mean2;
    public final double[] specializedStationaryMeanScratch;
    public final double[] specializedCompressedGradientScratch;
    public final double[] specializedNativeGradientScratch;
    public final double[] specializedRotationGradientFlatScratch;
    public final double[][] specializedRotationGradientScratch;
    public final double[][] covariance;
    public final double[][] covariance2;
    public final double[] transitionMatrixFlat;
    public final double[] covarianceAdjointFlat;
    public final double[] matrixProductFlat;
    public final double[] varianceFlat;
    public final double[] precisionFlat;
    public final double[] reducedPrecisionFlatScratch;
    public final double[] reducedCovarianceFlatScratch;
    public final double[] reducedInformationScratch;
    public final double[] reducedMeanScratch;
    public final double[] localGradientA;
    public final double[] localGradientQ;
    public final double[] localGradientMuVector;
    public final double[] localGradientMuScalar;
    private final EnumSet<WorkspaceCapability> capabilities;

    public BranchGradientWorkspace(final int dim) {
        this(dim, EnumSet.allOf(WorkspaceCapability.class));
    }

    public BranchGradientWorkspace(final int dim, final EnumSet<WorkspaceCapability> capabilities) {
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
        this.observationWorkspace      = traversal == null ? null : traversal.observationWorkspace;
        this.tipParentMessageWorkspace = traversal == null ? null : traversal.tipParentMessageWorkspace;

        this.combinedState = adjoint == null ? null : adjoint.combinedState;
        this.parentPosterior = adjoint == null ? null : adjoint.parentPosterior;
        this.pairState = adjoint == null ? null : adjoint.pairState;
        this.contributionWorkspace = adjoint == null ? null : adjoint.contributionWorkspace;
        this.transitionAdjointWorkspace = adjoint == null ? null : adjoint.transitionAdjointWorkspace;
        this.contribution = adjoint == null ? null : adjoint.contribution;
        this.adjoints = adjoint == null ? null : adjoint.adjoints;
        this.partialObservationWorkspace = adjoint == null ? null : adjoint.partialObservationWorkspace;
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

        this.specializedStationaryMeanScratch =
                gradient == null ? null : gradient.specializedStationaryMeanScratch;
        this.specializedCompressedGradientScratch =
                gradient == null ? null : gradient.specializedCompressedGradientScratch;
        this.specializedNativeGradientScratch =
                gradient == null ? null : gradient.specializedNativeGradientScratch;
        this.specializedRotationGradientFlatScratch =
                gradient == null ? null : gradient.specializedRotationGradientFlatScratch;
        this.specializedRotationGradientScratch =
                gradient == null ? null : gradient.specializedRotationGradientScratch;
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

    public void validate(final int dim) {
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
        if (hasCapability(WorkspaceCapability.SPECIALIZED_GRADIENT)) {
            requireLength(specializedCompressedGradientScratch, dim + 2 * (dim / 2),
                    "specializedCompressedGradientScratch");
            requireLength(specializedNativeGradientScratch, ((dim & 1) == 1 ? 1 : 0) + 3 * (dim / 2),
                    "specializedNativeGradientScratch");
            requireLength(specializedRotationGradientFlatScratch, dim * dim,
                    "specializedRotationGradientFlatScratch");
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

    public double[] localGradientMu(final int gradientLength, final int dim) {
        return denseGradient().localGradientMu(gradientLength, dim);
    }

    public void clearLocalGradientBuffers(final int gradALength,
                                   final int gradMuLength,
                                   final int dim,
                                   final boolean orthogonalSelection,
                                   final int compressedGradientLength) {
        gradient.clearLocalGradientBuffers(
                gradALength, gradMuLength, dim, orthogonalSelection, compressedGradientLength);
    }

    public DenseGradientWorkspace denseGradient() {
        if (gradient == null || gradient.dense == null) {
            throw new IllegalStateException("Dense gradient workspace capability was not requested");
        }
        return gradient.dense;
    }

    public SpecializedGradientWorkspace specializedGradient() {
        if (gradient == null || gradient.specialized == null) {
            throw new IllegalStateException("Specialized gradient workspace capability was not requested");
        }
        return gradient.specialized;
    }

    public CanonicalBranchWorkspace
    ensureSpecializedBranchWorkspace(final CanonicalPreparedTransitionCapability selection) {
        return specializedGradient().ensureSpecializedBranchWorkspace(selection);
    }

    public boolean hasCapability(final WorkspaceCapability capability) {
        return capabilities.contains(capability);
    }

    private static boolean hasGradientCapability(final EnumSet<WorkspaceCapability> capabilities) {
        return capabilities.contains(WorkspaceCapability.DENSE_GRADIENT)
                || capabilities.contains(WorkspaceCapability.SPECIALIZED_GRADIENT);
    }
}
