package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.PartialIdentityTipObservationWorkspace;

final class BranchAdjointWorkspace {
    final CanonicalGaussianState combinedState;
    final CanonicalGaussianState parentPosterior;
    final CanonicalGaussianState pairState;
    final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
    final CanonicalTransitionAdjointUtils.Workspace transitionAdjointWorkspace;
    final CanonicalBranchMessageContribution contribution;
    final CanonicalLocalTransitionAdjoints adjoints;
    final PartialIdentityTipObservationWorkspace partialObservationWorkspace;
    final int[] observedIndexScratch;
    final int[] missingIndexScratch;
    final int[] reducedIndexByTraitScratch;
    final double[] mean;
    final double[] mean2;
    final double[][] covariance;
    final double[] covarianceFlat;        // was double[][] covariance; row-major, length dim*dim
    final double[] varianceFlat;
    final double[] precisionFlat;
    final double[] reducedPrecisionFlatScratch;
    final double[] reducedCovarianceFlatScratch;
    final double[] reducedInformationScratch;
    final double[] reducedMeanScratch;
    final GaussianFormConverter.Workspace converterWorkspace;

    BranchAdjointWorkspace(final int dim) {
        if (dim < 1) throw new IllegalArgumentException("dim must be >= 1");
        this.combinedState = new CanonicalGaussianState(dim);
        this.parentPosterior = new CanonicalGaussianState(dim);
        this.pairState = new CanonicalGaussianState(2 * dim);
        this.contributionWorkspace = new CanonicalBranchMessageContributionUtils.Workspace(dim);
        this.transitionAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(dim);
        this.contribution = new CanonicalBranchMessageContribution(dim);
        this.adjoints = new CanonicalLocalTransitionAdjoints(dim);
        this.partialObservationWorkspace = new PartialIdentityTipObservationWorkspace(dim);
        this.observedIndexScratch = partialObservationWorkspace.partition.observedIndices;
        this.missingIndexScratch = partialObservationWorkspace.partition.missingIndices;
        this.reducedIndexByTraitScratch = partialObservationWorkspace.partition.reducedIndexByTrait;
        this.mean = new double[dim];
        this.mean2 = partialObservationWorkspace.missingMean;
        this.covariance = new double[dim][dim];
        this.covarianceFlat = new double[dim * dim];
        this.varianceFlat = new double[dim * dim];
        this.precisionFlat = new double[dim * dim];
        this.reducedPrecisionFlatScratch = partialObservationWorkspace.reducedPrecision;
        this.reducedCovarianceFlatScratch = partialObservationWorkspace.reducedCovariance;
        this.reducedInformationScratch = partialObservationWorkspace.reducedInformation;
        this.reducedMeanScratch = partialObservationWorkspace.reducedMean;
        this.converterWorkspace = new GaussianFormConverter.Workspace();
        this.converterWorkspace.ensureDim(dim);
    }

    void fillMomentsFromCanonical(final CanonicalGaussianState state,
                                  final double[] meanOut,
                                  final double[][] covarianceOut,
                                  final int dim) {
        GaussianFormConverter.fillMomentsFromState(
                state,
                meanOut,
                covarianceFlat,
                dim,
                converterWorkspace);
        MatrixOps.fromFlat(covarianceFlat, covarianceOut, dim);
    }
}
