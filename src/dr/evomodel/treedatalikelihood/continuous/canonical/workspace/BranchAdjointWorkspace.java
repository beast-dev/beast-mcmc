package dr.evomodel.treedatalikelihood.continuous.canonical.workspace;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.PartialIdentityTipObservationWorkspace;

public final class BranchAdjointWorkspace {
    public final CanonicalGaussianState combinedState;
    public final CanonicalGaussianState parentPosterior;
    public final CanonicalGaussianState pairState;
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
    public final double[] covarianceFlat;
    public final double[] varianceFlat;
    public final double[] precisionFlat;
    public final double[] reducedPrecisionFlatScratch;
    public final double[] reducedCovarianceFlatScratch;
    public final double[] reducedInformationScratch;
    public final double[] reducedMeanScratch;
    public final GaussianFormConverter.Workspace converterWorkspace;

    public BranchAdjointWorkspace(final int dim) {
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

    public void fillMomentsFromCanonical(final CanonicalGaussianState state,
                                  final double[] meanOut,
                                  final double[] covarianceOut,
                                  final int dim) {
        GaussianFormConverter.fillMomentsFromState(
                state,
                meanOut,
                covarianceOut,
                dim,
                converterWorkspace);
    }
}
