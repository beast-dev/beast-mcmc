package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.treedatalikelihood.continuous.gaussian.message.GaussianMatrixOps;

import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalTransitionAdjointUtils;

import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.gaussian.DiffusionMatrixParameterization;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Diffusion gradient formula driven by branch-local canonical posterior states.
 */
public final class CanonicalDiffusionMatrixGradientFormula implements CanonicalGradientFormula {

    private final DiffusionMatrixParameterization diffusionParameterization;
    private final int stateDimension;
    private final OUProcessModel processModel;

    private final CanonicalBranchPosterior branchPosterior;
    private final CanonicalBranchMessageContribution localContribution;
    private final CanonicalLocalTransitionAdjoints localAdjoints;
    private final CanonicalTransitionAdjointUtils.Workspace canonicalAdjointWorkspace;
    private final double[][] covarianceAdjointScratch;

    public CanonicalDiffusionMatrixGradientFormula(final DiffusionMatrixParameterization diffusionParameterization,
                                                   final int stateDimension) {
        this(null, diffusionParameterization, stateDimension);
    }

    public CanonicalDiffusionMatrixGradientFormula(final OUProcessModel processModel,
                                                   final DiffusionMatrixParameterization diffusionParameterization,
                                                   final int stateDimension) {
        if (diffusionParameterization == null) {
            throw new IllegalArgumentException("diffusionParameterization must not be null");
        }
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        this.processModel = processModel;
        this.diffusionParameterization = diffusionParameterization;
        this.stateDimension = stateDimension;

        this.branchPosterior = new CanonicalBranchPosterior(stateDimension);
        this.localContribution = new CanonicalBranchMessageContribution(stateDimension);
        this.localAdjoints = new CanonicalLocalTransitionAdjoints(stateDimension);
        this.canonicalAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(stateDimension);
        this.covarianceAdjointScratch = new double[stateDimension][stateDimension];
    }

    @Override
    public boolean supportsParameter(final Parameter parameter) {
        return diffusionParameterization.supportsParameter(parameter);
    }

    @Override
    public double[] computeGradient(final Parameter parameter,
                                    final CanonicalForwardTrajectory trajectory,
                                    final GaussianTransitionRepresentation repr,
                                    final TimeGrid timeGrid) {
        final int d = stateDimension;
        final int T = trajectory.timeCount;
        final double[] gradientAccumulator = new double[d * d];
        final OrthogonalBlockCanonicalParameterization orthogonalParameterization =
                getOrthogonalParameterizationIfAvailable();

        for (int t = 0; t < T - 1; ++t) {
            branchPosterior.fillFromCanonicalPairState(trajectory.branchPairStates[t]);
            branchPosterior.fillLocalMessageContribution(localContribution);
            CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                    trajectory.transitions[t],
                    localContribution,
                    canonicalAdjointWorkspace,
                    localAdjoints);
            GaussianMatrixOps.copyFlatToMatrix(localAdjoints.dLogL_dOmega, covarianceAdjointScratch, d);
            if (orthogonalParameterization != null) {
                orthogonalParameterization.accumulateDiffusionGradient(
                        processModel.getDiffusionMatrix(),
                        timeGrid.getDelta(t, t + 1),
                        covarianceAdjointScratch,
                        gradientAccumulator);
            } else {
                repr.accumulateDiffusionGradient(
                        t, t + 1, timeGrid, covarianceAdjointScratch, gradientAccumulator);
            }
        }

        if (parameter == diffusionParameterization.getMatrixParameter()) {
            return gradientAccumulator;
        }
        return diffusionParameterization.pullBackGradient(parameter, gradientAccumulator);
    }

    private OrthogonalBlockCanonicalParameterization getOrthogonalParameterizationIfAvailable() {
        if (processModel == null) {
            return null;
        }
        if (processModel.getSelectionMatrixParameterization()
                instanceof OrthogonalBlockCanonicalParameterization) {
            return (OrthogonalBlockCanonicalParameterization)
                    processModel.getSelectionMatrixParameterization();
        }
        return null;
    }
}
