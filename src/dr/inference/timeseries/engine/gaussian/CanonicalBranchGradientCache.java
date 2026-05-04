package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Per-smoother cache of branch-local canonical adjoints shared by gradient formulas.
 */
public final class CanonicalBranchGradientCache {

    private final CanonicalBranchMessageContribution[] contributions;
    private final CanonicalLocalTransitionAdjoints[] adjoints;
    private final double[][] transitionMatrices;
    private final GaussianTransitionRepresentation transitionRepresentation;
    private final TimeGrid timeGrid;
    private final double[][] transitionMatrixWorkspace;
    private final double[][] transitionCovarianceWorkspace;
    private final double[] transitionCovarianceFlatWorkspace;
    private final double[] transitionOffsetWorkspace;
    private final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
    private final CanonicalTransitionAdjointUtils.Workspace adjointWorkspace;
    private SharedCanonicalTimeSeriesSchedule sharedSchedule;
    private boolean known;
    private long buildCount;

    CanonicalBranchGradientCache(final int timeCount,
                                 final int stateDimension,
                                 final GaussianTransitionRepresentation transitionRepresentation,
                                 final TimeGrid timeGrid) {
        final int branchCount = Math.max(0, timeCount - 1);
        this.contributions = new CanonicalBranchMessageContribution[branchCount];
        this.adjoints = new CanonicalLocalTransitionAdjoints[branchCount];
        this.transitionMatrices = new double[branchCount][stateDimension * stateDimension];
        this.transitionRepresentation = transitionRepresentation;
        this.timeGrid = timeGrid;
        this.transitionMatrixWorkspace = new double[stateDimension][stateDimension];
        this.transitionCovarianceWorkspace = new double[stateDimension][stateDimension];
        this.transitionCovarianceFlatWorkspace = new double[stateDimension * stateDimension];
        this.transitionOffsetWorkspace = new double[stateDimension];
        for (int i = 0; i < branchCount; ++i) {
            contributions[i] = new CanonicalBranchMessageContribution(stateDimension);
            adjoints[i] = new CanonicalLocalTransitionAdjoints(stateDimension);
        }
        this.contributionWorkspace = new CanonicalBranchMessageContributionUtils.Workspace(stateDimension);
        this.adjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(stateDimension);
        this.sharedSchedule = null;
        this.known = false;
        this.buildCount = 0L;
    }

    public void setSharedSchedule(final SharedCanonicalTimeSeriesSchedule sharedSchedule) {
        this.sharedSchedule = sharedSchedule;
        makeDirty();
    }

    public void makeDirty() {
        known = false;
    }

    public void ensure(final CanonicalForwardTrajectory trajectory) {
        if (known) {
            return;
        }
        ++buildCount;
        for (int t = 0; t < contributions.length; ++t) {
            if (sharedSchedule == null
                    || !sharedSchedule.fillContribution(
                    t,
                    trajectory.branchPairStates[t],
                    contributionWorkspace,
                    contributions[t])) {
                CanonicalBranchMessageContributionUtils.fillFromPairState(
                        trajectory.branchPairStates[t],
                        contributionWorkspace,
                        contributions[t]);
            }
            if (transitionRepresentation == null || timeGrid == null) {
                CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                        trajectory.transitions[t],
                        contributions[t],
                        adjointWorkspace,
                        adjoints[t]);
            } else {
                transitionRepresentation.getTransitionMatrix(t, t + 1, timeGrid, transitionMatrixWorkspace);
                transitionRepresentation.getTransitionOffset(t, t + 1, timeGrid, transitionOffsetWorkspace);
                transitionRepresentation.getTransitionCovariance(t, t + 1, timeGrid, transitionCovarianceWorkspace);
                GaussianMatrixOps.copyMatrixToFlat(
                        transitionMatrixWorkspace,
                        adjointWorkspace.transitionMatrix,
                        adjointWorkspace.getDimension());
                GaussianMatrixOps.copyMatrixToFlat(
                        transitionCovarianceWorkspace,
                        transitionCovarianceFlatWorkspace,
                        adjointWorkspace.getDimension());
                CanonicalTransitionAdjointUtils.fillFromMoments(
                        trajectory.transitions[t].precisionYY,
                        transitionCovarianceFlatWorkspace,
                        adjointWorkspace.transitionMatrix,
                        transitionOffsetWorkspace,
                        contributions[t],
                        adjointWorkspace,
                        adjoints[t]);
            }
            System.arraycopy(adjointWorkspace.transitionMatrix, 0,
                    transitionMatrices[t], 0, transitionMatrices[t].length);
        }
        known = true;
    }

    public CanonicalBranchMessageContribution getContribution(final int branchIndex) {
        return contributions[branchIndex];
    }

    public CanonicalLocalTransitionAdjoints getAdjoints(final int branchIndex) {
        return adjoints[branchIndex];
    }

    public double[] getTransitionMatrix(final int branchIndex) {
        return transitionMatrices[branchIndex];
    }

    public long getBuildCount() {
        return buildCount;
    }
}
