package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.representation.CachedGaussianTransitionRepresentation;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import dr.inference.timeseries.representation.TransitionMomentsView;

import java.util.Arrays;

/**
 * Per-smoother cache of branch-local canonical adjoints shared by gradient formulas.
 */
public final class CanonicalBranchGradientCache {

    private final CanonicalBranchMessageContribution[] contributions;
    private final boolean[] contributionKnown;
    private final CanonicalLocalTransitionAdjoints[] adjoints;
    private final double[][] transitionMatrices;
    private final double[][] ownedTransitionMatrices;
    private final GaussianTransitionRepresentation transitionRepresentation;
    private final CachedGaussianTransitionRepresentation cachedTransitionRepresentation;
    private final TimeGrid timeGrid;
    private final double[][] transitionMatrixWorkspace;
    private final double[][] transitionCovarianceWorkspace;
    private final double[] transitionCovarianceFlatWorkspace;
    private final double[] transitionOffsetWorkspace;
    private final TransitionMomentsView transitionMomentsView;
    private final CanonicalBranchMessageContribution workingContribution;
    private final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
    private final CanonicalTransitionAdjointUtils.Workspace adjointWorkspace;
    private SharedCanonicalTimeSeriesSchedule sharedSchedule;
    private CanonicalForwardTrajectory knownTrajectory;
    private boolean known;
    private long buildCount;

    CanonicalBranchGradientCache(final int timeCount,
                                 final int stateDimension,
                                 final GaussianTransitionRepresentation transitionRepresentation,
                                 final TimeGrid timeGrid) {
        final int branchCount = Math.max(0, timeCount - 1);
        this.contributions = new CanonicalBranchMessageContribution[branchCount];
        this.contributionKnown = new boolean[branchCount];
        this.adjoints = new CanonicalLocalTransitionAdjoints[branchCount];
        this.transitionMatrices = new double[branchCount][];
        this.ownedTransitionMatrices = new double[branchCount][stateDimension * stateDimension];
        this.transitionRepresentation = transitionRepresentation;
        this.cachedTransitionRepresentation =
                transitionRepresentation instanceof CachedGaussianTransitionRepresentation
                        ? (CachedGaussianTransitionRepresentation) transitionRepresentation
                        : null;
        this.timeGrid = timeGrid;
        this.transitionMatrixWorkspace = new double[stateDimension][stateDimension];
        this.transitionCovarianceWorkspace = new double[stateDimension][stateDimension];
        this.transitionCovarianceFlatWorkspace = new double[stateDimension * stateDimension];
        this.transitionOffsetWorkspace = new double[stateDimension];
        this.transitionMomentsView = new TransitionMomentsView();
        this.workingContribution = new CanonicalBranchMessageContribution(stateDimension);
        for (int i = 0; i < branchCount; ++i) {
            adjoints[i] = new CanonicalLocalTransitionAdjoints(stateDimension);
            transitionMatrices[i] = ownedTransitionMatrices[i];
        }
        this.contributionWorkspace = new CanonicalBranchMessageContributionUtils.Workspace(stateDimension);
        this.adjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(stateDimension);
        this.sharedSchedule = null;
        this.knownTrajectory = null;
        this.known = false;
        this.buildCount = 0L;
    }

    public void setSharedSchedule(final SharedCanonicalTimeSeriesSchedule sharedSchedule) {
        this.sharedSchedule = sharedSchedule;
        makeDirty();
    }

    public void makeDirty() {
        known = false;
        knownTrajectory = null;
        Arrays.fill(contributionKnown, false);
    }

    public void ensure(final CanonicalForwardTrajectory trajectory) {
        if (known) {
            return;
        }
        ++buildCount;
        knownTrajectory = trajectory;
        Arrays.fill(contributionKnown, false);
        for (int t = 0; t < adjoints.length; ++t) {
            fillContribution(t, trajectory, workingContribution);
            if (transitionRepresentation == null || timeGrid == null) {
                CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                        trajectory.transitions[t],
                        workingContribution,
                        adjointWorkspace,
                        adjoints[t]);
                System.arraycopy(adjointWorkspace.transitionMatrix, 0,
                        ownedTransitionMatrices[t], 0, ownedTransitionMatrices[t].length);
                transitionMatrices[t] = ownedTransitionMatrices[t];
            } else if (cachedTransitionRepresentation != null
                    && cachedTransitionRepresentation.getTransitionMomentsView(
                    t,
                    t + 1,
                    timeGrid,
                    transitionMomentsView)) {
                transitionMatrices[t] = transitionMomentsView.getTransitionMatrix();
                CanonicalTransitionAdjointUtils.fillFromMoments(
                        trajectory.transitions[t].precisionYY,
                        transitionMomentsView.getTransitionCovariance(),
                        transitionMomentsView.getTransitionMatrix(),
                        transitionMomentsView.getTransitionOffset(),
                        workingContribution,
                        adjointWorkspace,
                        adjoints[t]);
            } else {
                transitionRepresentation.getTransitionMatrix(t, t + 1, timeGrid, transitionMatrixWorkspace);
                transitionRepresentation.getTransitionOffset(t, t + 1, timeGrid, transitionOffsetWorkspace);
                transitionRepresentation.getTransitionCovariance(t, t + 1, timeGrid, transitionCovarianceWorkspace);
                GaussianMatrixOps.copyMatrixToFlat(
                        transitionMatrixWorkspace,
                        ownedTransitionMatrices[t],
                        adjointWorkspace.getDimension());
                GaussianMatrixOps.copyMatrixToFlat(
                        transitionCovarianceWorkspace,
                        transitionCovarianceFlatWorkspace,
                        adjointWorkspace.getDimension());
                CanonicalTransitionAdjointUtils.fillFromMoments(
                        trajectory.transitions[t].precisionYY,
                        transitionCovarianceFlatWorkspace,
                        ownedTransitionMatrices[t],
                        transitionOffsetWorkspace,
                        workingContribution,
                        adjointWorkspace,
                        adjoints[t]);
                transitionMatrices[t] = ownedTransitionMatrices[t];
            }
        }
        known = true;
    }

    public CanonicalBranchMessageContribution getContribution(final int branchIndex) {
        if (!known || knownTrajectory == null) {
            throw new IllegalStateException("Canonical branch gradient cache is not initialized");
        }
        if (contributions[branchIndex] == null) {
            contributions[branchIndex] =
                    new CanonicalBranchMessageContribution(adjoints[branchIndex].getDimension());
        }
        if (!contributionKnown[branchIndex]) {
            fillContribution(branchIndex, knownTrajectory, contributions[branchIndex]);
            contributionKnown[branchIndex] = true;
        }
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

    private void fillContribution(final int branchIndex,
                                  final CanonicalForwardTrajectory trajectory,
                                  final CanonicalBranchMessageContribution out) {
        if (sharedSchedule == null
                || !sharedSchedule.fillContribution(
                branchIndex,
                trajectory.branchPairStates[branchIndex],
                contributionWorkspace,
                out)) {
            CanonicalBranchMessageContributionUtils.fillFromPairState(
                    trajectory.branchPairStates[branchIndex],
                    contributionWorkspace,
                    out);
        }
    }
}
