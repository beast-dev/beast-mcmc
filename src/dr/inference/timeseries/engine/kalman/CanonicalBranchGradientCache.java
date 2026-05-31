package dr.inference.timeseries.engine.kalman;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
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
    private final double[] transitionMatrices;
    private final GaussianTransitionRepresentation transitionRepresentation;
    private final CachedGaussianTransitionRepresentation cachedTransitionRepresentation;
    private final TimeGrid timeGrid;
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
    private final int matrixLength;

    CanonicalBranchGradientCache(final int timeCount,
                                 final int stateDimension,
                                 final GaussianTransitionRepresentation transitionRepresentation,
                                 final TimeGrid timeGrid) {
        final int branchCount = Math.max(0, timeCount - 1);
        this.matrixLength = stateDimension * stateDimension;
        this.contributions = new CanonicalBranchMessageContribution[branchCount];
        this.contributionKnown = new boolean[branchCount];
        this.adjoints = new CanonicalLocalTransitionAdjoints[branchCount];
        this.transitionMatrices = new double[branchCount * matrixLength];
        this.transitionRepresentation = transitionRepresentation;
        this.cachedTransitionRepresentation =
                transitionRepresentation instanceof CachedGaussianTransitionRepresentation
                        ? (CachedGaussianTransitionRepresentation) transitionRepresentation
                        : null;
        this.timeGrid = timeGrid;
        this.transitionCovarianceFlatWorkspace = new double[matrixLength];
        this.transitionOffsetWorkspace = new double[stateDimension];
        this.transitionMomentsView = new TransitionMomentsView();
        this.workingContribution = new CanonicalBranchMessageContribution(stateDimension);
        for (int i = 0; i < branchCount; ++i) {
            adjoints[i] = new CanonicalLocalTransitionAdjoints(stateDimension);
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
            final int transitionMatrixOffset = transitionMatrixOffset(t);
            fillContribution(t, trajectory, workingContribution);
            if (transitionRepresentation == null || timeGrid == null) {
                CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                        trajectory.transitions[t],
                        workingContribution,
                        adjointWorkspace,
                        adjoints[t]);
                System.arraycopy(adjointWorkspace.transitionMatrix, 0,
                        transitionMatrices, transitionMatrixOffset, matrixLength);
            } else if (cachedTransitionRepresentation != null
                    && cachedTransitionRepresentation.getTransitionMomentsView(
                    t,
                    t + 1,
                    timeGrid,
                    transitionMomentsView)) {
                System.arraycopy(transitionMomentsView.getTransitionMatrix(), 0,
                        transitionMatrices, transitionMatrixOffset, matrixLength);
                CanonicalTransitionAdjointUtils.fillFromMoments(
                        trajectory.transitions[t].precisionYY,
                        transitionMomentsView.getTransitionCovariance(),
                        transitionMomentsView.getTransitionMatrix(),
                        transitionMomentsView.getTransitionOffset(),
                        workingContribution,
                        adjointWorkspace,
                        adjoints[t]);
            } else {
                transitionRepresentation.getTransitionMatrixFlat(
                        t, t + 1, timeGrid, adjointWorkspace.transitionMatrix);
                System.arraycopy(adjointWorkspace.transitionMatrix, 0,
                        transitionMatrices, transitionMatrixOffset, matrixLength);
                transitionRepresentation.getTransitionOffset(t, t + 1, timeGrid, transitionOffsetWorkspace);
                transitionRepresentation.getTransitionCovarianceFlat(
                        t, t + 1, timeGrid, transitionCovarianceFlatWorkspace);
                CanonicalTransitionAdjointUtils.fillFromMoments(
                        trajectory.transitions[t].precisionYY,
                        transitionCovarianceFlatWorkspace,
                        adjointWorkspace.transitionMatrix,
                        transitionOffsetWorkspace,
                        workingContribution,
                        adjointWorkspace,
                        adjoints[t]);
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

    public double[] getTransitionMatrixData() {
        return transitionMatrices;
    }

    public int getTransitionMatrixOffset(final int branchIndex) {
        return transitionMatrixOffset(branchIndex);
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

    private int transitionMatrixOffset(final int branchIndex) {
        return branchIndex * matrixLength;
    }
}
