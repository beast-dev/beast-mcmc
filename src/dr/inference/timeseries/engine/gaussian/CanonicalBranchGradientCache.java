package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;

/**
 * Per-smoother cache of branch-local canonical adjoints shared by gradient formulas.
 */
public final class CanonicalBranchGradientCache {

    private final CanonicalBranchMessageContribution[] contributions;
    private final CanonicalLocalTransitionAdjoints[] adjoints;
    private final double[][] transitionMatrices;
    private final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
    private final CanonicalTransitionAdjointUtils.Workspace adjointWorkspace;
    private boolean known;
    private long buildCount;

    CanonicalBranchGradientCache(final int timeCount, final int stateDimension) {
        final int branchCount = Math.max(0, timeCount - 1);
        this.contributions = new CanonicalBranchMessageContribution[branchCount];
        this.adjoints = new CanonicalLocalTransitionAdjoints[branchCount];
        this.transitionMatrices = new double[branchCount][stateDimension * stateDimension];
        for (int i = 0; i < branchCount; ++i) {
            contributions[i] = new CanonicalBranchMessageContribution(stateDimension);
            adjoints[i] = new CanonicalLocalTransitionAdjoints(stateDimension);
        }
        this.contributionWorkspace = new CanonicalBranchMessageContributionUtils.Workspace(stateDimension);
        this.adjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(stateDimension);
        this.known = false;
        this.buildCount = 0L;
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
            CanonicalBranchMessageContributionUtils.fillFromPairState(
                    trajectory.branchPairStates[t],
                    contributionWorkspace,
                    contributions[t]);
            CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                    trajectory.transitions[t],
                    contributions[t],
                    adjointWorkspace,
                    adjoints[t]);
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
