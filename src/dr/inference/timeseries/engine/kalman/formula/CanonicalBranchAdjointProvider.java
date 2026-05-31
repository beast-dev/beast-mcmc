package dr.inference.timeseries.engine.kalman.formula;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
import dr.inference.timeseries.engine.kalman.CanonicalBranchGradientCache;
import dr.inference.timeseries.engine.kalman.CanonicalForwardTrajectory;

/**
 * Branch-local canonical adjoint access with an optional shared cache.
 */
final class CanonicalBranchAdjointProvider {

    private final CanonicalBranchMessageContribution localContribution;
    private final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
    private final CanonicalLocalTransitionAdjoints localAdjoints;
    private final CanonicalTransitionAdjointUtils.Workspace transitionWorkspace;

    CanonicalBranchAdjointProvider(final int stateDimension) {
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        this.localContribution = new CanonicalBranchMessageContribution(stateDimension);
        this.contributionWorkspace = new CanonicalBranchMessageContributionUtils.Workspace(stateDimension);
        this.localAdjoints = new CanonicalLocalTransitionAdjoints(stateDimension);
        this.transitionWorkspace = new CanonicalTransitionAdjointUtils.Workspace(stateDimension);
    }

    void ensure(final CanonicalForwardTrajectory trajectory,
                final CanonicalBranchGradientCache branchGradientCache) {
        if (branchGradientCache != null) {
            branchGradientCache.ensure(trajectory);
        }
    }

    CanonicalBranchMessageContribution localContribution(final int branchIndex,
                                                         final CanonicalForwardTrajectory trajectory,
                                                         final CanonicalBranchGradientCache branchGradientCache) {
        if (branchGradientCache != null) {
            return branchGradientCache.getContribution(branchIndex);
        }
        CanonicalBranchMessageContributionUtils.fillFromPairState(
                trajectory.branchPairStates[branchIndex],
                contributionWorkspace,
                localContribution);
        return localContribution;
    }

    CanonicalLocalTransitionAdjoints localAdjoints(final int branchIndex,
                                                   final CanonicalForwardTrajectory trajectory,
                                                   final CanonicalBranchGradientCache branchGradientCache) {
        if (branchGradientCache != null) {
            return branchGradientCache.getAdjoints(branchIndex);
        }
        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                trajectory.transitions[branchIndex],
                localContribution(branchIndex, trajectory, null),
                transitionWorkspace,
                localAdjoints);
        return localAdjoints;
    }

    double[] transitionMatrixData(final int branchIndex,
                                  final CanonicalBranchGradientCache branchGradientCache) {
        return branchGradientCache == null
                ? transitionWorkspace.transitionMatrix
                : branchGradientCache.getTransitionMatrixData();
    }

    int transitionMatrixOffset(final int branchIndex,
                               final CanonicalBranchGradientCache branchGradientCache) {
        return branchGradientCache == null
                ? 0
                : branchGradientCache.getTransitionMatrixOffset(branchIndex);
    }
}
