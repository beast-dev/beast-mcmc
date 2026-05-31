package dr.evomodel.continuous.ou.canonical;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.inference.model.MatrixParameterInterface;

/**
 * Shared access point for prepared canonical OU transition implementations.
 *
 * <p>Tree and time-series caches have different invalidation keys, but once a
 * cache decides to use a prepared branch they perform the same OU-layer calls:
 * create a handle, create a branch workspace, prepare for a branch length, and
 * fill either canonical parameters or transition moments.</p>
 */
public final class CanonicalPreparedTransitionSupport {

    private final OUProcessModel processModel;
    private final CanonicalPreparedTransitionCapability preparedTransition;
    private final ThreadLocal<CanonicalBranchWorkspace> workspaces;

    public static CanonicalPreparedTransitionSupport create(final OUProcessModel processModel) {
        if (processModel == null
                || !(processModel.getSelectionMatrixParameterization()
                instanceof CanonicalPreparedTransitionCapability)) {
            return null;
        }
        return new CanonicalPreparedTransitionSupport(
                processModel,
                (CanonicalPreparedTransitionCapability)
                        processModel.getSelectionMatrixParameterization());
    }

    private CanonicalPreparedTransitionSupport(
            final OUProcessModel processModel,
            final CanonicalPreparedTransitionCapability preparedTransition) {
        this.processModel = processModel;
        this.preparedTransition = preparedTransition;
        this.workspaces = ThreadLocal.withInitial(preparedTransition::createBranchWorkspace);
    }

    public CanonicalPreparedBranchHandle createPreparedBranchHandle() {
        return preparedTransition.createPreparedBranchHandle();
    }

    public void prepareBranch(final double branchLength,
                              final double[] stationaryMean,
                              final CanonicalPreparedBranchHandle prepared) {
        preparedTransition.prepareBranch(branchLength, stationaryMean, prepared);
    }

    public void fillCanonicalTransitionPrepared(
            final CanonicalPreparedBranchHandle prepared,
            final CanonicalGaussianTransition out) {
        preparedTransition.fillCanonicalTransitionPrepared(
                prepared,
                processModel.getDiffusionMatrix(),
                workspaces.get(),
                out);
    }

    public boolean fillTransitionMomentsPreparedFlat(
            final CanonicalPreparedBranchHandle prepared,
            final double[] transitionMatrixOut,
            final double[] transitionOffsetOut,
            final double[] transitionCovarianceOut) {
        return preparedTransition.fillTransitionMomentsPreparedFlat(
                prepared,
                processModel.getDiffusionMatrix(),
                workspaces.get(),
                transitionMatrixOut,
                transitionOffsetOut,
                transitionCovarianceOut);
    }

    public MatrixParameterInterface getDiffusionMatrix() {
        return processModel.getDiffusionMatrix();
    }
}
