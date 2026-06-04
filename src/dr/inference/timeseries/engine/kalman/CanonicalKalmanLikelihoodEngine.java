package dr.inference.timeseries.engine.kalman;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;

import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.engine.LikelihoodEngine;
import dr.inference.timeseries.model.gaussian.LinearGaussianObservationModel;
import dr.inference.timeseries.representation.CachedGaussianTransitionRepresentation;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianBranchTransitionKernel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;

/**
 * Forward-only Gaussian likelihood engine that propagates canonical-form states.
 *
 * <p>This engine is intentionally parallel to {@link MomentKalmanLikelihoodEngine}: it
 * computes the same marginal likelihood, but uses precision/information algebra
 * throughout the latent-state propagation. It is the first step toward a selectable
 * canonical backend without disturbing the working moment-form smoother and
 * gradient code.
 */
public final class CanonicalKalmanLikelihoodEngine implements LikelihoodEngine {

    private final CanonicalGaussianBranchTransitionKernel transitionKernel;
    private final CachedGaussianTransitionRepresentation cachedTransitionRepresentation;
    private final TimeGrid timeGrid;

    private final int stateDimension;
    private final int timeCount;

    private final CanonicalGaussianState filteredState;
    private final CanonicalGaussianState predictedState;
    private final CanonicalGaussianTransition transition;
    private final CanonicalGaussianMessageOps.Workspace messageWorkspace;
    private final CanonicalObservationUpdate observationUpdate;

    private boolean likelihoodKnown;
    private double logLikelihood;

    public CanonicalKalmanLikelihoodEngine(final CanonicalGaussianBranchTransitionKernel transitionKernel,
                                           final LinearGaussianObservationModel observationModel,
                                           final TimeGrid timeGrid) {
        this(transitionKernel, null, observationModel, timeGrid);
    }

    public CanonicalKalmanLikelihoodEngine(final CanonicalGaussianBranchTransitionKernel transitionKernel,
                                           final CachedGaussianTransitionRepresentation cachedTransitionRepresentation,
                                           final LinearGaussianObservationModel observationModel,
                                           final TimeGrid timeGrid) {
        if (transitionKernel == null) {
            throw new IllegalArgumentException("transitionKernel must not be null");
        }
        if (observationModel == null) {
            throw new IllegalArgumentException("observationModel must not be null");
        }
        if (timeGrid == null) {
            throw new IllegalArgumentException("timeGrid must not be null");
        }
        this.transitionKernel = transitionKernel;
        this.cachedTransitionRepresentation = cachedTransitionRepresentation;
        this.timeGrid = timeGrid;
        this.stateDimension = transitionKernel.getStateDimension();
        this.timeCount = timeGrid.getTimeCount();

        if (observationModel.getTimeCount() != timeCount) {
            throw new IllegalArgumentException(
                    "Observation time count (" + observationModel.getTimeCount() +
                            ") must match time grid count (" + timeCount + ")");
        }

        this.filteredState = new CanonicalGaussianState(stateDimension);
        this.predictedState = new CanonicalGaussianState(stateDimension);
        this.transition = new CanonicalGaussianTransition(stateDimension);
        this.messageWorkspace = new CanonicalGaussianMessageOps.Workspace(stateDimension);
        this.observationUpdate = new CanonicalObservationUpdate(observationModel, stateDimension);
        if (cachedTransitionRepresentation != null) {
            cachedTransitionRepresentation.prepareTimeGrid(timeGrid);
        }
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = computeLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }

    private double computeLogLikelihood() {
        observationUpdate.refreshStaticMatrices();

        transitionKernel.fillInitialCanonicalState(filteredState);

        double value = 0.0;
        for (int timeIndex = 0; timeIndex < timeCount; ++timeIndex) {
            if (timeIndex == 0) {
                copyState(filteredState, predictedState);
            } else {
                fillCanonicalTransition(timeIndex - 1, timeIndex, transition);
                predict(filteredState, transition, predictedState);
            }

            if (!observationUpdate.prepareForTime(timeIndex)) {
                copyState(predictedState, filteredState);
                continue;
            }

            value += observationUpdate.applyTo(predictedState, filteredState, messageWorkspace);
        }
        return value;
    }

    private void predict(final CanonicalGaussianState previous,
                         final CanonicalGaussianTransition transition,
                         final CanonicalGaussianState out) {
        CanonicalGaussianMessageOps.pushForward(previous, transition, messageWorkspace, out);
    }

    private double validatedDelta(final int fromIndex, final int toIndex) {
        final double dt = timeGrid.getDelta(fromIndex, toIndex);
        if (!(dt > 0.0)) {
            throw new IllegalArgumentException("Time increments must be strictly positive");
        }
        return dt;
    }

    private void fillCanonicalTransition(final int fromIndex,
                                         final int toIndex,
                                         final CanonicalGaussianTransition out) {
        if (cachedTransitionRepresentation != null) {
            cachedTransitionRepresentation.getCanonicalTransition(fromIndex, toIndex, timeGrid, out);
        } else {
            transitionKernel.fillCanonicalTransition(validatedDelta(fromIndex, toIndex), out);
        }
    }

    private static void copyState(final CanonicalGaussianState source,
                                  final CanonicalGaussianState target) {
        System.arraycopy(source.precision, 0, target.precision, 0, source.precision.length);
        System.arraycopy(source.information, 0, target.information, 0, source.information.length);
        target.logNormalizer = source.logNormalizer;
    }
}
