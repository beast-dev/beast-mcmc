package dr.inference.timeseries.engine.gaussian;

import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.engine.LikelihoodEngine;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianBranchTransitionKernel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import dr.inference.timeseries.representation.RepresentableProcess;

/**
 * Factory for forward Gaussian likelihood engines.
 *
 * <p>This centralizes the choice between expectation-form and canonical-form
 * forward inference while leaving smoothing and analytical gradients on the
 * established expectation-form path for now.
 */
public final class GaussianLikelihoodEngineFactory {

    private GaussianLikelihoodEngineFactory() {
        // no instances
    }

    public static LikelihoodEngine createForwardEngine(final RepresentableProcess process,
                                                       final GaussianObservationModel observationModel,
                                                       final TimeGrid timeGrid,
                                                       final GaussianForwardComputationMode mode) {
        if (process == null) {
            throw new IllegalArgumentException("process must not be null");
        }
        if (observationModel == null) {
            throw new IllegalArgumentException("observationModel must not be null");
        }
        if (timeGrid == null) {
            throw new IllegalArgumentException("timeGrid must not be null");
        }
        if (mode == null) {
            throw new IllegalArgumentException("mode must not be null");
        }

        switch (mode) {
            case EXPECTATION:
                return new KalmanLikelihoodEngine(
                        process.getRepresentation(GaussianTransitionRepresentation.class),
                        observationModel,
                        timeGrid);
            case CANONICAL:
                if (!process.supportsRepresentation(CanonicalGaussianBranchTransitionKernel.class)) {
                    throw new IllegalArgumentException(
                            "Process does not support canonical Gaussian branch transitions: "
                                    + process.getClass().getName());
                }
                return new CanonicalKalmanLikelihoodEngine(
                        process.getRepresentation(CanonicalGaussianBranchTransitionKernel.class),
                        observationModel,
                        timeGrid);
            default:
                throw new IllegalArgumentException("Unsupported forward computation mode: " + mode);
        }
    }
}
