package dr.inference.timeseries.engine.gaussian;

import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Canonical-path analogue of {@link GradientFormula}.
 */
public interface CanonicalGradientFormula {

    boolean supportsParameter(Parameter parameter);

    double[] computeGradient(Parameter parameter,
                             CanonicalForwardTrajectory trajectory,
                             GaussianTransitionRepresentation repr,
                             TimeGrid timeGrid);
}
