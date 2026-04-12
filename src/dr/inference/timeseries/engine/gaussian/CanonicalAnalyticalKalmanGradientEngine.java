package dr.inference.timeseries.engine.gaussian;

import dr.inference.model.Parameter;
import dr.inference.timeseries.engine.GradientEngine;

/**
 * Analytical gradient engine driven by the canonical-form smoother backend.
 *
 * <p>This engine consumes the canonical smoother results directly through
 * {@link CanonicalSmootherResults} and dispatches to {@link CanonicalGradientFormula}
 * implementations that operate on {@link CanonicalForwardTrajectory}.
 */
public final class CanonicalAnalyticalKalmanGradientEngine implements GradientEngine {

    private final CanonicalSmootherResults smootherEngine;
    private final CanonicalGradientFormula[] formulas;

    public CanonicalAnalyticalKalmanGradientEngine(final CanonicalSmootherResults smootherEngine,
                                                   final CanonicalGradientFormula... formulas) {
        if (smootherEngine == null) {
            throw new IllegalArgumentException("smootherEngine must not be null");
        }
        if (formulas == null || formulas.length == 0) {
            throw new IllegalArgumentException("at least one GradientFormula must be provided");
        }
        this.smootherEngine = smootherEngine;
        this.formulas = formulas.clone();
    }

    @Override
    public boolean supportsGradientWrt(final Parameter parameter) {
        if (parameter == null) {
            return false;
        }
        for (final CanonicalGradientFormula formula : formulas) {
            if (formula.supportsParameter(parameter)) {
                return true;
            }
        }
        return false;
    }

    @Override
    public double[] getGradientWrt(final Parameter parameter) {
        if (!supportsGradientWrt(parameter)) {
            throw new IllegalArgumentException(
                    "No registered GradientFormula supports parameter: " +
                            (parameter == null ? "null" : parameter.getId()));
        }

        final CanonicalForwardTrajectory trajectory = smootherEngine.getCanonicalTrajectory();

        for (final CanonicalGradientFormula formula : formulas) {
            if (formula.supportsParameter(parameter)) {
                return formula.computeGradient(
                        parameter,
                        trajectory,
                        smootherEngine.getTransitionRepresentation(),
                        smootherEngine.getTimeGrid());
            }
        }

        throw new AssertionError("unreachable");
    }

    @Override
    public void makeDirty() {
        smootherEngine.makeDirty();
    }
}
