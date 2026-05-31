package dr.inference.timeseries.engine.kalman;

import dr.inference.model.Parameter;
import dr.inference.timeseries.engine.GradientEngine;
import dr.inference.timeseries.engine.kalman.formula.ExpectationGradientFormula;

/**
 * Analytical gradient engine for the linear-Gaussian Kalman path.
 *
 * <p>Gradients are computed by combining an RTS smoother pass (via
 * {@link ExpectationKalmanSmootherEngine}) with one or more {@link ExpectationGradientFormula} implementations
 * that accumulate per-branch contributions. The smoother result is shared across all
 * formulas registered at construction time, so the O(T·d³) forward + backward pass is
 * paid only once per parameter update regardless of how many gradient parameters are
 * queried.
 *
 * <h3>Adding support for a new parameter</h3>
 * Implement {@link ExpectationGradientFormula} and pass an instance to the constructor. No other
 * changes are required.
 *
 * <h3>Caching</h3>
 * The smoother engine caches its results behind a shared dirty flag. Calling
 * {@link #makeDirty()} invalidates both the likelihood and the smoother statistics so
 * that the next gradient query triggers a fresh pass.
 */
public class ExpectationAnalyticalKalmanGradientEngine implements GradientEngine {

    private final GaussianSmootherResults smootherEngine;
    private final ExpectationGradientFormula[]    formulas;

    /**
     * @param smootherEngine the smoother engine providing shared forward + backward statistics
     * @param formulas       one or more gradient formulas to register; queried in order for
     *                       {@link #supportsGradientWrt}
     */
    public ExpectationAnalyticalKalmanGradientEngine(final ExpectationKalmanSmootherEngine smootherEngine,
                                          final ExpectationGradientFormula... formulas) {
        this((GaussianSmootherResults) smootherEngine, formulas);
    }

    public ExpectationAnalyticalKalmanGradientEngine(final GaussianSmootherResults smootherEngine,
                                          final ExpectationGradientFormula... formulas) {
        if (smootherEngine == null) {
            throw new IllegalArgumentException("smootherEngine must not be null");
        }
        if (formulas == null || formulas.length == 0) {
            throw new IllegalArgumentException("at least one ExpectationGradientFormula must be provided");
        }
        this.smootherEngine = smootherEngine;
        this.formulas       = formulas.clone();
    }

    @Override
    public boolean supportsGradientWrt(final Parameter parameter) {
        if (parameter == null) {
            return false;
        }
        for (final ExpectationGradientFormula formula : formulas) {
            if (formula.supportsParameter(parameter)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Computes the gradient of the log-likelihood with respect to {@code parameter}.
     *
     * <p>The smoother statistics are computed (or retrieved from cache) first; then the
     * matching {@link ExpectationGradientFormula} runs its branch-sum loop over the stored
     * trajectory and smoother stats.
     *
     * @throws IllegalArgumentException if no registered formula supports {@code parameter}
     */
    @Override
    public double[] getGradientWrt(final Parameter parameter) {
        if (!supportsGradientWrt(parameter)) {
            throw new IllegalArgumentException(
                    "No registered ExpectationGradientFormula supports parameter: " +
                    (parameter == null ? "null" : parameter.getId()));
        }

        // One call ensures both forward trajectory and smoother stats are up to date.
        final BranchSmootherStats[] stats       = smootherEngine.getSmootherStats();
        final ForwardTrajectory     trajectory  = smootherEngine.getTrajectory();

        for (final ExpectationGradientFormula formula : formulas) {
            if (formula.supportsParameter(parameter)) {
                return formula.computeGradient(
                        parameter,
                        stats,
                        trajectory,
                        smootherEngine.getTransitionRepresentation(),
                        smootherEngine.getTimeGrid());
            }
        }

        // Unreachable given the supportsGradientWrt guard above.
        throw new AssertionError("unreachable");
    }

    @Override
    public void prepareGradient() {
        smootherEngine.getTrajectory();
    }

    @Override
    public void makeDirty() {
        smootherEngine.makeDirty();
    }
}
