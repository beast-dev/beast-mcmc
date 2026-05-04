package dr.inference.timeseries.engine.gaussian;

import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Strategy interface for computing ∂logL/∂θ for a specific model parameter using the
 * sufficient statistics produced by the Kalman smoother.
 *
 * <h3>Separation of concerns</h3>
 * <ol>
 *   <li><b>Smoother statistics</b> ({@link BranchSmootherStats}): smoothed means,
 *       covariances, and smoother gains — produced once by the RTS backward pass,
 *       independent of any particular parameter.</li>
 *   <li><b>Transition-space sensitivity</b>: ∂logL/∂F_t and ∂logL/∂f_t — derived
 *       inside {@link #computeGradient} from smoother statistics, independent of how
 *       F and f are parametrised by θ.</li>
 *   <li><b>Parameter Jacobian</b>: the chain rule ∂(F_t, f_t)/∂θ — delegated to
 *       {@link GaussianTransitionRepresentation#accumulateSelectionGradient}, which is
 *       the sole place that changes when switching parametrisations (Euler ↔ exact
 *       matrix-exponential ↔ any future discretisation).</li>
 * </ol>
 *
 * <p>Adding gradient support for a new parameter (e.g. diffusion matrix Q or
 * stationary mean μ) requires only a new implementation of this interface; the smoother
 * infrastructure and the {@link AnalyticalKalmanGradientEngine} loop remain unchanged.
 */
public interface GradientFormula {

    /**
     * Returns {@code true} if this formula can compute the gradient with respect to
     * {@code parameter}.
     */
    boolean supportsParameter(Parameter parameter);

    /**
     * Computes the gradient of the log-likelihood with respect to the supported parameter.
     *
     * <p>Implementations should iterate over branches t = 0 … T−2, accumulate
     * per-branch contributions using the smoother statistics and the forward trajectory,
     * and delegate the parameter-specific chain rule to
     * {@link GaussianTransitionRepresentation#accumulateSelectionGradient}.
     *
     * @param smootherStats per-step RTS smoother statistics, indexed 0 … T−1
     * @param trajectory    forward-pass snapshot (filtered/predicted moments, transitions)
     * @param repr          transition representation; provides the parameter Jacobian
     * @param timeGrid      time discretisation
     * @return gradient array of length equal to the supported parameter's dimension
     */
    double[] computeGradient(Parameter parameter,
                             BranchSmootherStats[] smootherStats,
                             ForwardTrajectory trajectory,
                             GaussianTransitionRepresentation repr,
                             TimeGrid timeGrid);
}
