package dr.inference.timeseries.representation;

import dr.inference.timeseries.core.TimeGrid;

/**
 * Computational representation for linear Gaussian state-space dynamics.
 *
 * <h3>Gradient extensibility</h3>
 * Analytical gradients with respect to the selection (drift) matrix A are computed in
 * two parametrisation-independent steps:
 * <ol>
 *   <li>The RTS smoother produces transition-space sensitivities ∂logL/∂F_t and
 *       ∂logL/∂f_t that do not depend on how F and f are parametrised by A.</li>
 *   <li>{@link #accumulateSelectionGradient} applies the Jacobian ∂(F_t, f_t)/∂A
 *       specific to this representation and accumulates the result. This is the
 *       <em>only</em> point that changes when switching parametrisations — for example
 *       from a first-order Euler approximation (F = I − dt·A) to an exact
 *       matrix-exponential discretisation (F = expm(−dt·A)), where the Fréchet
 *       derivative of the matrix exponential is used instead.</li>
 * </ol>
 */
public interface GaussianTransitionRepresentation {

    int getStateDimension();

    void getInitialMean(double[] out);

    void getInitialCovariance(double[][] out);

    void getTransitionMatrix(int fromIndex, int toIndex, TimeGrid timeGrid, double[][] out);

    void getTransitionOffset(int fromIndex, int toIndex, TimeGrid timeGrid, double[] out);

    void getTransitionCovariance(int fromIndex, int toIndex, TimeGrid timeGrid, double[][] out);

    /**
     * Accumulates the chain-rule contribution of branch [fromIndex → toIndex] to the
     * gradient ∂logL/∂A of the selection (drift) matrix A coming from the
     * dependence of the transition matrix F and offset f on A.
     *
     * <p>The caller provides the transition-space sensitivities ∂logL/∂F_t and
     * ∂logL/∂f_t, obtained from the RTS smoother statistics independently of how A
     * parametrises (F_t, f_t). This method applies the Jacobian ∂(F_t, f_t)/∂A and
     * accumulates the result into {@code gradientAccumulator}.
     *
     * <p>For the first-order Euler approximation (F = I − dt·A, f = dt·A·μ):
     * <pre>
     *   ∂logL/∂A_{kl}  +=  dt · (−(∂logL/∂F)_{kl}  +  (∂logL/∂f)_k · μ_l)
     * </pre>
     * For an exact matrix-exponential discretisation the same interface is implemented
     * with Fréchet-derivative machinery without changing any caller.
     *
     * @param fromIndex           start time index of the branch
     * @param toIndex             end time index of the branch
     * @param timeGrid            time discretisation providing the interval width
     * @param dLogL_dF            ∂logL/∂F_t, shape [d × d] (read-only)
     * @param dLogL_df            ∂logL/∂f_t, shape [d]     (read-only)
     * @param gradientAccumulator row-major flattened ∂logL/∂A, shape [d·d]; accumulated in-place
     */
    void accumulateSelectionGradient(int fromIndex, int toIndex, TimeGrid timeGrid,
                                     double[][] dLogL_dF, double[] dLogL_df,
                                     double[] gradientAccumulator);

    /**
     * Accumulates the chain-rule contribution of branch [fromIndex → toIndex] to the
     * gradient ∂logL/∂A coming from the dependence of the step noise covariance V on A.
     *
     * <p>The caller provides ∂logL/∂V_t computed from the RTS smoother sufficient
     * statistics (expected squared prediction-error second moment). This method applies
     * the Jacobian ∂V_t/∂A and accumulates the result.
     *
     * <p>For models where V does not depend on A (e.g. Euler: V = dt·Q), the default
     * no-op is correct.  For exact matrix-exponential discretisations, the Van Loan
     * block-exponential adjoint is used:
     * <pre>
     *   M = [[-A dt, Q dt], [0, A^T dt]]
     *   expm(M) = [[F, F·V], [0, F^{-T}]]
     *   G_E     = [[0,  F^{-T}·(∂logL/∂V)], [0, 0]]
     *   G_M     = adjointExp(M, G_E)
     *   ∂logL/∂A += -dt·G_M_{top-left} + dt·G_M_{bottom-right}^T
     * </pre>
     *
     * @param fromIndex           start time index of the branch
     * @param toIndex             end time index of the branch
     * @param timeGrid            time discretisation providing the interval width
     * @param dLogL_dV            ∂logL/∂V_t, shape [d × d] (read-only)
     * @param gradientAccumulator row-major flattened ∂logL/∂A, shape [d·d]; accumulated in-place
     */
    default void accumulateSelectionGradientFromCovariance(int fromIndex, int toIndex,
                                                            TimeGrid timeGrid,
                                                            double[][] dLogL_dV,
                                                            double[] gradientAccumulator) {
        // no-op: V does not depend on A for this parametrisation (e.g., Euler: V = dt·Q)
    }

    /**
     * Accumulates the chain-rule contribution of branch [fromIndex → toIndex] to the
     * gradient with respect to the diffusion matrix Q coming from the dependence of the
     * step covariance V on Q.
     */
    default void accumulateDiffusionGradient(int fromIndex,
                                             int toIndex,
                                             TimeGrid timeGrid,
                                             double[][] dLogL_dV,
                                             double[] gradientAccumulator) {
        // no-op by default
    }
}
