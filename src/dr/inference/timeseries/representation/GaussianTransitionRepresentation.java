package dr.inference.timeseries.representation;

import dr.inference.timeseries.core.TimeGrid;

/**
 * Computational representation for linear Gaussian state-space dynamics.
 *
 * <h3>Gradient extensibility</h3>
 * Analytical gradients with respect to the selection (drift) matrix A are computed in
 * two parametrisation-independent steps:
 * <ol>
 *   <li>The RTS smoother produces transition-space sensitivities dlogL/dF_t and
 *       dlogL/df_t that do not depend on how F and f are parametrised by A.</li>
 *   <li>{@link #accumulateSelectionGradientFlat} applies the Jacobian d(F_t, f_t)/dA
 *       specific to this representation and accumulates the result. This is the
 *       only point that changes when switching parametrisations, for example from a
 *       first-order Euler approximation (F = I - dt*A) to an exact
 *       matrix-exponential discretisation where the Frechet derivative is used.</li>
 * </ol>
 *
 * <p>Matrices are row-major flattened arrays. This keeps the Kalman and canonical
 * time-series paths on contiguous storage and avoids boundary copies.
 */
public interface GaussianTransitionRepresentation {

    int getStateDimension();

    void getInitialMean(double[] out);

    void getInitialCovarianceFlat(double[] out);

    void getTransitionMatrixFlat(int fromIndex, int toIndex, TimeGrid timeGrid, double[] out);

    void getTransitionOffset(int fromIndex, int toIndex, TimeGrid timeGrid, double[] out);

    void getTransitionCovarianceFlat(int fromIndex, int toIndex, TimeGrid timeGrid, double[] out);

    /**
     * Accumulates the chain-rule contribution of branch [fromIndex -> toIndex] to the
     * gradient dlogL/dA of the selection (drift) matrix A coming from the dependence
     * of the transition matrix F and offset f on A.
     *
     * <p>The caller provides transition-space sensitivities dlogL/dF_t and dlogL/df_t,
     * obtained from the RTS smoother statistics independently of how A parametrises
     * (F_t, f_t). This method applies the Jacobian d(F_t, f_t)/dA and accumulates the
     * result into {@code gradientAccumulator}.
     *
     * <p>For the first-order Euler approximation (F = I - dt*A, f = dt*A*mu):
     * <pre>
     *   dlogL/dA_kl += dt * (-(dlogL/dF)_kl + (dlogL/df)_k * mu_l)
     * </pre>
     * For an exact matrix-exponential discretisation the same interface is implemented
     * with Frechet-derivative machinery without changing any caller.
     *
     * @param fromIndex           start time index of the branch
     * @param toIndex             end time index of the branch
     * @param timeGrid            time discretisation providing the interval width
     * @param dLogL_dF            row-major flattened dlogL/dF_t, shape [d*d] (read-only)
     * @param dLogL_df            dlogL/df_t, shape [d] (read-only)
     * @param gradientAccumulator row-major flattened dlogL/dA, shape [d*d]; accumulated in-place
     */
    void accumulateSelectionGradientFlat(int fromIndex, int toIndex, TimeGrid timeGrid,
                                         double[] dLogL_dF, double[] dLogL_df,
                                         double[] gradientAccumulator);

    /**
     * Accumulates the chain-rule contribution of branch [fromIndex -> toIndex] to the
     * gradient dlogL/dA coming from the dependence of the step noise covariance V on A.
     *
     * <p>For models where V does not depend on A, for example Euler with V = dt*Q,
     * the default no-op is correct. Exact OU implementations override this path with
     * a Lyapunov or Van Loan adjoint.
     *
     * @param fromIndex           start time index of the branch
     * @param toIndex             end time index of the branch
     * @param timeGrid            time discretisation providing the interval width
     * @param dLogL_dV            row-major flattened dlogL/dV_t, shape [d*d] (read-only)
     * @param gradientAccumulator row-major flattened dlogL/dA, shape [d*d]; accumulated in-place
     */
    default void accumulateSelectionGradientFromCovarianceFlat(final int fromIndex,
                                                               final int toIndex,
                                                               final TimeGrid timeGrid,
                                                               final double[] dLogL_dV,
                                                               final double[] gradientAccumulator) {
        // no-op by default
    }

    /**
     * Accumulates the chain-rule contribution of branch [fromIndex -> toIndex] to the
     * gradient with respect to the diffusion matrix Q coming from the dependence of the
     * step covariance V on Q.
     */
    default void accumulateDiffusionGradientFlat(final int fromIndex,
                                                 final int toIndex,
                                                 final TimeGrid timeGrid,
                                                 final double[] dLogL_dV,
                                                 final double[] gradientAccumulator) {
        // no-op by default
    }
}
