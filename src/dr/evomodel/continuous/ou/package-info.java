/**
 * Ornstein-Uhlenbeck process kernels and selection-parameterization boundaries.
 *
 * <p>{@code OUProcessModel} is the BEAST model wrapper: it owns live variables,
 * model events, and time-series representation compatibility. Canonical tree
 * code should prefer the narrower
 * {@code dr.evomodel.continuous.ou.canonical.CanonicalOUKernel} and capability
 * interfaces when it only needs transition, covariance, or gradient math.</p>
 *
 * <p>Flat matrices in canonical hot paths use row-major order:
 * {@code index = row * dimension + column}. Conversions to BEAST parameter
 * layout should be explicit at adapter boundaries.</p>
 */
package dr.evomodel.continuous.ou;
