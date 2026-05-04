/**
 * Canonical tree state and post-order/pre-order propagation.
 *
 * <p>This package owns tree traversal mechanics only: state storage, upward
 * message propagation, and downward conditioning context. Gradient target math,
 * BEAST/CDI lifecycle, and OU backend specialization belong outside this
 * package.</p>
 */
package dr.evomodel.treedatalikelihood.continuous.canonical.traversal;
