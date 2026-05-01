/**
 * Pure flat-array linear algebra and Gaussian conversion utilities.
 *
 * <p>Primary APIs use row-major {@code double[]} storage. Boundary helpers convert
 * BEAST {@code double[][]} and EJML matrices into flat form. New canonical
 * hot-path code should use this package instead of
 * {@code dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps}.</p>
 */
package dr.evomodel.treedatalikelihood.continuous.canonical.math;
