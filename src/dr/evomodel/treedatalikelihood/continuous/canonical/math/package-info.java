/**
 * Pure flat-array linear algebra and Gaussian conversion utilities.
 *
 * <p>This package provides the mathematical primitives that underpin the canonical
 * diffusion framework. All types operate on {@code double[]} arrays in row-major
 * order; no {@code double[][]} appears in any primary method signature.
 *
 * <h3>Main classes</h3>
 * <ul>
 *   <li>{@link dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps} —
 *       stateless flat matrix operations (arithmetic, Cholesky, EJML bridge, conversion helpers).</li>
 *   <li>{@link dr.evomodel.treedatalikelihood.continuous.canonical.math.GaussianFormConverter} —
 *       moment ↔ canonical form conversions with explicit, reusable {@code Workspace} objects.</li>
 * </ul>
 *
 * <h3>Dependency rule</h3>
 * <p>This package may depend on {@code canonical.message} (for {@code CanonicalGaussianState}
 * and {@code CanonicalGaussianTransition}) but must not depend on {@code canonical},
 * {@code canonical.adapter}, or any BEAST model layer.
 */
package dr.evomodel.treedatalikelihood.continuous.canonical.math;
