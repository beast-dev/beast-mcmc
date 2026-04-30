/**
 * BEAST-facing adapters for the canonical OU tree likelihood path.
 *
 * <p>This package owns model lifecycle, XML-facing gradient adapters, tip-data
 * synchronization, root-prior adaptation, transition-cache diagnostics, and
 * conversion from BEAST diffusion precision parameters into canonical OU
 * covariance snapshots. Numerical traversal and gradient accumulation live in
 * the integration package; reusable OU transition math lives in
 * {@code dr.evomodel.continuous.ou}.</p>
 */
package dr.evomodel.treedatalikelihood.continuous.adapter;
