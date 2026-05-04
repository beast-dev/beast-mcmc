/**
 * Tree-level canonical OU adjoint preparation and gradient pullbacks.
 *
 * <p>This package owns the conversion from frozen branch-local adjoints into
 * selection, diffusion, stationary-mean, root-mean, and branch-length gradient
 * targets. It should depend on canonical provider/capability contracts, not on
 * BEAST/CDI lifecycle adapters or specialized backend implementation packages.</p>
 */
package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;
