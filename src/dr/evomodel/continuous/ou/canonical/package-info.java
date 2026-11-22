/**
 * Canonical OU kernel contracts and backend capability interfaces.
 *
 * <p>This package is the narrow numerical surface between the core OU process
 * model and canonical tree likelihood code. It contains transition kernels,
 * prepared-branch handles, and optional native gradient capabilities. BEAST
 * model ownership and covariance strategy selection remain in
 * {@code dr.evomodel.continuous.ou}; specialized implementation packages remain
 * behind these capability interfaces.</p>
 */
package dr.evomodel.continuous.ou.canonical;
