/**
 * Canonical tree integration internals.
 *
 * <p>The package is organized around four responsibilities:
 * traversal builds post-order and pre-order canonical messages; state storage
 * owns the reusable tree buffers and MCMC lifecycle; branch contribution
 * assembly converts local parent/child messages into transition adjoints; and
 * gradient engines pull those branch adjoints back to OU parameters.</p>
 */
package dr.evomodel.treedatalikelihood.continuous.integration;
