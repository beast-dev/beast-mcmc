/**
 * Small public contracts for the canonical OU tree likelihood path.
 *
 * <p>Classes outside the integration package should depend on these interfaces
 * rather than concrete traversal, cache, or gradient implementations. This
 * package is the boundary for transition providers, root priors, tip
 * observations, prepared branch snapshots, and cache diagnostic phases.</p>
 */
package dr.evomodel.treedatalikelihood.continuous.framework;
