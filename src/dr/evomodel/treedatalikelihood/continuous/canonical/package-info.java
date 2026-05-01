/**
 * Canonical tree integration internals.
 *
 * <p>The package is organized around a fixed pipeline:</p>
 *
 * <ol>
 *   <li>Project tip observations into canonical Gaussian messages.</li>
 *   <li>Run post-order traversal from tips to root.</li>
 *   <li>Run pre-order traversal from root to descendants.</li>
 *   <li>Assemble branch-local parent/child contributions.</li>
 *   <li>Freeze branch adjoints in {@code BranchGradientInputs}.</li>
 *   <li>Pull frozen adjoints back to selection, diffusion, stationary mean,
 *       root mean, or branch length targets.</li>
 * </ol>
 *
 * <p>Traversal code should stay separate from gradient target code. Shared
 * provider checks, diagnostic phase names, and fallback policy belong in small
 * helpers so the message passer remains orchestration.</p>
 */
package dr.evomodel.treedatalikelihood.continuous.canonical;
