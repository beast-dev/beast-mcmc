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
 * <p>The root package intentionally holds the public facade, transition/root
 * provider contracts, cache diagnostics, and orchestration. Implementation
 * details live in narrower subpackages:</p>
 *
 * <ul>
 *   <li>{@code traversal}: post-order/pre-order tree walks and tree state.</li>
 *   <li>{@code contribution}: branch-local canonical contribution assembly.</li>
 *   <li>{@code gradient}: frozen branch adjoints and target pullbacks.</li>
 *   <li>{@code workspace}: reusable buffers and workspace factories.</li>
 *   <li>{@code scheduling}: parallel chunk sizing policies.</li>
 *   <li>{@code message}, {@code math}, and {@code adapter}: canonical Gaussian
 *       algebra, matrix conversion, and model/integrator adapters.</li>
 * </ul>
 *
 * <p>Tree-layer code in this package must talk to specialized OU backends only
 * through {@code dr.evomodel.continuous.ou.canonical} capability interfaces.
 * Backend implementation packages such as {@code orthogonalblockdiagonal}
 * remain behind those capabilities.</p>
 */
package dr.evomodel.treedatalikelihood.continuous.canonical;
