/**
 * Private block-diagonal OU selection backend.
 *
 * <p>This package owns the canonical backend for selection matrices of the form
 * {@code A = R D R^{-1}}, where {@code D} is assembled from scalar and
 * two-by-two stable blocks. The package still contains the orthogonal
 * specialization for charts with {@code R^{-1} = R^T}, but shared transition
 * caching, prepared branch data, covariance assembly, and canonical-message
 * assembly are general block-diagonal code. BEAST XML/model wiring should stay
 * outside this package, and tree traversal should depend only on the canonical
 * capability interfaces exposed by the OU layer. Public callers should not
 * depend on this package as an architectural boundary.</p>
 */
package dr.evomodel.continuous.ou.blockdiagonal;
