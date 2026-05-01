/**
 * Private orthogonal block-diagonal OU selection backend.
 *
 * <p>This package contains the specialized canonical backend for selection
 * matrices of the form {@code R D R^T}, where {@code R} is parameterized by an
 * orthogonal chart and {@code D} is assembled from scalar and two-by-two stable
 * blocks. The package owns basis caching, prepared branch data, transition and
 * covariance assembly, and native gradient pullbacks. BEAST XML/model wiring
 * should stay outside this package, and tree traversal should depend only on
 * the canonical capability interfaces exposed by the OU layer. Public callers
 * should not depend on this package as an architectural boundary.</p>
 */
package dr.evomodel.continuous.ou.orthogonalblockdiagonal;
