/**
 * Canonical Gaussian message primitives and numerics.
 *
 * <p>This package owns the low-level canonical representation used by the tree
 * traversal: precision/information states, branch transition factors, local
 * transition adjoints, and Gaussian message algebra. Classes here should remain
 * tree-topology agnostic; traversal, branch selection, BEAST model lifecycle, and
 * general matrix utilities belong in parent/sibling packages. New canonical
 * hot-path matrix code should use
 * {@code dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps}.</p>
 */
package dr.evomodel.treedatalikelihood.continuous.canonical.message;
