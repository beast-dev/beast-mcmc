/**
 * Canonical Gaussian message primitives and numerics.
 *
 * <p>This package owns the low-level canonical representation used by the tree
 * traversal: precision/information states, branch transition factors, local
 * transition adjoints, and matrix utilities needed to combine, condition, and
 * convert Gaussian messages. Classes here should remain tree-topology agnostic;
 * traversal, branch selection, and BEAST model lifecycle belong in parent
 * packages.</p>
 */
package dr.evomodel.treedatalikelihood.continuous.canonical.message;
