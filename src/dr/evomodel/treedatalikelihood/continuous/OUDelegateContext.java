package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.MultivariateElasticModel;

/**
 * Minimal context that an OU actualization strategy needs from the delegate.
 * The delegate (e.g. OUDiffusionModelDelegate) should implement this.
 */
public interface OUDelegateContext {

    Tree getTree();

    MultivariateElasticModel getElasticModel();

    int getDim();

    /** Map a node index to the matrix-buffer offset used by the CDI. */
    int getMatrixBufferOffsetIndex(int nodeIndex);

    /** Map an eigen-buffer index (usually 0) to the CDI eigen buffer offset. */
    int getEigenBufferOffsetIndex(int eigenIndex);

    /** Per-node drift rate (length dim). */
    double[] getDriftRate(NodeRef node);
}

