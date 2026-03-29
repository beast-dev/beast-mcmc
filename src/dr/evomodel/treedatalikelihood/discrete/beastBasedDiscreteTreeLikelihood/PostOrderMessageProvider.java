package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood;

/**
 * Narrow interface exposing post-order branch-side messages in the STANDARD basis.
 *
 * This keeps the pre-order machinery independent from the internal post-order
 * representation (standard, spectral, etc.).
 */
/*
 * @author Filippo Monti
 */
public interface PostOrderMessageProvider {

    int getStateCount();

    int getPatternCount();

    int getCategoryCount();

    /**
     * Get the post-order message at the TOP of the branch leading into childNodeNumber,
     * in the standard basis, for a single category/pattern slice.
     *
     * This is the upward message from the child subtree, propagated to the parent side.
     */
    void getPostOrderBranchTopInto(int childNodeNumber, int category, int pattern, double[] outStandardPartial);

    /**
     * Get the post-order message at the BOTTOM of the branch leading into childNodeNumber,
     * in the standard basis, for a single category/pattern slice.
     *
     * This is the child-node post-order message before propagation along the branch.
     */
    void getPostOrderBranchBottomInto(int childNodeNumber, int category, int pattern, double[] outStandardPartial);

    void getPostOrderBranchScalesInto(int nodeNumber, double[] dest);
}
