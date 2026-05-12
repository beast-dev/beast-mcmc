package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood;

/**
 * Narrow interface exposing post-order branch-side messages.
 *
 * Algorithmic accessors return messages in the post-order representation's
 * internal coordinates. Export accessors are for reporting/debugging.
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
     * in the post-order representation's internal basis, for a single
     * category/pattern slice.
     *
     * This is the upward message from the child subtree, propagated to the parent side.
     */
    void getPostOrderBranchTopInto(int childNodeNumber, int category, int pattern, double[] outPartial);

    /**
     * Get the post-order message at the TOP of the branch in standard data-type
     * coordinates. Traversals that combine messages in standard basis can use
     * this to avoid repeated internal-to-standard conversion.
     */
    void getPostOrderBranchTopStandardInto(int childNodeNumber, int category, int pattern, double[] outPartial);

    /**
     * Get the post-order message at the BOTTOM of the branch leading into childNodeNumber,
     * in exported coordinates, for a single category/pattern slice.
     *
     * This is the child-node post-order message before propagation along the branch.
     */
    void getPostOrderBranchBottomInto(int childNodeNumber, int category, int pattern, double[] outPartial);

    /**
     * Get the post-order message at the TOP of the branch in exported coordinates.
     * This is intended for report/export code, not traversal coordination.
     */
    void getPostOrderBranchTopExportInto(int childNodeNumber, int category, int pattern, double[] outPartial);

    void getPostOrderBranchScalesInto(int nodeNumber, double[] dest);
}
