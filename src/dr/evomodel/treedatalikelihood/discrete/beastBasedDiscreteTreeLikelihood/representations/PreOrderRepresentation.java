package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.representations;

import dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.DiscretePreOrderDelegate;

/**
 * Defines the internal algebra of pre-order partials used by
 * {@link DiscretePreOrderDelegate}.
 *
 * The internal representation may be:
 * - standard basis: q
 * - spectral dual basis: y = R^T q
 * - etc.
 *
 * All methods are allocation-free with caller-owned buffers.
 */
/*
 * @author Filippo Monti
 */
public interface PreOrderRepresentation {

    String getName();

    int getStateCount();

    void markDirty();

    void storeState();

    void restoreState();

    void updateForLikelihood();

    /**
     * Initialize the root pre-order message from root frequencies (standard basis).
     */
    void initializeRootPartial(double[] rootFrequencies, double[] outRootPreOrder);

    /**
     * Combine the parent-node pre-order message with the sibling branch-top post-order
     * message to obtain the pre-order message at the TOP of the child's branch.
     *
     * All inputs/outputs are in the internal pre-order representation except the sibling
     * post-order message, which is given in standard basis.
     */
    void combineParentAndSibling(double[] parentNodePreOrder,
                                 double[] siblingBranchTopPostOrderStandard,
                                 double[] outChildBranchTopPreOrder);

    /**
     * Propagate a pre-order message from the TOP of the child's branch to the BOTTOM
     * of the child's branch (i.e. to the child node).
     */
    void propagateToBranchBottom(int childNodeNumber,
                                 double branchLength,
                                 double[] childBranchTopPreOrder,
                                 double[] outChildNodePreOrder);

    /**
     * Convert one internal pre-order slice to standard basis.
     */
    void toStandard(double[] preOrderPartial, double[] outStandardPartial);
}
