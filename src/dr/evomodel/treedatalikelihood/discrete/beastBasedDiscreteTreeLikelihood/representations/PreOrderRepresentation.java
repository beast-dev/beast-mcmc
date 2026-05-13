package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.representations;

import dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.DiscretePreOrderDelegate;

/**
 * Defines the internal algebra of pre-order partials used by
 * {@link DiscretePreOrderDelegate}.
 *
 * The internal representation is fully owned by the implementation.
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
     * The parent input/output are in the internal pre-order representation. The sibling
     * input is in the paired post-order representation's internal basis.
     */
    void combineParentAndSibling(double[] parentNodePreOrder,
                                 double[] siblingBranchTopPostOrder,
                                 double[] outChildBranchTopPreOrder);

    /**
     * Whether the internal pre-order partials are already expressed in the
     * standard data-type basis.
     */
    default boolean storesPartialsInStandardBasis() {
        return true;
    }

    /**
     * Convert a standard-basis pre-order partial to this representation's
     * internal basis.
     */
    default void importPreOrderPartialFromStandard(double[] standardPartial, double[] outPreOrderPartial) {
        System.arraycopy(standardPartial, 0, outPreOrderPartial, 0, standardPartial.length);
    }

    /**
     * Export one internal pre-order partial slice to the standard data-type basis.
     */
    default void exportPreOrderPartialToStandard(double[] preOrderPartial, double[] outStandardPartial) {
        exportPreOrderPartial(preOrderPartial, outStandardPartial);
    }

    /**
     * Offset-aware variant — reads src[srcOff..srcOff+K] and writes dst[dstOff..dstOff+K].
     * Avoids round-trip copies when src/dst are flat multi-pattern buffers.
     * Subclasses with non-trivial basis transforms should override for efficiency.
     */
    default void exportPreOrderPartialToStandard(double[] src, int srcOff, double[] dst, int dstOff) {
        final int K = getStateCount();
        final double[] tmpIn  = new double[K];
        final double[] tmpOut = new double[K];
        System.arraycopy(src, srcOff, tmpIn, 0, K);
        exportPreOrderPartialToStandard(tmpIn, tmpOut);
        System.arraycopy(tmpOut, 0, dst, dstOff, K);
    }

    /**
     * Propagate a pre-order message from the TOP of the child's branch to the BOTTOM
     * of the child's branch (i.e. to the child node).
     */
    void propagateToBranchBottom(int childNodeNumber,
                                 double branchLength,
                                 double[] childBranchTopPreOrder,
                                 double[] outChildNodePreOrder);

    /**
     * Export one internal pre-order slice to the representation's external/reporting coordinates.
     */
    void exportPreOrderPartial(double[] preOrderPartial, double[] outPartial);
}
