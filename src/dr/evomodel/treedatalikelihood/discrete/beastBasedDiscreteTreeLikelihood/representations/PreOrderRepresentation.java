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
     * Combine parent pre-order with sibling post-order, reading the parent slice
     * from {@code parentNodePreOrder[parentOff..parentOff+K]}.
     * Default: copies the slice into a temp buffer, then calls the zero-offset variant.
     */
    default void combineParentAndSibling(double[] parentNodePreOrder, int parentOff,
                                         double[] siblingBranchTopPostOrder,
                                         double[] outChildBranchTopPreOrder) {
        final int K = getStateCount();
        final double[] tmp = new double[K];
        System.arraycopy(parentNodePreOrder, parentOff, tmp, 0, K);
        combineParentAndSibling(tmp, siblingBranchTopPostOrder, outChildBranchTopPreOrder);
    }

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
     * Convert the element-wise product of two standard-basis partials to the
     * internal pre-order representation.
     */
    default void importPreOrderProductFromStandard(double[] leftStandard, int leftOff,
                                                   double[] rightStandard,
                                                   double[] outPreOrderPartial) {
        final int K = getStateCount();
        final double[] tmp = new double[K];
        for (int s = 0; s < K; s++) {
            tmp[s] = leftStandard[leftOff + s] * rightStandard[s];
        }
        importPreOrderPartialFromStandard(tmp, outPreOrderPartial);
    }

    /**
     * Export one internal pre-order partial slice to the standard data-type basis.
     */
    default void exportPreOrderPartialToStandard(double[] preOrderPartial, double[] outStandardPartial) {
        exportPreOrderPartial(preOrderPartial, outStandardPartial);
    }

    /**
     * Offset-aware export: reads src[srcOff..srcOff+K] and writes dst[dstOff..dstOff+K].
     * Default: copies slices into temp buffers, then calls the zero-offset variant.
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
     * Normalize one internal pre-order slice in place and export the normalized
     * value to the standard data-type basis.
     *
     * @return the log scale added by normalization.
     */
    default double normalizeAndExportPreOrderPartialToStandard(double[] src, int srcOff,
                                                               double[] dst, int dstOff,
                                                               double scalingFloor,
                                                               double scalingCeiling) {
        final int K = getStateCount();
        double max = 0.0;
        for (int s = 0; s < K; s++) {
            max = Math.max(max, Math.abs(src[srcOff + s]));
        }

        if (max == 0.0) {
            final double uniform = 1.0 / K;
            for (int s = 0; s < K; s++) {
                src[srcOff + s] = uniform;
            }
            exportPreOrderPartialToStandard(src, srcOff, dst, dstOff);
            return Double.NEGATIVE_INFINITY;
        }

        if (max < scalingFloor || max > scalingCeiling) {
            for (int s = 0; s < K; s++) {
                src[srcOff + s] /= max;
            }
            exportPreOrderPartialToStandard(src, srcOff, dst, dstOff);
            return Math.log(max);
        }

        exportPreOrderPartialToStandard(src, srcOff, dst, dstOff);
        return 0.0;
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
     * Offset-aware variant: writes the result to {@code out[outOff..outOff+K]}.
     * Default: calls the zero-offset variant into a temp buffer, then copies.
     */
    default void propagateToBranchBottom(int childNodeNumber,
                                         double branchLength,
                                         double[] childBranchTopPreOrder,
                                         double[] out, int outOff) {
        final int K = getStateCount();
        final double[] tmp = new double[K];
        propagateToBranchBottom(childNodeNumber, branchLength, childBranchTopPreOrder, tmp);
        System.arraycopy(tmp, 0, out, outOff, K);
    }

    /**
     * Export one internal pre-order slice to the representation's external/reporting coordinates.
     */
    void exportPreOrderPartial(double[] preOrderPartial, double[] outPartial);
}
