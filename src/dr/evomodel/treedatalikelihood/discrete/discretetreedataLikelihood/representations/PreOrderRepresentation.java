package dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.representations;

import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.DiscretePreOrderDelegate;

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

    default void initializeRootPartial(double[] rootFrequencies, double[] outRootPreOrder, int outOffset) {
        final int K = getStateCount();
        for (int s = 0; s < K; s++) {
            outRootPreOrder[outOffset + s] = rootFrequencies[s];
        }
    }


    /**
     * Combine parent pre-order with sibling post-order, reading the parent slice
     * from {@code parentNodePreOrder[parentOff..parentOff+K]}.
     */
    default void combineParentAndSibling(double[] parentNodePreOrder, int parentOff,
                                         double[] siblingBranchTopPostOrder,
                                         double[] outChildBranchTopPreOrder) {
        combineParentAndSibling(parentNodePreOrder, parentOff,
                siblingBranchTopPostOrder, outChildBranchTopPreOrder, 0);
    }

    default void combineParentAndSibling(double[] parentNodePreOrder, int parentOff,
                                         double[] siblingBranchTopPostOrder,
                                         double[] outChildBranchTopPreOrder,
                                         int outOff) {
        combineParentAndSibling(parentNodePreOrder, parentOff,
                siblingBranchTopPostOrder, 0,
                outChildBranchTopPreOrder, outOff);
    }

    default void combineParentAndSibling(double[] parentNodePreOrder, int parentOff,
                                         double[] siblingBranchTopPostOrder, int siblingOff,
                                         double[] outChildBranchTopPreOrder,
                                         int outOff) {
        final int K = getStateCount();
        for (int s = 0; s < K; s++) {
            outChildBranchTopPreOrder[outOff + s] =
                    parentNodePreOrder[parentOff + s] * siblingBranchTopPostOrder[siblingOff + s];
        }
    }

    /**
     * Whether the internal pre-order partials are already expressed in the
     * standard data-type basis.
     */
    default boolean storesPartialsInStandardBasis() {
        return true;
    }

    /**
     * Whether traversal should cache only branch-top pre-order messages.
     * Implementations may still compute branch-bottom messages transiently.
     */
    default boolean cacheOnlyBranchTopPreOrder() {
        return false;
    }

    /**
     * Convert the element-wise product of two standard-basis partials to the
     * internal pre-order representation.
     */
    default void importPreOrderProductFromStandard(double[] leftStandard, int leftOff,
                                                   double[] rightStandard,
                                                   double[] outPreOrderPartial) {
        importPreOrderProductFromStandard(leftStandard, leftOff, rightStandard, outPreOrderPartial, 0);
    }

    default void importPreOrderProductFromStandard(double[] leftStandard, int leftOff,
                                                   double[] rightStandard,
                                                   double[] outPreOrderPartial,
                                                   int outOff) {
        importPreOrderProductFromStandard(leftStandard, leftOff,
                rightStandard, 0, outPreOrderPartial, outOff);
    }

    default void importPreOrderProductFromStandard(double[] leftStandard, int leftOff,
                                                   double[] rightStandard, int rightOff,
                                                   double[] outPreOrderPartial,
                                                   int outOff) {
        final int K = getStateCount();
        for (int s = 0; s < K; s++) {
            outPreOrderPartial[outOff + s] =
                    leftStandard[leftOff + s] * rightStandard[rightOff + s];
        }
    }

    /**
     * Offset-aware export: reads src[srcOff..srcOff+K] and writes dst[dstOff..dstOff+K].
     */
    default void exportPreOrderPartialToStandard(double[] src, int srcOff, double[] dst, int dstOff) {
        final int K = getStateCount();
        for (int s = 0; s < K; s++) {
            dst[dstOff + s] = src[srcOff + s];
        }
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


    void propagateToBranchBottom(int childNodeNumber,
                                 double branchLength,
                                 double[] childBranchTopPreOrder,
                                 int childBranchTopOffset,
                                 double[] out, int outOff);

    /**
     * Export one internal pre-order slice to the representation's external/reporting coordinates.
     */
    void exportPreOrderPartial(double[] preOrderPartial, double[] outPartial);

    default void exportPreOrderPartial(double[] preOrderPartial, int preOrderOffset,
                                       double[] outPartial, int outOffset) {
        final int K = getStateCount();
        for (int s = 0; s < K; s++) {
            outPartial[outOffset + s] = preOrderPartial[preOrderOffset + s];
        }
    }
}
