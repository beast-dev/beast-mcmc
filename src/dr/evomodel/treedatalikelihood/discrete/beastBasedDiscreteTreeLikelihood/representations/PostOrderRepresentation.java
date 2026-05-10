package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.representations;

import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;

/**
 * Defines the internal algebra of post-order partials used by
 * {@link DiscreteDataLikelihoodDelegate}.
 *
 * Implementations may store partials in any coordinate system, as long as they provide:
 *
 * 1. tip initialization,
 * 2. branch propagation to the branch top,
 * 3. child combination at internal nodes,
 * 4. root evaluation,
 * 5. export for caches/debugging.
 *
 * All methods are allocation-free with caller-owned buffers.
 */
/*
 * @author Filippo Monti
 */
public interface PostOrderRepresentation {

    String getName();

    int getStateCount();

    void markDirty();

    void storeState();

    void restoreState();

    void updateForLikelihood();

    /**
     * Initialize one tip partial from a tip vector expressed in the standard basis.
     *
     * @param standardTip tip vector in standard basis
     * @param outPartial  output partial in the internal representation
     */
    void initializeTipPartial(double[] standardTip, double[] outPartial);

    /**
     * Propagate a child partial from the child node to the top of its incident branch.
     *
     * @param nodeNumber          child node number identifying the branch
     * @param branchLength        effective branch length
     * @param childPartial        child partial in internal representation
     * @param outBranchTopPartial propagated partial at branch top in internal representation
     */
    void propagateToBranchTop(int nodeNumber,
                              double branchLength,
                              double[] childPartial,
                              double[] outBranchTopPartial);

    /**
     * Combine the two child branch-top partials into the parent node partial.
     *
     * @param leftBranchTopPartial  left child contribution at parent side
     * @param rightBranchTopPartial right child contribution at parent side
     * @param outParentPartial      parent partial in internal representation
     */
    void combineBranchTopPartials(double[] leftBranchTopPartial,
                                  double[] rightBranchTopPartial,
                                  double[] outParentPartial);

    /**
     * Combine child branch-top partials and optionally expose their standard-basis
     * forms to the caller. Implementations that already compute standard-basis
     * intermediates can override this to avoid duplicate basis transforms.
     */
    default void combineBranchTopPartials(double[] leftBranchTopPartial,
                                          double[] rightBranchTopPartial,
                                          double[] outParentPartial,
                                          double[] outLeftStandard,
                                          double[] outRightStandard) {
        combineBranchTopPartials(leftBranchTopPartial, rightBranchTopPartial, outParentPartial);
        if (outLeftStandard != null) {
            exportPostOrderPartialToStandard(leftBranchTopPartial, outLeftStandard);
        }
        if (outRightStandard != null) {
            exportPostOrderPartialToStandard(rightBranchTopPartial, outRightStandard);
        }
    }

    /**
     * Evaluate the root contribution for one root partial slice.
     *
     * @param rootFrequencies root frequencies in standard basis
     * @param rootPartial     root partial in internal representation
     * @return root contribution
     */
    double rootContribution(double[] rootFrequencies, double[] rootPartial);

    /**
     * Export one partial slice from the internal representation.
     *
     * @param partial    input partial in internal representation
     * @param outPartial output partial in the representation's external/reporting coordinates
     */
    void exportPostOrderPartial(double[] partial, double[] outPartial);

    /**
     * Whether the internal post-order partials are already expressed in the
     * standard data-type basis. Implementations with spectral/internal bases can
     * override this so traversal code can maintain standard-basis side caches.
     */
    default boolean storesPartialsInStandardBasis() {
        return true;
    }

    /**
     * Export one partial slice to the standard data-type basis. For standard-basis
     * implementations this is identical to the usual export operation.
     */
    default void exportPostOrderPartialToStandard(double[] partial, double[] outPartial) {
        exportPostOrderPartial(partial, outPartial);
    }

    /**
     * Whether it is valid to apply scaling directly to the internal representation.
     * This is true for the standard and spectral implementations here.
     */
    default boolean supportsInternalScaling() {
        return true;
    }
}
