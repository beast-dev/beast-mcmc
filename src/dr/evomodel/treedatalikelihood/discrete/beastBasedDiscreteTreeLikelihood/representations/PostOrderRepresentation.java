package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.representations;

import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;

/**
 * Defines the internal algebra of post-order partials used by
 * {@link DiscreteDataLikelihoodDelegate}.
 *
 * Implementations may store partials in the standard basis, a spectral basis,
 * or any other representation, as long as they provide:
 *
 * 1. tip initialization,
 * 2. branch propagation to the branch top,
 * 3. child combination at internal nodes,
 * 4. root evaluation,
 * 5. export to standard basis for caches/debugging.
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
     * Evaluate the root contribution for one root partial slice.
     *
     * @param rootFrequencies root frequencies in standard basis
     * @param rootPartial     root partial in internal representation
     * @return root contribution
     */
    double rootContribution(double[] rootFrequencies, double[] rootPartial);

    /**
     * Convert one partial slice from the internal representation to the standard basis.
     *
     * @param partial             input partial in internal representation
     * @param outStandardPartial  output partial in standard basis
     */
    void toStandard(double[] partial, double[] outStandardPartial);

    /**
     * Whether it is valid to apply scaling directly to the internal representation.
     * This is true for the standard and spectral implementations here.
     */
    default boolean supportsInternalScaling() {
        return true;
    }
}
