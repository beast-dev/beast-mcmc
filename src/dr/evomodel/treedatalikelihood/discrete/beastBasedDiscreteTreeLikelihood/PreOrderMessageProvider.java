package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood;

import beagle.BeaglePreorderType;

/**
 * Narrow interface exposing pre-order branch-side messages in the STANDARD basis.
 */
/*
 * @author Filippo Monti
 */
public interface PreOrderMessageProvider {

    int getStateCount();

    int getPatternCount();

    int getCategoryCount();

    /**
     * Get the pre-order message at the TOP of the branch leading into childNodeNumber,
     * in the standard basis, for a single category/pattern slice.
     */
    void getPreOrderBranchTopInto(int childNodeNumber, int category, int pattern, double[] outStandardPartial);

    /**
     * Get the pre-order message at the BOTTOM of the branch leading into childNodeNumber,
     * in the standard basis, for a single category/pattern slice.
     */
    void getPreOrderBranchBottomInto(int childNodeNumber, int category, int pattern, double[] outStandardPartial);

    void getPreorderPartials(int nodeNumber, BeaglePreorderType type, double[] out);
}
