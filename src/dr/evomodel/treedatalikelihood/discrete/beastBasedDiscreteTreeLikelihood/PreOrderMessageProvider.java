package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood;

import dr.evomodel.treedatalikelihood.preorder.DiscretePreOrderType;

/**
 * Narrow interface exposing exported pre-order branch-side messages.
 */
/*
 * @author Filippo Monti
 */
public interface PreOrderMessageProvider {

    int getStateCount();

    int getPatternCount();

    int getCategoryCount();

    /**
     * Get the exported pre-order message at the TOP of the branch leading into childNodeNumber,
     * for a single category/pattern slice.
     */
    default void getPreOrderBranchTopInto(int childNodeNumber, int category, int pattern, double[] outPartial) {
        throw new RuntimeException("Not implemented");
    }

    /**
     * Get the exported pre-order message at the BOTTOM of the branch leading into childNodeNumber,
     * for a single category/pattern slice.
     */
    default void getPreOrderBranchBottomInto(int childNodeNumber, int category, int pattern, double[] outPartial) {
        throw new RuntimeException("Not implemented");
    }

    void getPreorderPartials(int nodeNumber, DiscretePreOrderType type, double[] out);
}
