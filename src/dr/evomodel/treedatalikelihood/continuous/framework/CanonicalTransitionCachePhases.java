package dr.evomodel.treedatalikelihood.continuous.framework;

/**
 * Shared transition-cache diagnostic phase names for canonical tree traversal
 * and gradient preparation.
 */
public final class CanonicalTransitionCachePhases {

    public static final String POSTORDER = "postorder";
    public static final String PREORDER = "preorder";
    public static final String GRADIENT_PREP = "gradientPrep";
    public static final String BRANCH_LENGTH_GRADIENT = "branchLengthGradient";

    private CanonicalTransitionCachePhases() {
        // no instances
    }
}
