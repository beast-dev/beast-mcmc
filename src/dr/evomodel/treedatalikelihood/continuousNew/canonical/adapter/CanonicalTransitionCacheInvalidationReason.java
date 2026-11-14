package dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter;

enum CanonicalTransitionCacheInvalidationReason {
    DIFFUSION_CHANGED,
    STATIONARY_MEAN_CHANGED,
    SELECTION_CHANGED,
    BRANCH_LENGTH_CHANGED,
    MODEL_CHANGED,
    RESTORE_STATE,
    EXPLICIT
}
