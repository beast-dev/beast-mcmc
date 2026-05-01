package dr.evomodel.continuous.ou.canonical;

/**
 * Opaque, reusable branch-local data prepared by a canonical OU backend.
 */
public interface CanonicalPreparedBranchHandle {

    void invalidateCovariance();
}
