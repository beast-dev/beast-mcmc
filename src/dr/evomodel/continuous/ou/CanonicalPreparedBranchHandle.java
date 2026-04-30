package dr.evomodel.continuous.ou;

/**
 * Opaque, reusable branch-local data prepared by a canonical OU backend.
 */
public interface CanonicalPreparedBranchHandle {

    void invalidateCovariance();
}
