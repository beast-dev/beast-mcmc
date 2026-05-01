package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;

/**
 * Opaque public handle for an orthogonal block-diagonal prepared branch basis.
 */
public final class OrthogonalBlockPreparedBranchHandle implements CanonicalPreparedBranchHandle {

    private final OrthogonalBlockPreparedBranchBasis basis;

    OrthogonalBlockPreparedBranchHandle(final OrthogonalBlockPreparedBranchBasis basis) {
        this.basis = basis;
    }

    public OrthogonalBlockPreparedBranchBasis getBasis() {
        return basis;
    }

    @Override
    public void invalidateCovariance() {
        basis.invalidateCovariance();
    }
}
