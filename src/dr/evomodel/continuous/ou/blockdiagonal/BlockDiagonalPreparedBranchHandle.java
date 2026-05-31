package dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;

/**
 * Opaque public handle for a block-diagonal prepared branch basis.
 */
public class BlockDiagonalPreparedBranchHandle implements CanonicalPreparedBranchHandle {

    private final BlockDiagonalPreparedBranchBasis basis;

    BlockDiagonalPreparedBranchHandle(final BlockDiagonalPreparedBranchBasis basis) {
        this.basis = basis;
    }

    public BlockDiagonalPreparedBranchBasis getBasis() {
        return basis;
    }

    @Override
    public void invalidateCovariance() {
        basis.invalidateCovariance();
    }
}
