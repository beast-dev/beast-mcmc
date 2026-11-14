package dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;

/**
 * Handle for a block-diagonal prepared branch basis, implementing
 * {@link CanonicalPreparedBranchHandle} for use by callers that only need the
 * capability interface. {@link #getBasis()} exposes the wrapped, mutable
 * {@link BlockDiagonalPreparedBranchBasis} directly; it is intended for
 * same-package production code and white-box tests, not as a general
 * external accessor.
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
