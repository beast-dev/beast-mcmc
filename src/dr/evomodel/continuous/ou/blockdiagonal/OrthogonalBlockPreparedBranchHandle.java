package dr.evomodel.continuous.ou.blockdiagonal;

/**
 * Handle for an orthogonal block-diagonal prepared branch basis. See
 * {@link BlockDiagonalPreparedBranchHandle} for the {@link #getBasis()}
 * access contract.
 */
public final class OrthogonalBlockPreparedBranchHandle extends BlockDiagonalPreparedBranchHandle {

    OrthogonalBlockPreparedBranchHandle(final OrthogonalBlockPreparedBranchBasis basis) {
        super(basis);
    }

    @Override
    public OrthogonalBlockPreparedBranchBasis getBasis() {
        return (OrthogonalBlockPreparedBranchBasis) super.getBasis();
    }
}
