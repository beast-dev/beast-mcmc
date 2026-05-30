package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

/**
 * Opaque public handle for an orthogonal block-diagonal prepared branch basis.
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
