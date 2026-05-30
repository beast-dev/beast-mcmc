package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

/**
 * Canonical OU capability exposed by orthogonal block-diagonal selection parametrizations.
 */
public interface OrthogonalBlockCanonicalParameterization
        extends BlockDiagonalNativeCanonicalParameterization {

    @Override
    OrthogonalBlockPreparedBranchBasis createPreparedBranchBasis();
}
