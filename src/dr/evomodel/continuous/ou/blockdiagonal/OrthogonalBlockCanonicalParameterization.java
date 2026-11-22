package dr.evomodel.continuous.ou.blockdiagonal;

/**
 * Canonical OU capability exposed by orthogonal block-diagonal selection parametrizations.
 */
public interface OrthogonalBlockCanonicalParameterization
        extends BlockDiagonalNativeCanonicalParameterization {

    @Override
    OrthogonalBlockPreparedBranchBasis createPreparedBranchBasis();
}
