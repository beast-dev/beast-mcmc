package dr.evomodel.continuous.ou.kernel;

/**
 * Marker for branch-local state prepared by an {@link OUBranchKernel}.
 */
public interface OUPreparedBranch {

    void markInvalid();
}
