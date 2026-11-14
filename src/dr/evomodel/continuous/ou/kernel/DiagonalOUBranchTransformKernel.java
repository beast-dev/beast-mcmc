package dr.evomodel.continuous.ou.kernel;

/**
 * Capability for OU kernels that expose a diagonal actualization branch transform.
 *
 * <p>The arrays returned here are owned by the supplied prepared branch and must
 * not be retained beyond that branch's cache lifetime.</p>
 */
public interface DiagonalOUBranchTransformKernel extends OUBranchKernel {

    double[] actualizationDiagonal(OUPreparedBranch preparedBranch);

    double[] displacement(OUPreparedBranch preparedBranch);

    double[] precision(OUPreparedBranch preparedBranch);

    double logNormalizer(OUPreparedBranch preparedBranch);

    boolean isCurrent(OUPreparedBranch preparedBranch, double branchLength);
}
