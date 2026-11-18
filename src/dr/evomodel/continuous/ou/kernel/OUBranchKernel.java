package dr.evomodel.continuous.ou.kernel;

import dr.inference.model.MatrixParameterInterface;

/**
 * Branch-local OU algebra independent of tree traversal and tip observation
 * storage. Implementations fill caller-owned prepared state and workspaces.
 */
public interface OUBranchKernel {

    int getDimension();

    OUPreparedBranch createPreparedBranch();

    OUBranchKernelWorkspace createWorkspace();

    default <T> T getCapability(final Class<T> capabilityType) {
        return capabilityType.isInstance(this) ? capabilityType.cast(this) : null;
    }

    void prepareBranch(double effectiveBranchLength,
                       MatrixParameterInterface diffusionCovariance,
                       double[] stationaryMean,
                       OUPreparedBranch preparedBranch,
                       OUBranchKernelWorkspace workspace);
}
