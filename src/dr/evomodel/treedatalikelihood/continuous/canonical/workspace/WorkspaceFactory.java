package dr.evomodel.treedatalikelihood.continuous.canonical.workspace;

import java.util.EnumSet;

public final class WorkspaceFactory {

    private WorkspaceFactory() {
        // no instances
    }

    public static BranchGradientWorkspace branchGradientWorkspace(final int dimension) {
        return branchGradientWorkspace(
                dimension,
                EnumSet.allOf(WorkspaceCapability.class));
    }

    public static BranchGradientWorkspace branchGradientWorkspace(
            final int dimension,
            final EnumSet<WorkspaceCapability> capabilities) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be >= 1");
        }
        if (capabilities == null || capabilities.isEmpty()) {
            throw new IllegalArgumentException("capabilities must not be empty");
        }
        validateCapabilities(capabilities);
        final BranchGradientWorkspace workspace = new BranchGradientWorkspace(dimension, capabilities);
        workspace.validate(dimension);
        return workspace;
    }

    private static void validateCapabilities(final EnumSet<WorkspaceCapability> capabilities) {
        if (capabilities.contains(WorkspaceCapability.PARTIAL_OBSERVATION)
                && !capabilities.contains(WorkspaceCapability.ADJOINTS)) {
            throw new IllegalArgumentException(
                    "PARTIAL_OBSERVATION requires ADJOINTS workspace capability");
        }
        if (capabilities.contains(WorkspaceCapability.ORTHOGONAL_BLOCK_GRADIENT)
                && !capabilities.contains(WorkspaceCapability.DENSE_GRADIENT)) {
            throw new IllegalArgumentException(
                    "ORTHOGONAL_BLOCK_GRADIENT requires DENSE_GRADIENT for shared gradient accumulators");
        }
    }
}
