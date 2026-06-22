package dr.evomodel.continuous.ou.kernel;

/**
 * Optional capability for OU selection parametrizations that can expose a
 * specialized branch kernel without depending on tree-likelihood storage.
 */
public interface OUBranchKernelCapability {

    OUBranchKernel createBranchKernel();
}
