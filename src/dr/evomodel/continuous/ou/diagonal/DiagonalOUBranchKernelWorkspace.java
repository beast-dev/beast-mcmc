package dr.evomodel.continuous.ou.diagonal;

import dr.evomodel.continuous.ou.kernel.OUBranchKernelWorkspace;

public final class DiagonalOUBranchKernelWorkspace implements OUBranchKernelWorkspace {

    final double[] oneMinusActualization;
    final double[] branchCovariance;
    final double[] branchCovarianceScratch;
    final double[] branchCholeskyScratch;
    final double[] branchLowerInverseScratch;

    DiagonalOUBranchKernelWorkspace(final int dimension) {
        final int d2 = dimension * dimension;
        this.oneMinusActualization = new double[dimension];
        this.branchCovariance = new double[d2];
        this.branchCovarianceScratch = new double[d2];
        this.branchCholeskyScratch = new double[d2];
        this.branchLowerInverseScratch = new double[d2];
    }
}
