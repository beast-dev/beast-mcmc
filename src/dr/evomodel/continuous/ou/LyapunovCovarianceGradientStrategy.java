package dr.evomodel.continuous.ou;

final class LyapunovCovarianceGradientStrategy implements OUCovarianceGradientStrategy {

    private final ThreadLocal<OUCovarianceGradientWorkspace> workspaceLocal =
            new ThreadLocal<OUCovarianceGradientWorkspace>() {
                @Override
                protected OUCovarianceGradientWorkspace initialValue() {
                    return null;
                }
            };

    @Override
    public void accumulate(final double dt,
                           final int dimension,
                           final double[][] selection,
                           final double[][] diffusion,
                           final double[][] covarianceAdjoint,
                           final double[] gradientAccumulator) {
        final OUCovarianceGradientWorkspace workspace = workspace(dimension);
        final double[][] gSym = workspace.squareMatrices[0];
        OUCovarianceGradientMath.fillSymmetricFromMatrix(covarianceAdjoint, dimension, gSym);
        accumulateSymmetric(dt, dimension, selection, diffusion, gSym, workspace, gradientAccumulator);
    }

    @Override
    public void accumulateFlat(final double dt,
                               final int dimension,
                               final double[][] selection,
                               final double[][] diffusion,
                               final double[] covarianceAdjoint,
                               final boolean transposeAdjoint,
                               final double[] gradientAccumulator) {
        final OUCovarianceGradientWorkspace workspace = workspace(dimension);
        final double[][] gSym = workspace.squareMatrices[0];
        OUCovarianceGradientMath.fillSymmetricFromFlat(covarianceAdjoint, transposeAdjoint, dimension, gSym);
        accumulateSymmetric(dt, dimension, selection, diffusion, gSym, workspace, gradientAccumulator);
    }

    private void accumulateSymmetric(final double dt,
                                     final int dimension,
                                     final double[][] selection,
                                     final double[][] diffusion,
                                     final double[][] gSym,
                                     final OUCovarianceGradientWorkspace workspace,
                                     final double[] gradientAccumulator) {
        final double[][] fRemaining = workspace.squareMatrices[1];
        final double[][] fRemainingT = workspace.squareMatrices[2];
        final double[][] psi = workspace.squareMatrices[3];
        final double[][] tempDxD = workspace.squareMatrices[4];
        final double[][] vS = workspace.squareMatrices[5];
        final double[][] expScratch = workspace.squareMatrices[6];

        for (int idx = 0; idx < OUCovarianceGradientMath.GL5_NODES.length; ++idx) {
            final double s = 0.5 * dt * (OUCovarianceGradientMath.GL5_NODES[idx] + 1.0);
            final double scaledWeight = 0.5 * dt * OUCovarianceGradientMath.GL5_WEIGHTS[idx];

            OUCovarianceGradientMath.buildExpmMinusAs(
                    dt - s, selection, dimension, fRemaining, expScratch);
            MatrixExponentialUtils.transpose(fRemaining, fRemainingT);
            MatrixExponentialUtils.multiply(fRemainingT, gSym, tempDxD);
            MatrixExponentialUtils.multiply(tempDxD, fRemaining, psi);
            OUCovarianceGradientMath.buildVanLoanCovariance(
                    s, selection, diffusion, dimension, vS,
                    workspace.blockMatrices[0], workspace.blockMatrices[1]);

            MatrixExponentialUtils.multiply(vS, psi, tempDxD);
            for (int i = 0; i < dimension; ++i) {
                for (int j = 0; j < dimension; ++j) {
                    gradientAccumulator[i * dimension + j] += -2.0 * scaledWeight * tempDxD[i][j];
                }
            }
        }
    }

    private OUCovarianceGradientWorkspace workspace(final int dimension) {
        OUCovarianceGradientWorkspace workspace = workspaceLocal.get();
        if (workspace == null || workspace.dimension != dimension) {
            workspace = new OUCovarianceGradientWorkspace(dimension, 7, 2, false);
            workspaceLocal.set(workspace);
        }
        return workspace;
    }
}
