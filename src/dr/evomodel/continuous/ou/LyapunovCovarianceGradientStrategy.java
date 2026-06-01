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
    public void accumulateFlat(final double dt,
                               final int dimension,
                               final double[] selection,
                               final double[] diffusion,
                               final double[] covarianceAdjoint,
                               final boolean transposeAdjoint,
                               final double[] gradientAccumulator) {
        final OUCovarianceGradientWorkspace workspace = workspace(dimension);
        final double[] square = workspace.squareMatrices;
        final int gSym = workspace.squareOffset(0);
        OUCovarianceGradientMath.fillSymmetricFromFlat(
                covarianceAdjoint, transposeAdjoint, dimension, square, gSym);
        accumulateSymmetric(
                dt,
                dimension,
                selection,
                0,
                diffusion,
                0,
                square,
                gSym,
                workspace,
                gradientAccumulator);
    }

    private void accumulateSymmetric(final double dt,
                                     final int dimension,
                                     final double[] selection,
                                     final int selectionOffset,
                                     final double[] diffusion,
                                     final int diffusionOffset,
                                     final double[] gSym,
                                     final int gSymOffset,
                                     final OUCovarianceGradientWorkspace workspace,
                                     final double[] gradientAccumulator) {
        final double[] square = workspace.squareMatrices;
        final double[] block = workspace.blockMatrices;
        final int fRemaining = workspace.squareOffset(1);
        final int fRemainingT = workspace.squareOffset(2);
        final int psi = workspace.squareOffset(3);
        final int tempDxD = workspace.squareOffset(4);
        final int vS = workspace.squareOffset(5);
        final int expScratch = workspace.squareOffset(6);
        final int vanLoan = workspace.blockOffset(0);
        final int vanLoanExp = workspace.blockOffset(1);

        for (int idx = 0; idx < OUCovarianceGradientMath.GL5_NODES.length; ++idx) {
            final double s = 0.5 * dt * (OUCovarianceGradientMath.GL5_NODES[idx] + 1.0);
            final double scaledWeight = 0.5 * dt * OUCovarianceGradientMath.GL5_WEIGHTS[idx];

            OUCovarianceGradientMath.buildExpmMinusAs(
                    dt - s,
                    selection,
                    selectionOffset,
                    dimension,
                    square,
                    fRemaining,
                    square,
                    expScratch);
            MatrixExponentialUtils.transposeFlat(square, fRemaining, square, fRemainingT, dimension);
            MatrixExponentialUtils.multiplyFlat(square, fRemainingT, gSym, gSymOffset, square, tempDxD, dimension);
            MatrixExponentialUtils.multiplyFlat(square, tempDxD, square, fRemaining, square, psi, dimension);
            OUCovarianceGradientMath.buildVanLoanCovariance(
                    s,
                    selection,
                    selectionOffset,
                    diffusion,
                    diffusionOffset,
                    dimension,
                    square,
                    vS,
                    block,
                    vanLoan,
                    block,
                    vanLoanExp);

            MatrixExponentialUtils.multiplyFlat(square, vS, square, psi, square, tempDxD, dimension);
            for (int i = 0; i < dimension; ++i) {
                final int rowOffset = i * dimension;
                for (int j = 0; j < dimension; ++j) {
                    gradientAccumulator[rowOffset + j] +=
                            -2.0 * scaledWeight * square[tempDxD + rowOffset + j];
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
