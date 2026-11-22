package dr.evomodel.continuous.ou;

final class VanLoanCovarianceGradientStrategy implements OUCovarianceGradientStrategy {

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
        final double[] block = workspace.blockMatrices;
        final double[] square = workspace.squareMatrices;
        final int twoDim = 2 * dimension;
        final int vanLoanMatrix = workspace.blockOffset(0);
        final int vanLoanExp = workspace.blockOffset(1);

        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            final int topRow = vanLoanMatrix + i * twoDim;
            final int bottomRow = vanLoanMatrix + (i + dimension) * twoDim;
            for (int j = 0; j < dimension; ++j) {
                block[topRow + j] = -dt * selection[rowOffset + j];
                block[topRow + j + dimension] = dt * diffusion[rowOffset + j];
                block[bottomRow + j] = 0.0;
                block[bottomRow + j + dimension] = dt * selection[j * dimension + i];
            }
        }

        MatrixExponentialUtils.expmFlat(block, vanLoanMatrix, block, vanLoanExp, twoDim);

        final int gvF = workspace.squareOffset(0);
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += OUCovarianceGradientMath.flatSquareValue(
                            covarianceAdjoint, i, k, dimension, transposeAdjoint)
                            * block[vanLoanExp + k * twoDim + j];
                }
                square[gvF + rowOffset + j] = sum;
            }
        }

        accumulateFromVanLoanBlocks(
                dt,
                dimension,
                workspace,
                vanLoanMatrix,
                vanLoanExp,
                gvF,
                gradientAccumulator);
    }

    private void accumulateFromVanLoanBlocks(final double dt,
                                             final int dimension,
                                             final OUCovarianceGradientWorkspace workspace,
                                             final int vanLoanMatrix,
                                             final int vanLoanExp,
                                             final int gvF,
                                             final double[] gradientAccumulator) {
        final double[] square = workspace.squareMatrices;
        final double[] block = workspace.blockMatrices;
        final int twoDim = 2 * dimension;
        final int v = workspace.squareOffset(1);
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += block[vanLoanExp + i * twoDim + k + dimension]
                            * block[vanLoanExp + j * twoDim + k];
                }
                square[v + rowOffset + j] = sum;
            }
        }

        final int upstreamBlock = workspace.blockOffset(2);
        OUCovarianceGradientMath.clearMatrix(block, upstreamBlock, workspace.blockStride);
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            final int topRow = upstreamBlock + i * twoDim;
            final int bottomRow = upstreamBlock + (i + dimension) * twoDim;
            for (int j = 0; j < dimension; ++j) {
                block[topRow + j + dimension] = square[gvF + rowOffset + j];
                double vgvf = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    vgvf += square[v + rowOffset + k] * square[gvF + k * dimension + j];
                }
                block[bottomRow + j + dimension] = -vgvf;
            }
        }

        final int gradM = workspace.blockOffset(3);
        MatrixExponentialUtils.adjointExpFlat(
                block, vanLoanMatrix, block, upstreamBlock, block, gradM, twoDim);

        for (int k = 0; k < dimension; ++k) {
            final int rowOffset = k * dimension;
            for (int l = 0; l < dimension; ++l) {
                gradientAccumulator[rowOffset + l] +=
                        -dt * block[gradM + k * twoDim + l]
                                + dt * block[gradM + (dimension + l) * twoDim + dimension + k];
            }
        }
    }

    private OUCovarianceGradientWorkspace workspace(final int dimension) {
        OUCovarianceGradientWorkspace workspace = workspaceLocal.get();
        if (workspace == null || workspace.dimension != dimension) {
            workspace = new OUCovarianceGradientWorkspace(dimension, 2, 4, false);
            workspaceLocal.set(workspace);
        }
        return workspace;
    }
}
