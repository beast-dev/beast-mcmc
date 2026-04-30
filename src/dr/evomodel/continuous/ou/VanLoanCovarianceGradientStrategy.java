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
    public void accumulate(final double dt,
                           final int dimension,
                           final double[][] selection,
                           final double[][] diffusion,
                           final double[][] covarianceAdjoint,
                           final double[] gradientAccumulator) {
        final OUCovarianceGradientWorkspace workspace = workspace(dimension);

        final double[][] vanLoanMatrix = workspace.blockMatrices[0];
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                vanLoanMatrix[i][j] = -dt * selection[i][j];
                vanLoanMatrix[i][j + dimension] = dt * diffusion[i][j];
                vanLoanMatrix[i + dimension][j] = 0.0;
                vanLoanMatrix[i + dimension][j + dimension] = dt * selection[j][i];
            }
        }

        final double[][] vanLoanExp = workspace.blockMatrices[1];
        MatrixExponentialUtils.expm(vanLoanMatrix, vanLoanExp);

        final double[][] gvF = workspace.squareMatrices[0];
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += covarianceAdjoint[i][k] * vanLoanExp[k][j];
                }
                gvF[i][j] = sum;
            }
        }

        accumulateFromVanLoanBlocks(dt, dimension, vanLoanMatrix, vanLoanExp, gvF, workspace, gradientAccumulator);
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

        final double[][] vanLoanMatrix = workspace.blockMatrices[0];
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                vanLoanMatrix[i][j] = -dt * selection[i][j];
                vanLoanMatrix[i][j + dimension] = dt * diffusion[i][j];
                vanLoanMatrix[i + dimension][j] = 0.0;
                vanLoanMatrix[i + dimension][j + dimension] = dt * selection[j][i];
            }
        }

        final double[][] vanLoanExp = workspace.blockMatrices[1];
        MatrixExponentialUtils.expm(vanLoanMatrix, vanLoanExp);

        final double[][] gvF = workspace.squareMatrices[0];
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += OUCovarianceGradientMath.flatSquareValue(
                            covarianceAdjoint, i, k, dimension, transposeAdjoint) * vanLoanExp[k][j];
                }
                gvF[i][j] = sum;
            }
        }

        accumulateFromVanLoanBlocks(dt, dimension, vanLoanMatrix, vanLoanExp, gvF, workspace, gradientAccumulator);
    }

    private void accumulateFromVanLoanBlocks(final double dt,
                                             final int dimension,
                                             final double[][] vanLoanMatrix,
                                             final double[][] vanLoanExp,
                                             final double[][] gvF,
                                             final OUCovarianceGradientWorkspace workspace,
                                             final double[] gradientAccumulator) {
        final double[][] v = workspace.squareMatrices[1];
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += vanLoanExp[i][k + dimension] * vanLoanExp[j][k];
                }
                v[i][j] = sum;
            }
        }

        final double[][] upstreamBlock = workspace.blockMatrices[2];
        OUCovarianceGradientMath.clearSquare(upstreamBlock);
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                upstreamBlock[i][j + dimension] = gvF[i][j];
                double vgvf = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    vgvf += v[i][k] * gvF[k][j];
                }
                upstreamBlock[i + dimension][j + dimension] = -vgvf;
            }
        }

        final double[][] gradM = workspace.blockMatrices[3];
        MatrixExponentialUtils.adjointExp(vanLoanMatrix, upstreamBlock, gradM);

        for (int k = 0; k < dimension; ++k) {
            for (int l = 0; l < dimension; ++l) {
                gradientAccumulator[k * dimension + l] +=
                        -dt * gradM[k][l] + dt * gradM[dimension + l][dimension + k];
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
