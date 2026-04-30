package dr.evomodel.continuous.ou;

final class StationaryLyapunovCovarianceGradientStrategy implements OUCovarianceGradientStrategy {

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
        final double[][] gV = workspace.squareMatrices[0];
        OUCovarianceGradientMath.fillSymmetricFromMatrix(covarianceAdjoint, dimension, gV);
        accumulateSymmetric(dt, dimension, selection, diffusion, gV, workspace, gradientAccumulator);
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
        final double[][] gV = workspace.squareMatrices[0];
        OUCovarianceGradientMath.fillSymmetricFromFlat(covarianceAdjoint, transposeAdjoint, dimension, gV);
        accumulateSymmetric(dt, dimension, selection, diffusion, gV, workspace, gradientAccumulator);
    }

    private void accumulateSymmetric(final double dt,
                                     final int dimension,
                                     final double[][] selection,
                                     final double[][] diffusion,
                                     final double[][] gV,
                                     final OUCovarianceGradientWorkspace workspace,
                                     final double[] gradientAccumulator) {
        final double[][] sStat = workspace.squareMatrices[1];
        final double[][] f = workspace.squareMatrices[2];
        final double[][] tempDxD = workspace.squareMatrices[3];
        final double[][] gS = workspace.squareMatrices[4];
        final double[][] y = workspace.squareMatrices[5];
        final double[][] yT = workspace.squareMatrices[6];
        final double[][] gAStationary = workspace.squareMatrices[7];
        final double[][] sStatT = workspace.squareMatrices[8];
        final double[][] gFCov = workspace.squareMatrices[9];
        final double[][] gVT = workspace.squareMatrices[10];
        final double[][] tempDxD2 = workspace.squareMatrices[11];
        final double[][] gX = workspace.squareMatrices[12];
        final double[][] fT = workspace.squareMatrices[13];
        final double[][] expScratch = workspace.squareMatrices[14];
        final double[][] minusAdt = workspace.squareMatrices[15];
        final double[][] aT = workspace.squareMatrices[16];

        solveStationaryLyapunovDense(selection, diffusion, dimension, sStat, workspace);
        OUCovarianceGradientMath.buildExpmMinusAs(dt, selection, dimension, f, expScratch);

        MatrixExponentialUtils.transpose(f, fT);
        MatrixExponentialUtils.multiply(fT, gV, tempDxD);
        MatrixExponentialUtils.multiply(tempDxD, f, gS);
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                gS[i][j] = gV[i][j] - gS[i][j];
            }
        }

        solveAdjointStationaryLyapunovDense(selection, gS, dimension, y, aT, workspace);

        MatrixExponentialUtils.transpose(y, yT);
        MatrixExponentialUtils.multiply(y, sStat, gAStationary);
        MatrixExponentialUtils.multiply(yT, sStat, tempDxD);
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                gAStationary[i][j] = -(gAStationary[i][j] + tempDxD[i][j]);
            }
        }

        MatrixExponentialUtils.transpose(sStat, sStatT);
        MatrixExponentialUtils.multiply(gV, f, tempDxD);
        MatrixExponentialUtils.multiply(tempDxD, sStat, gFCov);
        MatrixExponentialUtils.transpose(gV, gVT);
        MatrixExponentialUtils.multiply(gVT, f, tempDxD2);
        MatrixExponentialUtils.multiply(tempDxD2, sStatT, tempDxD);
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                gFCov[i][j] = -(gFCov[i][j] + tempDxD[i][j]);
            }
        }

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                minusAdt[i][j] = -dt * selection[i][j];
            }
        }
        MatrixExponentialUtils.adjointExp(minusAdt, gFCov, gX);

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                gradientAccumulator[i * dimension + j] += gAStationary[j][i] - dt * gX[i][j];
            }
        }
    }

    private static void solveStationaryLyapunovDense(final double[][] selection,
                                                     final double[][] diffusion,
                                                     final int dimension,
                                                     final double[][] out,
                                                     final OUCovarianceGradientWorkspace workspace) {
        final int n = dimension * dimension;
        final double[][] augmented = workspace.lyapunovAugmented;
        if (augmented == null) {
            throw new IllegalStateException("Stationary Lyapunov workspace is not available");
        }

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                final int row = i * dimension + j;
                for (int p = 0; p < dimension; ++p) {
                    for (int r = 0; r < dimension; ++r) {
                        final int col = p * dimension + r;
                        double value = 0.0;
                        if (i == p) {
                            value += selection[j][r];
                        }
                        if (j == r) {
                            value += selection[i][p];
                        }
                        augmented[row][col] = value;
                    }
                }
                augmented[row][n] = diffusion[i][j];
            }
        }

        try {
            solveAugmentedLinearSystemInPlace(augmented);
        } catch (IllegalStateException e) {
            throw new IllegalArgumentException(
                    "Stationary covariance is only defined for stable drift matrices", e);
        }

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i][j] = augmented[i * dimension + j][n];
            }
        }
        MatrixExponentialUtils.symmetrize(out);
    }

    private static void solveAdjointStationaryLyapunovDense(final double[][] selection,
                                                            final double[][] gS,
                                                            final int dimension,
                                                            final double[][] yOut,
                                                            final double[][] aTranspose,
                                                            final OUCovarianceGradientWorkspace workspace) {
        MatrixExponentialUtils.transpose(selection, aTranspose);
        solveStationaryLyapunovDense(aTranspose, gS, dimension, yOut, workspace);
    }

    private static void solveAugmentedLinearSystemInPlace(final double[][] augmented) {
        final int n = augmented.length;
        for (int col = 0; col < n; ++col) {
            int pivot = col;
            double maxAbs = Math.abs(augmented[col][col]);
            for (int row = col + 1; row < n; ++row) {
                final double abs = Math.abs(augmented[row][col]);
                if (abs > maxAbs) {
                    maxAbs = abs;
                    pivot = row;
                }
            }

            if (maxAbs <= 1e-14) {
                throw new IllegalStateException("Singular system");
            }

            if (pivot != col) {
                final double[] tmp = augmented[col];
                augmented[col] = augmented[pivot];
                augmented[pivot] = tmp;
            }

            final double diag = augmented[col][col];
            for (int j = col; j <= n; ++j) {
                augmented[col][j] /= diag;
            }

            for (int row = 0; row < n; ++row) {
                if (row == col) {
                    continue;
                }
                final double factor = augmented[row][col];
                if (factor == 0.0) {
                    continue;
                }
                for (int j = col; j <= n; ++j) {
                    augmented[row][j] -= factor * augmented[col][j];
                }
            }
        }
    }

    private OUCovarianceGradientWorkspace workspace(final int dimension) {
        OUCovarianceGradientWorkspace workspace = workspaceLocal.get();
        if (workspace == null || workspace.dimension != dimension) {
            workspace = new OUCovarianceGradientWorkspace(dimension, 17, 0, true);
            workspaceLocal.set(workspace);
        }
        return workspace;
    }
}
