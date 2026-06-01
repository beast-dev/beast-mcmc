package dr.evomodel.continuous.ou;

public final class StationaryLyapunovCovarianceGradientStrategy implements OUCovarianceGradientStrategy {

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
        final int gV = workspace.squareOffset(0);
        OUCovarianceGradientMath.fillSymmetricFromFlat(
                covarianceAdjoint, transposeAdjoint, dimension, square, gV);
        accumulateSymmetric(
                dt,
                dimension,
                selection,
                0,
                diffusion,
                0,
                square,
                gV,
                workspace,
                gradientAccumulator);
    }

    private void accumulateSymmetric(final double dt,
                                     final int dimension,
                                     final double[] selection,
                                     final int selectionOffset,
                                     final double[] diffusion,
                                     final int diffusionOffset,
                                     final double[] gV,
                                     final int gVOffset,
                                     final OUCovarianceGradientWorkspace workspace,
                                     final double[] gradientAccumulator) {
        final double[] square = workspace.squareMatrices;
        final int sStat = workspace.squareOffset(1);
        final int f = workspace.squareOffset(2);
        final int tempDxD = workspace.squareOffset(3);
        final int gS = workspace.squareOffset(4);
        final int y = workspace.squareOffset(5);
        final int yT = workspace.squareOffset(6);
        final int gAStationary = workspace.squareOffset(7);
        final int sStatT = workspace.squareOffset(8);
        final int gFCov = workspace.squareOffset(9);
        final int gVT = workspace.squareOffset(10);
        final int tempDxD2 = workspace.squareOffset(11);
        final int gX = workspace.squareOffset(12);
        final int fT = workspace.squareOffset(13);
        final int expScratch = workspace.squareOffset(14);
        final int minusAdt = workspace.squareOffset(15);
        final int aT = workspace.squareOffset(16);

        solveStationaryLyapunovFlat(
                selection, selectionOffset, diffusion, diffusionOffset, dimension, square, sStat, workspace);
        OUCovarianceGradientMath.buildExpmMinusAs(
                dt, selection, selectionOffset, dimension, square, f, square, expScratch);

        MatrixExponentialUtils.transposeFlat(square, f, square, fT, dimension);
        MatrixExponentialUtils.multiplyFlat(square, fT, gV, gVOffset, square, tempDxD, dimension);
        MatrixExponentialUtils.multiplyFlat(square, tempDxD, square, f, square, gS, dimension);
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                square[gS + rowOffset + j] =
                        gV[gVOffset + rowOffset + j] - square[gS + rowOffset + j];
            }
        }

        solveAdjointStationaryLyapunovFlat(
                selection, selectionOffset, square, gS, dimension, square, y, square, aT, workspace);

        MatrixExponentialUtils.transposeFlat(square, y, square, yT, dimension);
        MatrixExponentialUtils.multiplyFlat(square, y, square, sStat, square, gAStationary, dimension);
        MatrixExponentialUtils.multiplyFlat(square, yT, square, sStat, square, tempDxD, dimension);
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                square[gAStationary + rowOffset + j] =
                        -(square[gAStationary + rowOffset + j] + square[tempDxD + rowOffset + j]);
            }
        }

        MatrixExponentialUtils.transposeFlat(square, sStat, square, sStatT, dimension);
        MatrixExponentialUtils.multiplyFlat(gV, gVOffset, square, f, square, tempDxD, dimension);
        MatrixExponentialUtils.multiplyFlat(square, tempDxD, square, sStat, square, gFCov, dimension);
        MatrixExponentialUtils.transposeFlat(gV, gVOffset, square, gVT, dimension);
        MatrixExponentialUtils.multiplyFlat(square, gVT, square, f, square, tempDxD2, dimension);
        MatrixExponentialUtils.multiplyFlat(square, tempDxD2, square, sStatT, square, tempDxD, dimension);
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                square[gFCov + rowOffset + j] =
                        -(square[gFCov + rowOffset + j] + square[tempDxD + rowOffset + j]);
            }
        }

        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                square[minusAdt + rowOffset + j] =
                        -dt * selection[selectionOffset + rowOffset + j];
            }
        }
        MatrixExponentialUtils.adjointExpFlat(square, minusAdt, square, gFCov, square, gX, dimension);

        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                gradientAccumulator[rowOffset + j] +=
                        square[gAStationary + j * dimension + i] - dt * square[gX + rowOffset + j];
            }
        }
    }

    private static void solveStationaryLyapunovFlat(final double[] selection,
                                                    final int selectionOffset,
                                                    final double[] diffusion,
                                                    final int diffusionOffset,
                                                    final int dimension,
                                                    final double[] out,
                                                    final int outOffset,
                                                    final OUCovarianceGradientWorkspace workspace) {
        final int n = dimension * dimension;
        final double[] augmented = workspace.lyapunovAugmented;
        final int columns = workspace.lyapunovAugmentedColumns;
        if (augmented == null) {
            throw new IllegalStateException("Stationary Lyapunov workspace is not available");
        }

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                final int row = i * dimension + j;
                final int augmentedRow = row * columns;
                for (int p = 0; p < dimension; ++p) {
                    for (int r = 0; r < dimension; ++r) {
                        final int col = p * dimension + r;
                        double value = 0.0;
                        if (i == p) {
                            value += selection[selectionOffset + j * dimension + r];
                        }
                        if (j == r) {
                            value += selection[selectionOffset + i * dimension + p];
                        }
                        augmented[augmentedRow + col] = value;
                    }
                }
                augmented[augmentedRow + n] = diffusion[diffusionOffset + i * dimension + j];
            }
        }

        try {
            solveAugmentedLinearSystemInPlace(augmented, n, columns);
        } catch (IllegalStateException e) {
            throw new IllegalArgumentException(
                    "Stationary covariance is only defined for stable drift matrices", e);
        }

        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                out[outOffset + rowOffset + j] = augmented[(rowOffset + j) * columns + n];
            }
        }
        MatrixExponentialUtils.symmetrizeFlat(out, outOffset, dimension);
    }

    private static void solveAdjointStationaryLyapunovFlat(final double[] selection,
                                                           final int selectionOffset,
                                                           final double[] gS,
                                                           final int gSOffset,
                                                           final int dimension,
                                                           final double[] yOut,
                                                           final int yOutOffset,
                                                           final double[] aTranspose,
                                                           final int aTransposeOffset,
                                                           final OUCovarianceGradientWorkspace workspace) {
        MatrixExponentialUtils.transposeFlat(selection, selectionOffset, aTranspose, aTransposeOffset, dimension);
        solveStationaryLyapunovFlat(
                aTranspose, aTransposeOffset, gS, gSOffset, dimension, yOut, yOutOffset, workspace);
    }

    private static void solveAugmentedLinearSystemInPlace(final double[] augmented,
                                                          final int n,
                                                          final int columns) {
        for (int col = 0; col < n; ++col) {
            int pivot = col;
            double maxAbs = Math.abs(augmented[col * columns + col]);
            for (int row = col + 1; row < n; ++row) {
                final double abs = Math.abs(augmented[row * columns + col]);
                if (abs > maxAbs) {
                    maxAbs = abs;
                    pivot = row;
                }
            }

            if (maxAbs <= 1e-14) {
                throw new IllegalStateException("Singular system");
            }

            if (pivot != col) {
                swapRows(augmented, columns, col, pivot);
            }

            final int pivotRow = col * columns;
            final double diag = augmented[pivotRow + col];
            for (int j = col; j <= n; ++j) {
                augmented[pivotRow + j] /= diag;
            }

            for (int row = 0; row < n; ++row) {
                if (row == col) {
                    continue;
                }
                final int rowOffset = row * columns;
                final double factor = augmented[rowOffset + col];
                if (factor == 0.0) {
                    continue;
                }
                for (int j = col; j <= n; ++j) {
                    augmented[rowOffset + j] -= factor * augmented[pivotRow + j];
                }
            }
        }
    }

    private static void swapRows(final double[] matrix,
                                 final int columns,
                                 final int rowA,
                                 final int rowB) {
        final int offsetA = rowA * columns;
        final int offsetB = rowB * columns;
        for (int j = 0; j < columns; ++j) {
            final double tmp = matrix[offsetA + j];
            matrix[offsetA + j] = matrix[offsetB + j];
            matrix[offsetB + j] = tmp;
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
