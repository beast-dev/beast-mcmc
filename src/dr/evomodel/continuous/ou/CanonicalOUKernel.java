package dr.evomodel.continuous.ou;

import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianBranchTransitionKernel;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianUtils;
import dr.evomodel.treedatalikelihood.continuous.gaussian.GaussianBranchTransitionKernel;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.timeseries.gaussian.DiffusionMatrixParameterization;

/**
 * Canonical OU transition and gradient kernel.
 *
 * <p>This class owns the reusable numerical workspace for canonical tree
 * transitions and gradients. BEAST model wrappers provide live parameter
 * access through small parameterization/provider interfaces.</p>
 */
public final class CanonicalOUKernel
        implements GaussianBranchTransitionKernel, CanonicalGaussianBranchTransitionKernel {

    private static final String DIFFUSION_GRADIENT_QUADRATURE_SUBSTEPS_PROPERTY =
            "beast.experimental.ouDiffusionGradientQuadratureSubsteps";

    public interface InitialMoments {
        void getInitialMean(double[] out);

        void getInitialCovariance(double[][] out);
    }

    private final int stateDimension;
    private final SelectionMatrixParameterization selectionMatrixParameterization;
    private final DiffusionMatrixParameterization diffusionMatrixParameterization;
    private final MatrixParameterInterface diffusionMatrix;
    private final OUCovarianceGradientStrategy covarianceGradientStrategy;
    private final InitialMoments initialMoments;
    private final ThreadLocal<Workspace> workspaceLocal;

    CanonicalOUKernel(final int stateDimension,
                      final SelectionMatrixParameterization selectionMatrixParameterization,
                      final DiffusionMatrixParameterization diffusionMatrixParameterization,
                      final MatrixParameterInterface diffusionMatrix,
                      final OUCovarianceGradientStrategy covarianceGradientStrategy,
                      final InitialMoments initialMoments) {
        this.stateDimension = stateDimension;
        this.selectionMatrixParameterization = selectionMatrixParameterization;
        this.diffusionMatrixParameterization = diffusionMatrixParameterization;
        this.diffusionMatrix = diffusionMatrix;
        this.covarianceGradientStrategy = covarianceGradientStrategy;
        this.initialMoments = initialMoments;
        this.workspaceLocal = new ThreadLocal<Workspace>() {
            @Override
            protected Workspace initialValue() {
                return new Workspace(stateDimension);
            }
        };
    }

    @Override
    public int getStateDimension() {
        return stateDimension;
    }

    @Override
    public void fillInitialCanonicalState(final CanonicalGaussianState out) {
        final Workspace workspace = workspace();
        final double[] mean = workspace.vector0;
        final double[][] covariance = workspace.squareMatrices[0];
        getInitialMean(mean);
        getInitialCovariance(covariance);
        CanonicalGaussianUtils.fillStateFromMoments(mean, covariance, out);
    }

    @Override
    public void fillCanonicalTransition(final double dt, final CanonicalGaussianTransition out) {
        final Workspace workspace = workspace();
        if (selectionMatrixParameterization instanceof CanonicalPreparedTransitionCapability) {
            final double[] mean = workspace.vector0;
            getInitialMean(mean);
            ((CanonicalPreparedTransitionCapability) selectionMatrixParameterization)
                    .fillCanonicalTransition(diffusionMatrix, mean, dt, out);
            return;
        }
        final double[][] transitionMatrix = workspace.squareMatrices[0];
        final double[] transitionOffset = workspace.vector1;
        final double[][] transitionCovariance = workspace.squareMatrices[1];
        fillTransitionMatrix(dt, transitionMatrix);
        fillTransitionOffset(dt, transitionOffset);
        fillTransitionCovariance(dt, transitionCovariance);
        CanonicalGaussianUtils.fillTransitionFromMoments(
                transitionMatrix,
                transitionOffset,
                transitionCovariance,
                out);
    }

    @Override
    public void getInitialMean(final double[] out) {
        checkVectorLength(out, stateDimension, "initial mean");
        initialMoments.getInitialMean(out);
    }

    @Override
    public void getInitialCovariance(final double[][] out) {
        checkSquareMatrix(out, stateDimension, "initial covariance");
        initialMoments.getInitialCovariance(out);
    }

    @Override
    public void fillTransitionMatrix(final double dt, final double[][] out) {
        checkSquareMatrix(out, stateDimension, "transition matrix");
        selectionMatrixParameterization.fillTransitionMatrix(dt, out);
    }

    public void fillTransitionMatrixFlat(final double dt, final double[] out) {
        checkFlatSquare(out, stateDimension, "transition matrix");
        selectionMatrixParameterization.fillTransitionMatrixFlat(dt, out);
    }

    @Override
    public void fillTransitionOffset(final double dt, final double[] out) {
        checkVectorLength(out, stateDimension, "transition offset");
        final Workspace workspace = workspace();
        final double[][] transitionMatrix = workspace.squareMatrices[2];
        final double[] mu = workspace.vector0;
        final double[] transformedMean = workspace.vector2;

        fillTransitionMatrix(dt, transitionMatrix);
        getInitialMean(mu);
        MatrixExponentialUtils.multiply(transitionMatrix, mu, transformedMean);

        for (int i = 0; i < stateDimension; ++i) {
            out[i] = mu[i] - transformedMean[i];
        }
    }

    @Override
    public void fillTransitionCovariance(final double dt, final double[][] out) {
        checkSquareMatrix(out, stateDimension, "transition covariance");
        final Workspace workspace = workspace();

        if (selectionMatrixParameterization instanceof CanonicalPreparedTransitionCapability) {
            ((CanonicalPreparedTransitionCapability) selectionMatrixParameterization)
                    .fillTransitionCovariance(diffusionMatrix, dt, out);
            return;
        }

        final double[][] a = workspace.squareMatrices[2];
        final double[][] q = workspace.squareMatrices[3];
        selectionMatrixParameterization.fillSelectionMatrix(a);
        diffusionMatrixParameterization.fillDiffusionMatrix(q);

        final double[][] vanLoan = workspace.blockMatrices[0];
        final double[][] vanLoanExp = workspace.blockMatrices[1];

        for (int i = 0; i < stateDimension; ++i) {
            for (int j = 0; j < stateDimension; ++j) {
                vanLoan[i][j] = -dt * a[i][j];
                vanLoan[i][j + stateDimension] = dt * q[i][j];
                vanLoan[i + stateDimension][j] = 0.0;
                vanLoan[i + stateDimension][j + stateDimension] = dt * a[j][i];
            }
        }

        MatrixExponentialUtils.expm(vanLoan, vanLoanExp);

        for (int i = 0; i < stateDimension; ++i) {
            for (int j = 0; j < stateDimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < stateDimension; ++k) {
                    sum += vanLoanExp[i][k + stateDimension] * vanLoanExp[j][k];
                }
                out[i][j] = sum;
            }
        }
        MatrixExponentialUtils.symmetrize(out);
    }

    @Override
    public void accumulateSelectionGradientFromCovariance(final double dt,
                                                          final double[][] dLogL_dV,
                                                          final double[] gradientAccumulator) {
        final int d = stateDimension;
        final Workspace workspace = workspace();

        final double[][] a = workspace.squareMatrices[0];
        final double[][] q = workspace.squareMatrices[1];
        selectionMatrixParameterization.fillSelectionMatrix(a);
        diffusionMatrixParameterization.fillDiffusionMatrix(q);

        covarianceGradientStrategy.accumulate(dt, d, a, q, dLogL_dV, gradientAccumulator);
    }

    public void accumulateSelectionGradientFromCovarianceFlat(final double dt,
                                                              final double[] dLogL_dV,
                                                              final boolean transposeAdjoint,
                                                              final double[] gradientAccumulator) {
        checkFlatSquare(dLogL_dV, stateDimension, "covariance adjoint");
        final int d = stateDimension;
        final Workspace workspace = workspace();

        final double[][] a = workspace.squareMatrices[0];
        final double[][] q = workspace.squareMatrices[1];
        selectionMatrixParameterization.fillSelectionMatrix(a);
        diffusionMatrixParameterization.fillDiffusionMatrix(q);

        covarianceGradientStrategy.accumulateFlat(
                dt, d, a, q, dLogL_dV, transposeAdjoint, gradientAccumulator);
    }

    @Override
    public void accumulateDiffusionGradient(final double dt,
                                            final double[][] dLogL_dV,
                                            final double[] gradientAccumulator) {
        final int d = stateDimension;
        final Workspace workspace = workspace();
        final double[][] a = workspace.squareMatrices[0];
        final double[][] gSym = workspace.squareMatrices[2];
        final double[][] fS = workspace.squareMatrices[3];
        final double[][] fST = workspace.squareMatrices[4];
        final double[][] tempDxD = workspace.squareMatrices[5];
        final double[][] contrib = workspace.squareMatrices[6];
        final double[][] expScratch = workspace.squareMatrices[7];

        selectionMatrixParameterization.fillSelectionMatrix(a);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                gSym[i][j] = 0.5 * (dLogL_dV[i][j] + dLogL_dV[j][i]);
            }
        }

        final int substeps = Math.max(1,
                Integer.getInteger(DIFFUSION_GRADIENT_QUADRATURE_SUBSTEPS_PROPERTY, 2));
        final double h = dt / substeps;
        for (int sub = 0; sub < substeps; ++sub) {
            final double t0 = sub * h;
            for (int idx = 0; idx < OUCovarianceGradientMath.GL5_NODES.length; ++idx) {
                final double s = t0 + 0.5 * h * (OUCovarianceGradientMath.GL5_NODES[idx] + 1.0);
                final double scaledWeight = 0.5 * h * OUCovarianceGradientMath.GL5_WEIGHTS[idx];

                OUCovarianceGradientMath.buildExpmMinusAs(s, a, d, fS, expScratch);
                MatrixExponentialUtils.transpose(fS, fST);
                MatrixExponentialUtils.multiply(fST, gSym, tempDxD);
                MatrixExponentialUtils.multiply(tempDxD, fS, contrib);

                for (int i = 0; i < d; ++i) {
                    for (int j = 0; j < d; ++j) {
                        gradientAccumulator[i * d + j] += scaledWeight * contrib[i][j];
                    }
                }
            }
        }
    }

    public void accumulateDiffusionGradientFlat(final double dt,
                                                final double[] dLogL_dV,
                                                final boolean transposeAdjoint,
                                                final double[] gradientAccumulator) {
        checkFlatSquare(dLogL_dV, stateDimension, "covariance adjoint");
        final int d = stateDimension;
        final Workspace workspace = workspace();
        final double[][] a = workspace.squareMatrices[0];
        final double[][] gSym = workspace.squareMatrices[2];
        final double[][] fS = workspace.squareMatrices[3];
        final double[][] fST = workspace.squareMatrices[4];
        final double[][] tempDxD = workspace.squareMatrices[5];
        final double[][] contrib = workspace.squareMatrices[6];
        final double[][] expScratch = workspace.squareMatrices[7];

        selectionMatrixParameterization.fillSelectionMatrix(a);
        fillSymmetricFromFlat(dLogL_dV, transposeAdjoint, d, gSym);

        final int substeps = Math.max(1,
                Integer.getInteger(DIFFUSION_GRADIENT_QUADRATURE_SUBSTEPS_PROPERTY, 2));
        final double h = dt / substeps;
        for (int sub = 0; sub < substeps; ++sub) {
            final double t0 = sub * h;
            for (int idx = 0; idx < OUCovarianceGradientMath.GL5_NODES.length; ++idx) {
                final double s = t0 + 0.5 * h * (OUCovarianceGradientMath.GL5_NODES[idx] + 1.0);
                final double scaledWeight = 0.5 * h * OUCovarianceGradientMath.GL5_WEIGHTS[idx];

                OUCovarianceGradientMath.buildExpmMinusAs(s, a, d, fS, expScratch);
                MatrixExponentialUtils.transpose(fS, fST);
                MatrixExponentialUtils.multiply(fST, gSym, tempDxD);
                MatrixExponentialUtils.multiply(tempDxD, fS, contrib);

                for (int i = 0; i < d; ++i) {
                    for (int j = 0; j < d; ++j) {
                        gradientAccumulator[i * d + j] += scaledWeight * contrib[i][j];
                    }
                }
            }
        }
    }

    public double contractBranchLengthGradientFlat(final double dt,
                                                   final double[] dLogL_dF,
                                                   final double[] dLogL_df,
                                                   final double[] dLogL_dV,
                                                   final boolean transposeCovarianceAdjoint) {
        checkFlatSquare(dLogL_dF, stateDimension, "transition adjoint");
        checkVectorLength(dLogL_df, stateDimension, "transition offset adjoint");
        checkFlatSquare(dLogL_dV, stateDimension, "covariance adjoint");

        final int d = stateDimension;
        final Workspace workspace = workspace();
        final double[][] a = workspace.squareMatrices[0];
        final double[][] q = workspace.squareMatrices[1];
        final double[][] f = workspace.squareMatrices[2];
        final double[][] dFdt = workspace.squareMatrices[3];
        final double[][] fq = workspace.squareMatrices[4];
        final double[] mu = workspace.vector0;

        selectionMatrixParameterization.fillSelectionMatrix(a);
        diffusionMatrixParameterization.fillDiffusionMatrix(q);
        fillTransitionMatrix(dt, f);
        getInitialMean(mu);

        double score = 0.0;
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                double value = 0.0;
                for (int k = 0; k < d; ++k) {
                    value -= a[i][k] * f[k][j];
                }
                dFdt[i][j] = value;
                score += dLogL_dF[iOff + j] * value;
            }
        }

        for (int i = 0; i < d; ++i) {
            double dfdt = 0.0;
            for (int j = 0; j < d; ++j) {
                dfdt -= dFdt[i][j] * mu[j];
            }
            score += dLogL_df[i] * dfdt;
        }

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double value = 0.0;
                for (int k = 0; k < d; ++k) {
                    value += f[i][k] * q[k][j];
                }
                fq[i][j] = value;
            }
        }

        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                double value = 0.0;
                for (int k = 0; k < d; ++k) {
                    value += fq[i][k] * f[j][k];
                }
                score += flatSquareValue(dLogL_dV, i, j, d, transposeCovarianceAdjoint) * value;
            }
        }

        return score;
    }

    @Override
    public void accumulateSelectionGradient(final double dt,
                                            final double[][] dLogL_dF,
                                            final double[] dLogL_df,
                                            final double[] gradientAccumulator) {
        final Workspace workspace = workspace();
        final double[] mu = workspace.vector0;
        getInitialMean(mu);
        accumulateSelectionGradient(dt, mu, dLogL_dF, dLogL_df, gradientAccumulator);
    }

    public void accumulateSelectionGradientFlat(final double dt,
                                                final double[] dLogL_dF,
                                                final double[] dLogL_df,
                                                final double[] gradientAccumulator) {
        checkFlatSquare(dLogL_dF, stateDimension, "transition adjoint");
        final Workspace workspace = workspace();
        final double[] mu = workspace.vector0;
        getInitialMean(mu);
        accumulateSelectionGradientFlat(dt, mu, dLogL_dF, dLogL_df, gradientAccumulator);
    }

    public void accumulateSelectionGradient(final double dt,
                                            final double[] mean,
                                            final double[][] dLogL_dF,
                                            final double[] dLogL_df,
                                            final double[] gradientAccumulator) {
        selectionMatrixParameterization.accumulateGradientFromTransition(
                dt, mean, dLogL_dF, dLogL_df, gradientAccumulator);
    }

    public void accumulateSelectionGradientFlat(final double dt,
                                                final double[] mean,
                                                final double[] dLogL_dF,
                                                final double[] dLogL_df,
                                                final double[] gradientAccumulator) {
        checkFlatSquare(dLogL_dF, stateDimension, "transition adjoint");
        selectionMatrixParameterization.accumulateGradientFromTransitionFlat(
                dt, mean, dLogL_dF, dLogL_df, gradientAccumulator);
    }

    private Workspace workspace() {
        return workspaceLocal.get();
    }

    private static void checkVectorLength(final double[] vector, final int expectedLength, final String label) {
        if (vector == null || vector.length != expectedLength) {
            throw new IllegalArgumentException(label + " must have length " + expectedLength);
        }
    }

    private static void checkSquareMatrix(final double[][] matrix, final int expectedSize, final String label) {
        if (matrix == null || matrix.length != expectedSize) {
            throw new IllegalArgumentException(label + " must be a " + expectedSize + "x" + expectedSize + " matrix");
        }
        for (double[] row : matrix) {
            if (row == null || row.length != expectedSize) {
                throw new IllegalArgumentException(label + " must be a " + expectedSize + "x" + expectedSize + " matrix");
            }
        }
    }

    private static void checkFlatSquare(final double[] matrix, final int expectedSize, final String label) {
        if (matrix == null || matrix.length != expectedSize * expectedSize) {
            throw new IllegalArgumentException(
                    label + " must have length " + (expectedSize * expectedSize)
                            + " for a " + expectedSize + "x" + expectedSize + " matrix");
        }
    }

    private static double flatSquareValue(final double[] source,
                                          final int row,
                                          final int col,
                                          final int dim,
                                          final boolean transpose) {
        return transpose ? source[col * dim + row] : source[row * dim + col];
    }

    private static void fillSymmetricFromFlat(final double[] source,
                                              final boolean transpose,
                                              final int dim,
                                              final double[][] out) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                out[i][j] = 0.5 * (
                        flatSquareValue(source, i, j, dim, transpose)
                                + flatSquareValue(source, j, i, dim, transpose));
            }
        }
    }

    private static final class Workspace {
        private static final int SQUARE_MATRIX_COUNT = 8;
        private static final int BLOCK_MATRIX_COUNT = 2;

        private final double[][][] squareMatrices;
        private final double[][][] blockMatrices;
        private final double[] vector0;
        private final double[] vector1;
        private final double[] vector2;

        private Workspace(final int dimension) {
            this.squareMatrices = allocateMatrixStack(SQUARE_MATRIX_COUNT, dimension, dimension);
            this.blockMatrices = allocateMatrixStack(BLOCK_MATRIX_COUNT, 2 * dimension, 2 * dimension);
            this.vector0 = new double[dimension];
            this.vector1 = new double[dimension];
            this.vector2 = new double[dimension];
        }
    }

    private static double[][][] allocateMatrixStack(final int count,
                                                    final int rows,
                                                    final int cols) {
        final double[][][] stack = new double[count][][];
        for (int i = 0; i < count; ++i) {
            stack[i] = allocateMatrix(rows, cols);
        }
        return stack;
    }

    private static double[][] allocateMatrix(final int rows,
                                             final int cols) {
        final double[][] matrix = new double[rows][];
        for (int i = 0; i < rows; ++i) {
            matrix[i] = new double[cols];
        }
        return matrix;
    }
}
