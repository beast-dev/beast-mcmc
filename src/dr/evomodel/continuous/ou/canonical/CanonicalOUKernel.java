package dr.evomodel.continuous.ou.canonical;

import dr.evomodel.continuous.ou.DiffusionMatrixParameterization;
import dr.evomodel.continuous.ou.MatrixExponentialUtils;
import dr.evomodel.continuous.ou.OUCovarianceGradientMath;
import dr.evomodel.continuous.ou.OUCovarianceGradientStrategy;
import dr.evomodel.continuous.ou.SelectionMatrixParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianBranchTransitionKernel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.inference.model.MatrixParameterInterface;

/**
 * Flat-array canonical OU transition and gradient kernel.
 */
public final class CanonicalOUKernel implements CanonicalGaussianBranchTransitionKernel {

    private static final String DIFFUSION_GRADIENT_QUADRATURE_SUBSTEPS_PROPERTY =
            "beast.experimental.ouDiffusionGradientQuadratureSubsteps";

    public interface InitialMoments {
        void getInitialMean(double[] out);

        void getInitialCovarianceFlat(double[] out);
    }

    private final int stateDimension;
    private final SelectionMatrixParameterization selectionMatrixParameterization;
    private final DiffusionMatrixParameterization diffusionMatrixParameterization;
    private final MatrixParameterInterface diffusionMatrix;
    private final OUCovarianceGradientStrategy covarianceGradientStrategy;
    private final InitialMoments initialMoments;
    private final ThreadLocal<Workspace> workspaceLocal;

    public CanonicalOUKernel(final int stateDimension,
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
        getInitialMean(workspace.vector0);
        initialMoments.getInitialCovarianceFlat(workspace.matrix0);
        GaussianFormConverter.fillStateFromMoments(
                workspace.vector0,
                workspace.matrix0,
                stateDimension,
                workspace.converterWorkspace,
                out);
    }

    @Override
    public void fillCanonicalTransition(final double dt, final CanonicalGaussianTransition out) {
        final Workspace workspace = workspace();
        if (selectionMatrixParameterization instanceof CanonicalPreparedTransitionCapability) {
            getInitialMean(workspace.vector0);
            ((CanonicalPreparedTransitionCapability) selectionMatrixParameterization)
                    .fillCanonicalTransition(diffusionMatrix, workspace.vector0, dt, out);
            return;
        }

        fillTransitionMomentsFlat(dt, workspace.matrix0, workspace.vector1, workspace.matrix1);
        GaussianFormConverter.fillTransitionFromMoments(
                workspace.matrix0,
                workspace.vector1,
                workspace.matrix1,
                stateDimension,
                workspace.converterWorkspace,
                out);
    }

    public void getInitialMean(final double[] out) {
        checkVectorLength(out, stateDimension, "initial mean");
        initialMoments.getInitialMean(out);
    }

    public void getInitialCovarianceFlat(final double[] out) {
        checkFlatSquare(out, stateDimension, "initial covariance");
        initialMoments.getInitialCovarianceFlat(out);
    }

    public void fillTransitionMatrixFlat(final double dt, final double[] out) {
        checkFlatSquare(out, stateDimension, "transition matrix");
        selectionMatrixParameterization.fillTransitionMatrixFlat(dt, out);
    }

    public void fillTransitionMomentsFlat(final double dt,
                                          final double[] transitionMatrixOut,
                                          final double[] transitionOffsetOut,
                                          final double[] transitionCovarianceOut) {
        fillTransitionMatrixFlat(dt, transitionMatrixOut);
        fillTransitionOffsetFromMatrix(transitionMatrixOut, transitionOffsetOut);
        fillTransitionCovarianceFlat(dt, transitionCovarianceOut);
    }

    public void fillTransitionCovarianceFlat(final double dt, final double[] out) {
        checkFlatSquare(out, stateDimension, "transition covariance");
        if (selectionMatrixParameterization instanceof CanonicalPreparedTransitionCapability) {
            ((CanonicalPreparedTransitionCapability) selectionMatrixParameterization)
                    .fillTransitionCovarianceFlat(diffusionMatrix, dt, out);
            return;
        }

        final int d = stateDimension;
        final Workspace workspace = workspace();
        final double[] a = workspace.matrix2;
        final double[] q = workspace.matrix3;
        final double[] vanLoan = workspace.block0;
        final double[] vanLoanExp = workspace.block1;

        selectionMatrixParameterization.fillSelectionMatrixFlat(a);
        diffusionMatrixParameterization.fillDiffusionMatrixFlat(q);

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                vanLoan[i * 2 * d + j] = -dt * a[i * d + j];
                vanLoan[i * 2 * d + j + d] = dt * q[i * d + j];
                vanLoan[(i + d) * 2 * d + j] = 0.0;
                vanLoan[(i + d) * 2 * d + j + d] = dt * a[j * d + i];
            }
        }

        MatrixExponentialUtils.expmFlat(vanLoan, vanLoanExp, 2 * d);

        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += vanLoanExp[i * 2 * d + k + d] * vanLoanExp[j * 2 * d + k];
                }
                out[iOff + j] = sum;
            }
        }
        MatrixExponentialUtils.symmetrizeFlat(out, d);
    }

    public void fillTransitionOffset(final double dt, final double[] out) {
        checkVectorLength(out, stateDimension, "transition offset");
        final Workspace workspace = workspace();

        fillTransitionMatrixFlat(dt, workspace.matrix0);
        fillTransitionOffsetFromMatrix(workspace.matrix0, out);
    }

    private void fillTransitionOffsetFromMatrix(final double[] transitionMatrix,
                                                final double[] out) {
        checkFlatSquare(transitionMatrix, stateDimension, "transition matrix");
        checkVectorLength(out, stateDimension, "transition offset");
        final Workspace workspace = workspace();
        getInitialMean(workspace.vector0);
        MatrixOps.matVec(transitionMatrix, workspace.vector0, workspace.vector2, stateDimension);

        for (int i = 0; i < stateDimension; ++i) {
            out[i] = workspace.vector0[i] - workspace.vector2[i];
        }
    }

    public void accumulateSelectionGradientFromCovarianceFlat(final double dt,
                                                              final double[] dLogL_dV,
                                                              final boolean transposeAdjoint,
                                                              final double[] gradientAccumulator) {
        checkFlatSquare(dLogL_dV, stateDimension, "covariance adjoint");
        final Workspace workspace = workspace();
        selectionMatrixParameterization.fillSelectionMatrixFlat(workspace.matrix0);
        diffusionMatrixParameterization.fillDiffusionMatrixFlat(workspace.matrix1);
        covarianceGradientStrategy.accumulateFlatMatrices(
                dt,
                stateDimension,
                workspace.matrix0,
                workspace.matrix1,
                dLogL_dV,
                transposeAdjoint,
                gradientAccumulator);
    }

    public void accumulateDiffusionGradientFlat(final double dt,
                                                final double[] dLogL_dV,
                                                final boolean transposeAdjoint,
                                                final double[] gradientAccumulator) {
        checkFlatSquare(dLogL_dV, stateDimension, "covariance adjoint");
        final int d = stateDimension;
        final Workspace workspace = workspace();
        final double[] a = workspace.matrix0;
        final double[] gSym = workspace.matrix2;
        final double[] fS = workspace.matrix3;
        final double[] fST = workspace.matrix4;
        final double[] temp = workspace.matrix5;
        final double[] contrib = workspace.matrix6;
        final double[] expScratch = workspace.matrix7;

        selectionMatrixParameterization.fillSelectionMatrixFlat(a);
        fillSymmetricFromFlat(dLogL_dV, transposeAdjoint, d, gSym);

        final int substeps = Math.max(1,
                Integer.getInteger(DIFFUSION_GRADIENT_QUADRATURE_SUBSTEPS_PROPERTY, 2));
        final double h = dt / substeps;
        for (int sub = 0; sub < substeps; ++sub) {
            final double t0 = sub * h;
            for (int idx = 0; idx < OUCovarianceGradientMath.GL5_NODES.length; ++idx) {
                final double s = t0 + 0.5 * h * (OUCovarianceGradientMath.GL5_NODES[idx] + 1.0);
                final double scaledWeight = 0.5 * h * OUCovarianceGradientMath.GL5_WEIGHTS[idx];

                buildExpmMinusAsFlat(s, a, d, fS, expScratch);
                MatrixExponentialUtils.transposeFlat(fS, fST, d);
                MatrixExponentialUtils.multiplyFlat(fST, gSym, temp, d);
                MatrixExponentialUtils.multiplyFlat(temp, fS, contrib, d);

                for (int i = 0; i < d * d; ++i) {
                    gradientAccumulator[i] += scaledWeight * contrib[i];
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
        final double[] a = workspace.matrix0;
        final double[] q = workspace.matrix1;
        final double[] f = workspace.matrix2;
        final double[] dFdt = workspace.matrix3;
        final double[] fq = workspace.matrix4;
        final double[] mu = workspace.vector0;

        selectionMatrixParameterization.fillSelectionMatrixFlat(a);
        diffusionMatrixParameterization.fillDiffusionMatrixFlat(q);
        fillTransitionMatrixFlat(dt, f);
        getInitialMean(mu);

        double score = 0.0;
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                double value = 0.0;
                for (int k = 0; k < d; ++k) {
                    value -= a[iOff + k] * f[k * d + j];
                }
                dFdt[iOff + j] = value;
                score += dLogL_dF[iOff + j] * value;
            }
        }

        for (int i = 0; i < d; ++i) {
            double dfdt = 0.0;
            for (int j = 0; j < d; ++j) {
                dfdt -= dFdt[i * d + j] * mu[j];
            }
            score += dLogL_df[i] * dfdt;
        }

        MatrixExponentialUtils.multiplyFlat(f, q, fq, d);
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                double value = 0.0;
                for (int k = 0; k < d; ++k) {
                    value += fq[iOff + k] * f[j * d + k];
                }
                score += flatSquareValue(dLogL_dV, i, j, d, transposeCovarianceAdjoint) * value;
            }
        }

        return score;
    }

    public void accumulateSelectionGradientFlat(final double dt,
                                                final double[] dLogL_dF,
                                                final double[] dLogL_df,
                                                final double[] gradientAccumulator) {
        checkFlatSquare(dLogL_dF, stateDimension, "transition adjoint");
        final Workspace workspace = workspace();
        getInitialMean(workspace.vector0);
        accumulateSelectionGradientFlat(dt, workspace.vector0, dLogL_dF, dLogL_df, gradientAccumulator);
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

    private static void buildExpmMinusAsFlat(final double t,
                                             final double[] a,
                                             final int d,
                                             final double[] out,
                                             final double[] scaledMinusA) {
        for (int i = 0; i < d * d; ++i) {
            scaledMinusA[i] = -t * a[i];
        }
        MatrixExponentialUtils.expmFlat(scaledMinusA, out, d);
    }

    private static void checkVectorLength(final double[] vector, final int expectedLength, final String label) {
        if (vector == null || vector.length != expectedLength) {
            throw new IllegalArgumentException(label + " must have length " + expectedLength);
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
                                              final double[] out) {
        for (int i = 0; i < dim; ++i) {
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                out[iOff + j] = 0.5 * (
                        flatSquareValue(source, i, j, dim, transpose)
                                + flatSquareValue(source, j, i, dim, transpose));
            }
        }
    }

    private static final class Workspace {
        private final double[] vector0;
        private final double[] vector1;
        private final double[] vector2;
        private final double[] matrix0;
        private final double[] matrix1;
        private final double[] matrix2;
        private final double[] matrix3;
        private final double[] matrix4;
        private final double[] matrix5;
        private final double[] matrix6;
        private final double[] matrix7;
        private final double[] block0;
        private final double[] block1;
        private final GaussianFormConverter.Workspace converterWorkspace;

        private Workspace(final int dimension) {
            this.vector0 = new double[dimension];
            this.vector1 = new double[dimension];
            this.vector2 = new double[dimension];
            this.matrix0 = new double[dimension * dimension];
            this.matrix1 = new double[dimension * dimension];
            this.matrix2 = new double[dimension * dimension];
            this.matrix3 = new double[dimension * dimension];
            this.matrix4 = new double[dimension * dimension];
            this.matrix5 = new double[dimension * dimension];
            this.matrix6 = new double[dimension * dimension];
            this.matrix7 = new double[dimension * dimension];
            this.block0 = new double[4 * dimension * dimension];
            this.block1 = new double[4 * dimension * dimension];
            this.converterWorkspace = new GaussianFormConverter.Workspace();
            this.converterWorkspace.ensureDim(dimension);
        }
    }
}
