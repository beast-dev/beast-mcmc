package test.dr.inference.timeseries;

import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.latent;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.representation;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.model.ScaledIdentityMatrixParameter;
import dr.inference.timeseries.core.BasicTimeSeriesModel;
import dr.inference.timeseries.core.IrregularTimeGrid;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.engine.kalman.GaussianForwardComputationMode;
import dr.inference.timeseries.likelihood.GaussianGradientComputationMode;
import dr.inference.timeseries.likelihood.GaussianSmootherComputationMode;
import dr.inference.timeseries.likelihood.GaussianTimeSeriesLikelihoodFactory;
import dr.inference.timeseries.likelihood.ParallelTimeSeriesLikelihood;
import dr.inference.timeseries.likelihood.TimeSeriesLikelihood;
import dr.inference.timeseries.model.gaussian.LinearGaussianObservationModel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import java.util.Arrays;

/**
 * Partially observed columns with mixed occupancy (0, 1, 2 or all coordinates observed)
 * on an irregular time grid.
 *
 * <p>The canonical likelihood is compared with a dense joint-Gaussian oracle that stacks
 * every latent state, marginalizes nothing analytically, and evaluates the observed cells
 * directly. Every analytical gradient (drift, diffusion, stationary mean, observation-noise
 * variance) is compared with fourth-order central finite differences of the same likelihood.
 * A two-series parallel aggregate with different missingness patterns checks the shared
 * canonical schedule fallback.</p>
 */
public class PartialObservationIrregularGridTest extends TestCase {

    private static final double NAN = Double.NaN;

    public PartialObservationIrregularGridTest(final String name) {
        super(name);
    }

    public static Test suite() {
        return new TestSuite(PartialObservationIrregularGridTest.class);
    }

    // ------------------------------------------------------------------ likelihood oracle

    public void testModerateNoiseLikelihoodMatchesDenseJointGaussianOracle() {
        final Fixture f = moderateFixture("oracle.moderate");
        final double oracle = denseJointGaussianLogLikelihood(f);
        assertEquals("canonical partial-observation likelihood vs dense oracle (moderate R)",
                oracle, f.likelihood.getLogLikelihood(), 1e-8);
    }

    public void testTinyNoiseTinyStepLikelihoodMatchesDenseJointGaussianOracle() {
        final Fixture f = bitmexLikeFixture("oracle.tiny");
        final double oracle = denseJointGaussianLogLikelihood(f);
        assertEquals("canonical partial-observation likelihood vs dense oracle (R = 1e-5 I, sub-second steps)",
                oracle, f.likelihood.getLogLikelihood(), 1e-7);
    }

    // ------------------------------------------------------------------ gradients, moderate

    public void testModerateNoiseDriftGradientMatchesFiniteDifference() {
        final Fixture f = moderateFixture("fd.moderate.drift");
        final double[] analytical = f.likelihood.getGradientWrt(f.drift).getGradientLogDensity(null);
        for (int i = 0; i < analytical.length; ++i) {
            assertEquals("drift gradient idx=" + i,
                    fdFourthOrder(f.likelihood, f.drift, i, 1e-3), analytical[i], 1e-6);
        }
    }

    public void testModerateNoiseDiffusionGradientMatchesFiniteDifference() {
        final Fixture f = moderateFixture("fd.moderate.diffusion");
        final double[] analytical = f.likelihood.getGradientWrt(f.diffusion).getGradientLogDensity(null);
        for (int i = 0; i < analytical.length; ++i) {
            assertEquals("diffusion gradient idx=" + i,
                    fdFourthOrder(f.likelihood, f.diffusion, i, 1e-4), analytical[i], 1e-6);
        }
    }

    public void testModerateNoiseMeanGradientMatchesFiniteDifference() {
        final Fixture f = moderateFixture("fd.moderate.mean");
        final double[] analytical = f.likelihood.getGradientWrt(f.mean).getGradientLogDensity(null);
        for (int i = 0; i < analytical.length; ++i) {
            assertEquals("mean gradient idx=" + i,
                    fdFourthOrder(f.likelihood, f.mean, i, 1e-3), analytical[i], 1e-7);
        }
    }

    public void testModerateNoiseDiagonalObservationNoiseGradientMatchesFiniteDifference() {
        final Fixture f = moderateFixture("fd.moderate.noise");
        final double[] analytical = f.likelihood.getGradientWrt(f.noiseVariance).getGradientLogDensity(null);
        assertEquals(3, analytical.length);
        for (int i = 0; i < analytical.length; ++i) {
            assertEquals("diagonal observation-noise gradient idx=" + i,
                    fdFourthOrder(f.likelihood, f.noiseVariance, i, 1e-4), analytical[i], 1e-6);
        }
    }

    // ------------------------------------------------------------------ gradients, tiny R / tiny dt

    public void testTinyNoiseDriftAndMeanGradientsMatchFiniteDifference() {
        final Fixture f = bitmexLikeFixture("fd.tiny.driftmean");
        final double[] drift = f.likelihood.getGradientWrt(f.drift).getGradientLogDensity(null);
        for (int i = 0; i < drift.length; ++i) {
            final double fd = fdFourthOrder(f.likelihood, f.drift, i, 1e-2);
            assertEquals("tiny-step drift gradient idx=" + i, fd, drift[i],
                    Math.max(1e-5, 1e-5 * Math.abs(fd)));
        }
        final double[] mean = f.likelihood.getGradientWrt(f.mean).getGradientLogDensity(null);
        for (int i = 0; i < mean.length; ++i) {
            final double fd = fdFourthOrder(f.likelihood, f.mean, i, 1e-2);
            assertEquals("tiny-step mean gradient idx=" + i, fd, mean[i],
                    Math.max(1e-5, 1e-5 * Math.abs(fd)));
        }
    }

    public void testTinyNoiseScalarObservationNoiseGradientMatchesLogScaleFiniteDifference() {
        final Fixture f = bitmexLikeFixture("fd.tiny.noise");
        final double[] analytical = f.likelihood.getGradientWrt(f.noiseVariance).getGradientLogDensity(null);
        assertEquals(1, analytical.length);
        final double analyticalLogScale = f.noiseVariance.getParameterValue(0) * analytical[0];
        final double fdLogScale = fdLogScaleFourthOrder(f.likelihood, f.noiseVariance, 0, 0.02);
        assertEquals("scalar observation-noise log-scale gradient (R = 1e-5)",
                fdLogScale, analyticalLogScale, Math.max(1e-6, 1e-4 * Math.abs(fdLogScale)));
    }

    // ------------------------------------------------------------------ parallel, different masks

    public void testParallelAggregateWithDifferentMissingnessPatternsSumsSeries() {
        final Fixture a = moderateFixture("parallel.a");
        // Second series: same process object, same grid length, different mask and values.
        final double[][] otherY = {
                {0.11, NAN, 0.02, NAN, NAN, 0.27, NAN, -0.05, NAN},
                {NAN, 0.19, NAN, NAN, -0.22, NAN, 0.14, NAN, 0.07},
                {NAN, -0.08, 0.16, NAN, 0.03, NAN, NAN, NAN, NAN}
        };
        final LinearGaussianObservationModel obsB = new LinearGaussianObservationModel(
                "parallel.b.obs", 3, a.design, a.noise, matrix("parallel.b.Y", otherY));
        final TimeSeriesLikelihood b = canonicalLikelihood("parallel.b", a.process, obsB, a.grid);

        final ParallelTimeSeriesLikelihood aggregate = new ParallelTimeSeriesLikelihood(
                "parallel.aggregate", 2, Arrays.asList(a.likelihood, b));

        for (int pass = 0; pass < 3; ++pass) {
            if (pass > 0) {
                // Move the shared drift so every cache is rebuilt and the schedule re-enabled.
                a.drift.setParameterValue(0, 0, a.drift.getParameterValue(0, 0) + 0.01 * pass);
                a.likelihood.makeDirty();
                b.makeDirty();
                aggregate.makeDirty();
            }
            assertEquals("aggregate log-likelihood pass " + pass,
                    a.likelihood.getLogLikelihood() + b.getLogLikelihood(),
                    aggregate.getLogLikelihood(), 1e-10);
            final double[] ga = a.likelihood.getGradientWrt(a.drift).getGradientLogDensity(null);
            final double[] gb = b.getGradientWrt(a.drift).getGradientLogDensity(null);
            final double[] gAgg = aggregate.getGradientWrt(a.drift).getGradientLogDensity(null);
            for (int i = 0; i < ga.length; ++i) {
                assertEquals("aggregate drift gradient pass " + pass + " idx=" + i,
                        ga[i] + gb[i], gAgg[i], 1e-9);
            }
        }
    }

    // ------------------------------------------------------------------ fixtures

    private static final class Fixture {
        OUProcessModel process;
        MatrixParameter drift;
        MatrixParameter diffusion;
        Parameter mean;
        Parameter noiseVariance;
        MatrixParameter design;
        MatrixParameter noise;
        MatrixParameter observations;
        TimeGrid grid;
        LinearGaussianObservationModel observation;
        TimeSeriesLikelihood likelihood;
    }

    /** 3D OU, non-identity H, diagonal R of moderate size, irregular grid, occupancy 0/1/2/3. */
    private static Fixture moderateFixture(final String p) {
        final Fixture f = new Fixture();
        f.drift = matrix(p + ".A", new double[][]{
                {0.90, -0.20, 0.10},
                {0.15, 1.30, -0.10},
                {0.05, 0.20, 0.70}
        });
        f.diffusion = matrix(p + ".Q", new double[][]{
                {0.90, 0.06, 0.02},
                {0.06, 0.82, 0.04},
                {0.02, 0.04, 1.04}
        });
        f.mean = new Parameter.Default(p + ".mu", new double[]{0.12, -0.07, 0.18});
        final MatrixParameter p0 = matrix(p + ".P0", scaledIdentity(3, 0.05));
        f.process = new OUProcessModel(p + ".ou", 3, f.drift, f.diffusion, f.mean, p0,
                OUProcessModel.CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        f.design = matrix(p + ".H", new double[][]{
                {1.00, 0.10, 0.00},
                {0.00, 1.00, -0.20},
                {0.05, 0.00, 1.00}
        });
        f.noiseVariance = new Parameter.Default(p + ".obsVar", new double[]{0.02, 0.05, 0.03});
        f.noise = new DiagonalMatrix(f.noiseVariance);
        // masks: t0 {}, t1 {0}, t2 {1,2}, t3 {0,1,2}, t4 {}, t5 {2}, t6 {0,2}, t7 {1}, t8 {0}
        f.observations = matrix(p + ".Y", new double[][]{
                {NAN, 0.21, NAN, 0.33, NAN, NAN, -0.12, NAN, 0.08},
                {NAN, NAN, -0.17, 0.05, NAN, NAN, NAN, 0.26, NAN},
                {NAN, NAN, 0.09, -0.24, NAN, 0.31, 0.14, NAN, NAN}
        });
        f.grid = new IrregularTimeGrid(new double[]{
                0.0, 0.0004, 0.0011, 0.05, 0.051, 0.3, 0.9, 0.9007, 1.4});
        f.observation = new LinearGaussianObservationModel(p + ".obs", 3, f.design, f.noise, f.observations);
        f.likelihood = canonicalLikelihood(p, f.process, f.observation, f.grid);
        return f;
    }

    /** 3D OU, H = I, R = 1e-5 I, sub-second (hour-unit) irregular steps, all-missing anchor at t = 0. */
    private static Fixture bitmexLikeFixture(final String p) {
        final Fixture f = new Fixture();
        f.drift = matrix(p + ".A", new double[][]{
                {0.90, -0.20, 0.10},
                {0.15, 1.30, -0.10},
                {0.05, 0.20, 0.70}
        });
        f.diffusion = matrix(p + ".Q", new double[][]{
                {0.40, 0.03, 0.00},
                {0.03, 0.30, 0.02},
                {0.00, 0.02, 0.50}
        });
        f.mean = new Parameter.Default(p + ".mu", new double[]{0.01, -0.02, 0.00});
        final MatrixParameter p0 = matrix(p + ".P0", scaledIdentity(3, 0.05));
        f.process = new OUProcessModel(p + ".ou", 3, f.drift, f.diffusion, f.mean, p0,
                OUProcessModel.CovarianceGradientMethod.LYAPUNOV_ADJOINT);

        f.design = matrix(p + ".H", scaledIdentity(3, 1.0));
        f.noiseVariance = new Parameter.Default(p + ".obsVar", new double[]{1.0e-5});
        f.noise = new ScaledIdentityMatrixParameter(p + ".R", 3, f.noiseVariance);
        // masks: t0 {}, t1 {0}, t2 {2}, t3 {0,1}, t4 {1}, t5 {2}, t6 {0}, t7 {1,2}
        f.observations = matrix(p + ".Y", new double[][]{
                {NAN, 0.013, NAN, 0.011, NAN, NAN, 0.009, NAN},
                {NAN, NAN, NAN, -0.021, -0.024, NAN, NAN, -0.019},
                {NAN, NAN, 0.004, NAN, NAN, 0.006, NAN, 0.002}
        });
        f.grid = new IrregularTimeGrid(new double[]{
                0.0, 1.6389e-4, 2.0e-4, 2.0278e-4, 9.1e-4, 1.3e-3, 1.3e-3 + 2.8e-7, 4.5e-3});
        f.observation = new LinearGaussianObservationModel(p + ".obs", 3, f.design, f.noise, f.observations);
        f.likelihood = canonicalLikelihood(p, f.process, f.observation, f.grid);
        return f;
    }

    private static TimeSeriesLikelihood canonicalLikelihood(final String name,
                                                            final OUProcessModel process,
                                                            final LinearGaussianObservationModel obs,
                                                            final TimeGrid grid) {
        final BasicTimeSeriesModel model = new BasicTimeSeriesModel(name + ".model", latent(process), obs, grid);
        return GaussianTimeSeriesLikelihoodFactory.create(
                name + ".likelihood",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);
    }

    // ------------------------------------------------------------------ dense oracle

    /**
     * Stacks x_0..x_{T-1} into one Gaussian via the transition representation, then evaluates
     * the observed cells y_{r,t} = H_r x_t + e_{r,t} directly as a dense multivariate normal.
     */
    private static double denseJointGaussianLogLikelihood(final Fixture f) {
        final int d = f.process.getStateDimension();
        final int T = f.grid.getTimeCount();
        final GaussianTransitionRepresentation repr =
                representation(f.process, GaussianTransitionRepresentation.class);

        final double[][] mean = new double[T][d];
        final double[][][][] cov = new double[T][T][][];
        {
            final double[] m0 = new double[d];
            final double[] p0 = new double[d * d];
            repr.getInitialMean(m0);
            repr.getInitialCovarianceFlat(p0);
            mean[0] = m0;
            cov[0][0] = unflatten(p0, d);
        }
        final double[] F = new double[d * d];
        final double[] c = new double[d];
        final double[] Q = new double[d * d];
        for (int t = 0; t < T - 1; ++t) {
            repr.getTransitionMatrixFlat(t, t + 1, f.grid, F);
            repr.getTransitionOffset(t, t + 1, f.grid, c);
            repr.getTransitionCovarianceFlat(t, t + 1, f.grid, Q);
            final double[][] Fm = unflatten(F, d);
            final double[][] Qm = unflatten(Q, d);
            mean[t + 1] = add(mul(Fm, mean[t]), c);
            for (int s = 0; s <= t; ++s) {
                cov[t + 1][s] = mul(Fm, cov[t][s]);
                cov[s][t + 1] = transpose(cov[t + 1][s]);
            }
            cov[t + 1][t + 1] = add(mul(mul(Fm, cov[t][t]), transpose(Fm)), Qm);
        }

        // observed cells
        final int obsDim = f.observation.getObservationDimension();
        final double[] H = new double[obsDim * d];
        final double[] R = new double[obsDim * obsDim];
        f.observation.fillDesignMatrixFlat(H, d);
        f.observation.fillNoiseCovarianceFlat(R);
        int n = 0;
        final int[] cellTime = new int[obsDim * T];
        final int[] cellRow = new int[obsDim * T];
        final double[] cellValue = new double[obsDim * T];
        for (int t = 0; t < T; ++t) {
            for (int r = 0; r < obsDim; ++r) {
                final double y = f.observations.getParameterValue(r, t);
                if (!Double.isNaN(y)) {
                    cellTime[n] = t;
                    cellRow[n] = r;
                    cellValue[n] = y;
                    ++n;
                }
            }
        }
        final double[] residual = new double[n];
        final double[] C = new double[n * n];
        for (int k = 0; k < n; ++k) {
            double mu = 0.0;
            for (int j = 0; j < d; ++j) {
                mu += H[cellRow[k] * d + j] * mean[cellTime[k]][j];
            }
            residual[k] = cellValue[k] - mu;
            for (int l = 0; l < n; ++l) {
                double v = 0.0;
                final double[][] S = cov[cellTime[k]][cellTime[l]];
                for (int i = 0; i < d; ++i) {
                    final double hk = H[cellRow[k] * d + i];
                    if (hk == 0.0) continue;
                    for (int j = 0; j < d; ++j) {
                        v += hk * S[i][j] * H[cellRow[l] * d + j];
                    }
                }
                if (cellTime[k] == cellTime[l]) {
                    v += R[cellRow[k] * obsDim + cellRow[l]];
                }
                C[k * n + l] = v;
            }
        }
        final double[] L = new double[n * n];
        final double[] Linv = new double[n * n];
        final double[] Cinv = new double[n * n];
        assertTrue("oracle covariance must be positive definite", MatrixOps.tryCholesky(C, L, n));
        final double logDet = MatrixOps.invertFromCholesky(L, Linv, Cinv, n);
        final double[] tmp = new double[n];
        final double quad = MatrixOps.quadraticForm(residual, Cinv, n, tmp);
        return -0.5 * (n * Math.log(2.0 * Math.PI) + logDet + quad);
    }

    // ------------------------------------------------------------------ helpers

    private static double fdFourthOrder(final dr.inference.model.Likelihood like,
                                        final Parameter param, final int index, final double h) {
        final double x = param.getParameterValue(index);
        final double p2 = evalAt(like, param, index, x + 2.0 * h);
        final double p1 = evalAt(like, param, index, x + h);
        final double m1 = evalAt(like, param, index, x - h);
        final double m2 = evalAt(like, param, index, x - 2.0 * h);
        evalAt(like, param, index, x);
        return (-p2 + 8.0 * p1 - 8.0 * m1 + m2) / (12.0 * h);
    }

    private static double fdLogScaleFourthOrder(final dr.inference.model.Likelihood like,
                                                final Parameter param, final int index, final double h) {
        final double x = param.getParameterValue(index);
        final double p2 = evalAt(like, param, index, x * Math.exp(2.0 * h));
        final double p1 = evalAt(like, param, index, x * Math.exp(h));
        final double m1 = evalAt(like, param, index, x * Math.exp(-h));
        final double m2 = evalAt(like, param, index, x * Math.exp(-2.0 * h));
        evalAt(like, param, index, x);
        return (-p2 + 8.0 * p1 - 8.0 * m1 + m2) / (12.0 * h);
    }

    private static double evalAt(final dr.inference.model.Likelihood like,
                                 final Parameter param, final int index, final double value) {
        param.setParameterValue(index, value);
        like.makeDirty();
        return like.getLogLikelihood();
    }

    private static MatrixParameter matrix(final String name, final double[][] values) {
        final MatrixParameter mp = new MatrixParameter(name, values.length, values[0].length);
        for (int i = 0; i < values.length; ++i) {
            for (int j = 0; j < values[0].length; ++j) {
                mp.setParameterValue(i, j, values[i][j]);
            }
        }
        return mp;
    }

    private static double[][] scaledIdentity(final int d, final double s) {
        final double[][] out = new double[d][d];
        for (int i = 0; i < d; ++i) out[i][i] = s;
        return out;
    }

    private static double[][] unflatten(final double[] flat, final int d) {
        final double[][] out = new double[d][d];
        for (int i = 0; i < d; ++i) System.arraycopy(flat, i * d, out[i], 0, d);
        return out;
    }

    private static double[][] mul(final double[][] A, final double[][] B) {
        final int n = A.length, m = B[0].length, k = B.length;
        final double[][] C = new double[n][m];
        for (int i = 0; i < n; ++i)
            for (int l = 0; l < k; ++l) {
                final double a = A[i][l];
                for (int j = 0; j < m; ++j) C[i][j] += a * B[l][j];
            }
        return C;
    }

    private static double[] mul(final double[][] A, final double[] x) {
        final double[] y = new double[A.length];
        for (int i = 0; i < A.length; ++i)
            for (int j = 0; j < x.length; ++j) y[i] += A[i][j] * x[j];
        return y;
    }

    private static double[][] add(final double[][] A, final double[][] B) {
        final double[][] C = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; ++i)
            for (int j = 0; j < A[0].length; ++j) C[i][j] = A[i][j] + B[i][j];
        return C;
    }

    private static double[] add(final double[] a, final double[] b) {
        final double[] c = new double[a.length];
        for (int i = 0; i < a.length; ++i) c[i] = a[i] + b[i];
        return c;
    }

    private static double[][] transpose(final double[][] A) {
        final double[][] T = new double[A[0].length][A.length];
        for (int i = 0; i < A.length; ++i)
            for (int j = 0; j < A[0].length; ++j) T[j][i] = A[i][j];
        return T;
    }
}
