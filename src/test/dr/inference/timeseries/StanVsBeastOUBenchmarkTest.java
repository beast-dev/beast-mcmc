package test.dr.inference.timeseries;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.BasicTimeSeriesModel;
import dr.inference.timeseries.core.IrregularTimeGrid;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.core.UniformTimeGrid;
import dr.inference.timeseries.engine.gaussian.GaussianForwardComputationMode;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.inference.timeseries.gaussian.OUTimeSeriesProcessAdapter;
import dr.inference.timeseries.likelihood.GaussianGradientComputationMode;
import dr.inference.timeseries.likelihood.GaussianSmootherComputationMode;
import dr.inference.timeseries.likelihood.GaussianTimeSeriesLikelihoodFactory;
import dr.inference.timeseries.likelihood.TimeSeriesLikelihood;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import java.util.Arrays;

/**
 * OU time-series likelihood benchmark designed to be compared with the R/Stan Lyapunov
 * parametrization benchmark in stanOuBenchMark.R.
 *
 * R settings: K=4, T=300, Delta=0.05, sigma_obs=0.10, n_warmup=20, n_timed=100.
 * BEAST uses the orthogonal block-diagonal polar-stable parametrization for the selection
 * matrix A.  Both compute the exact Kalman-filter log-likelihood and its gradient.
 */
public class StanVsBeastOUBenchmarkTest extends TestCase {

    // --- problem dimensions matching the R benchmark ---
    private static final int K       = 10;
    private static final int T_OBS   = 300;
    private static final double DELTA     = 0.05;
    private static final double SIGMA_OBS = 0.10;

    // --- timing settings matching the R benchmark ---
    private static final int N_WARMUP = 20;
    private static final int N_TIMED  = 100;

    public StanVsBeastOUBenchmarkTest(final String name) {
        super(name);
    }

    // -----------------------------------------------------------------------
    // Benchmark tests
    // -----------------------------------------------------------------------

    public void testBenchmarkLogLikelihoodAndGradient() {
        final OUProcessModel process = makeOrthogonalBlockProcess("bench.K" + K);
        final GaussianObservationModel obs = makeObservation("bench.obs");
        final Parameter nativeParam =
                ((OrthogonalBlockDiagonalPolarStableMatrixParameter) process.getDriftMatrix())
                        .getParameter();

        // ── Scenario A: fixed A (parameters don't change between evaluations) ──────────
        // Both grids have all transition matrices pre-cached by prepareTimeGrid();
        // makeDirty() only invalidates the Kalman filter, not the transition cache.
        // This models the cost of evaluating the likelihood at fixed A.

        final UniformTimeGrid uniformGrid = new UniformTimeGrid(T_OBS, 0.0, DELTA);
        final BenchmarkResult uniformFixed = runBenchmark(
                "bench.uniform.fixed", process, obs, uniformGrid, nativeParam);

        final IrregularTimeGrid irregularGrid = makeIrregularGrid();
        final BenchmarkResult irregularFixed = runBenchmark(
                "bench.irregular.fixed", process, obs, irregularGrid, nativeParam);

        // ── Scenario B: A changes every step (MCMC scenario) ────────────────────────
        // A parameter perturbation invalidates the transition-matrix cache before each
        // evaluation, forcing recomputation of exp(-A·Δt). Uniform benefits from caching
        // (1 expm per step); irregular must recompute all T-1 unique transitions.

        final BenchmarkResult uniformMcmc = runBenchmarkWithParamChange(
                "bench.uniform.mcmc", process, obs, uniformGrid, nativeParam);

        final BenchmarkResult irregularMcmc = runBenchmarkWithParamChange(
                "bench.irregular.mcmc", process, obs, irregularGrid, nativeParam);

        // ---- report ----
        System.out.println();
        System.out.println("=== OU benchmark: BEAST (orthogonal-block, canonical Kalman) ===");
        System.out.printf("K=%d  T=%d  sigma_obs=%.3f  n_warmup=%d  n_timed=%d%n",
                K, T_OBS, SIGMA_OBS, N_WARMUP, N_TIMED);
        System.out.println();

        System.out.println("Scenario A – fixed A (transition cache pre-built, A does not change):");
        System.out.printf("  %-32s  %12s  %12s  %10s%n",
                "Operation", "Uniform grid", "Irregular grid", "Slowdown");
        System.out.println("  " + repeat('-', 70));
        printRow("  log_prob (median ms)",
                uniformFixed.logLikNanos, irregularFixed.logLikNanos);
        printRow("  fwd+grad (median ms)",
                uniformFixed.combinedNanos, irregularFixed.combinedNanos);

        System.out.println();
        System.out.println("Scenario B – A changes every call (MCMC: transition cache busted each step):");
        System.out.printf("  %-32s  %12s  %12s  %10s%n",
                "Operation", "Uniform grid", "Irregular grid", "Slowdown");
        System.out.println("  " + repeat('-', 70));
        printRow("  log_prob (median ms)",
                uniformMcmc.logLikNanos, irregularMcmc.logLikNanos);
        printRow("  fwd+grad (median ms)",
                uniformMcmc.combinedNanos, irregularMcmc.combinedNanos);

        System.out.println();
        System.out.printf("  Uniform logL:   %.4f  (1 unique Δt = %.4f)%n",
                uniformFixed.logL, DELTA);
        System.out.printf("  Irregular logL: %.4f  (%d unique Δt values)%n",
                irregularFixed.logL, T_OBS - 1);
        System.out.println();
        System.out.println("  parametrization = orthogonal polar-stable block-diagonal A");
        System.out.println("  fwd+grad corresponds to Stan grad_log_prob");
    }

    // -----------------------------------------------------------------------
    // Core benchmark runner
    // -----------------------------------------------------------------------

    private static BenchmarkResult runBenchmark(final String tag,
                                                final OUProcessModel process,
                                                final GaussianObservationModel obs,
                                                final TimeGrid grid,
                                                final Parameter nativeParam) {
        final OUTimeSeriesProcessAdapter adapter = new OUTimeSeriesProcessAdapter(process);
        final BasicTimeSeriesModel model = new BasicTimeSeriesModel(tag + ".model", adapter, obs, grid);
        final TimeSeriesLikelihood likelihood = GaussianTimeSeriesLikelihoodFactory.create(
                tag + ".lik",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);

        // warm-up: log-likelihood
        for (int i = 0; i < N_WARMUP; ++i) {
            likelihood.makeDirty();
            likelihood.getLogLikelihood();
        }

        // time log-likelihood
        double logL = Double.NaN;
        final long[] logLikNanos = new long[N_TIMED];
        for (int i = 0; i < N_TIMED; ++i) {
            likelihood.makeDirty();
            final long t0 = System.nanoTime();
            logL = likelihood.getLogLikelihood();
            logLikNanos[i] = System.nanoTime() - t0;
        }
        assertTrue("log-likelihood must be finite", Double.isFinite(logL));

        // warm-up: combined fwd+grad
        for (int i = 0; i < N_WARMUP; ++i) {
            likelihood.makeDirty();
            likelihood.getLogLikelihood();
            likelihood.getGradientWrt(nativeParam).getGradientLogDensity(null);
        }

        // time combined forward + gradient
        final long[] combinedNanos = new long[N_TIMED];
        for (int i = 0; i < N_TIMED; ++i) {
            likelihood.makeDirty();
            final long t0 = System.nanoTime();
            likelihood.getLogLikelihood();
            likelihood.getGradientWrt(nativeParam).getGradientLogDensity(null);
            combinedNanos[i] = System.nanoTime() - t0;
        }

        return new BenchmarkResult(logL, logLikNanos, combinedNanos);
    }

    /** Builds a fully-irregular grid with the same T_OBS points and same mean spacing as the
     *  uniform grid, but with random-looking gaps derived from a deterministic sequence. */
    private static IrregularTimeGrid makeIrregularGrid() {
        final double[] times = new double[T_OBS];
        times[0] = 0.0;
        // Vary the gap ±40% around DELTA using a deterministic quasi-random sequence
        for (int i = 1; i < T_OBS; ++i) {
            final double jitter = 0.4 * DELTA * Math.sin(2.7183 * i + 1.4142);
            times[i] = times[i - 1] + DELTA + jitter;
        }
        return new IrregularTimeGrid(times);
    }

    /**
     * Benchmarks the cost when A changes on every call — the true MCMC scenario.
     * A parameter perturbation is applied before each evaluation to invalidate the
     * transition-matrix cache, forcing recomputation of exp(-A·Δt) from scratch.
     */
    private static BenchmarkResult runBenchmarkWithParamChange(final String tag,
                                                               final OUProcessModel process,
                                                               final GaussianObservationModel obs,
                                                               final TimeGrid grid,
                                                               final Parameter nativeParam) {
        final OUTimeSeriesProcessAdapter adapter = new OUTimeSeriesProcessAdapter(process);
        final BasicTimeSeriesModel model = new BasicTimeSeriesModel(tag + ".model", adapter, obs, grid);
        final TimeSeriesLikelihood likelihood = GaussianTimeSeriesLikelihoodFactory.create(
                tag + ".lik",
                model,
                GaussianForwardComputationMode.CANONICAL,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);

        // Use the first rho parameter to trigger cache invalidation — value stays near its original
        final double originalRho0 = nativeParam.getParameterValue(0);
        final double eps = 1e-9;  // negligible change, realistic values unchanged

        // warm-up with parameter perturbation
        for (int i = 0; i < N_WARMUP; ++i) {
            nativeParam.setParameterValue(0, originalRho0 + (i % 2 == 0 ? eps : -eps));
            likelihood.getLogLikelihood();
        }
        nativeParam.setParameterValue(0, originalRho0);

        // time log-likelihood with cache-busting param change
        double logL = Double.NaN;
        final long[] logLikNanos = new long[N_TIMED];
        for (int i = 0; i < N_TIMED; ++i) {
            nativeParam.setParameterValue(0, originalRho0 + (i % 2 == 0 ? eps : -eps));
            final long t0 = System.nanoTime();
            logL = likelihood.getLogLikelihood();
            logLikNanos[i] = System.nanoTime() - t0;
        }
        nativeParam.setParameterValue(0, originalRho0);
        assertTrue("log-likelihood must be finite (param-change, " + tag + ")", Double.isFinite(logL));

        // warm-up combined
        for (int i = 0; i < N_WARMUP; ++i) {
            nativeParam.setParameterValue(0, originalRho0 + (i % 2 == 0 ? eps : -eps));
            likelihood.getLogLikelihood();
            likelihood.getGradientWrt(nativeParam).getGradientLogDensity(null);
        }
        nativeParam.setParameterValue(0, originalRho0);

        // time combined with cache-busting
        final long[] combinedNanos = new long[N_TIMED];
        for (int i = 0; i < N_TIMED; ++i) {
            nativeParam.setParameterValue(0, originalRho0 + (i % 2 == 0 ? eps : -eps));
            final long t0 = System.nanoTime();
            likelihood.getLogLikelihood();
            likelihood.getGradientWrt(nativeParam).getGradientLogDensity(null);
            combinedNanos[i] = System.nanoTime() - t0;
        }
        nativeParam.setParameterValue(0, originalRho0);

        return new BenchmarkResult(logL, logLikNanos, combinedNanos);
    }

    public static Test suite() {
        return new TestSuite(StanVsBeastOUBenchmarkTest.class);
    }

    public static void main(final String[] args) {
        junit.textui.TestRunner.run(suite());
    }

    // -----------------------------------------------------------------------
    // Model construction
    // -----------------------------------------------------------------------

    /**
     * Builds an OU process model for even K using the orthogonal block-diagonal polar-stable
     * parametrization.
     *
     * For even K, the block structure is K/2 two-by-two blocks, no leading 1x1 scalar block.
     *   - Givens rotation angles: K*(K-1)/2
     *   - rho, theta, t: each K/2 values (one per 2x2 block)
     *   - scalarBlock: dimension 0 (not used for even K)
     */
    private static OUProcessModel makeOrthogonalBlockProcess(final String name) {
        if (K % 2 != 0) {
            throw new IllegalStateException("This benchmark requires even K; got K=" + K);
        }
        final int nAngles = K * (K - 1) / 2;
        final int nBlocks = K / 2;

        // Givens rotation angles: deterministic sequence
        final double[] anglesData = new double[nAngles];
        for (int i = 0; i < nAngles; ++i) {
            anglesData[i] = 0.15 * Math.sin(0.7 * i + 0.3);
        }
        final Parameter angles = new Parameter.Default(name + ".angles", anglesData);
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter(name + ".rotation", K, angles);

        // Scalar block: empty for even K
        final Parameter scalar = new Parameter.Default(0);

        // 2x2-block native parameters: deterministic sequences
        final double[] rhoData   = new double[nBlocks];
        final double[] thetaData = new double[nBlocks];
        final double[] tData     = new double[nBlocks];
        for (int b = 0; b < nBlocks; ++b) {
            rhoData[b]   = 0.60 + 0.15 * Math.sin(1.1 * b);
            thetaData[b] = 0.25 + 0.10 * Math.cos(0.9 * b + 0.5);
            tData[b]     = 0.05 * Math.sin(1.3 * b + 1.0);
        }
        final Parameter rho   = new Parameter.Default(name + ".rho",   rhoData);
        final Parameter theta = new Parameter.Default(name + ".theta", thetaData);
        final Parameter t     = new Parameter.Default(name + ".t",     tData);

        final OrthogonalBlockDiagonalPolarStableMatrixParameter drift =
                new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                        name + ".A", rotation, scalar, rho, theta, t);

        // Diffusion matrix Q and initial covariance P0: K x K SPD, diagonally dominant
        final MatrixParameter Q  = makeSpdMatrix(name + ".Q",  1.0, 0.12);
        final MatrixParameter P0 = makeSpdMatrix(name + ".P0", 1.0, 0.12);

        // Stationary mean: deterministic
        final double[] muData = new double[K];
        for (int i = 0; i < K; ++i) {
            muData[i] = 0.10 * Math.sin(0.5 * i + 0.2);
        }
        final Parameter mu = new Parameter.Default(name + ".mu", muData);

        return new OUProcessModel(name, K, drift, Q, mu, P0);
    }

    /**
     * Creates a K x K symmetric positive-definite matrix with diagonal value {@code diagVal}
     * and off-diagonal entries {@code offScale * sin(i+j+1)} for reproducibility.
     * The off-diagonal entries are small enough to ensure positive definiteness.
     */
    private static MatrixParameter makeSpdMatrix(final String name,
                                                 final double diagVal,
                                                 final double offScale) {
        final double[][] data = new double[K][K];
        for (int i = 0; i < K; ++i) {
            data[i][i] = diagVal;
            for (int j = i + 1; j < K; ++j) {
                final double v = offScale * Math.sin(i + j + 1.0);
                data[i][j] = v;
                data[j][i] = v;
            }
        }
        return makeMatrix(name, data);
    }

    /**
     * Creates a Gaussian observation model with identity design matrix H=I_K and
     * isotropic noise covariance R = sigma_obs^2 * I_K.
     * Observations are a deterministic sinusoidal pattern (K x T_OBS) for reproducibility.
     */
    private static GaussianObservationModel makeObservation(final String name) {
        // H = I_K
        final double[][] hData = new double[K][K];
        for (int i = 0; i < K; ++i) {
            hData[i][i] = 1.0;
        }

        // R = sigma_obs^2 * I_K
        final double[][] rData = new double[K][K];
        for (int i = 0; i < K; ++i) {
            rData[i][i] = SIGMA_OBS * SIGMA_OBS;
        }

        // Observation matrix: K rows, T_OBS columns – deterministic sinusoidal pattern
        final double[][] y = new double[K][T_OBS];
        for (int k = 0; k < K; ++k) {
            for (int t = 0; t < T_OBS; ++t) {
                y[k][t] = 0.25 * Math.sin(0.08 * t + 0.5 * k)
                        + 0.06 * Math.cos(0.21 * t - 0.3 * k);
            }
        }

        return new GaussianObservationModel(name, K,
                makeMatrix(name + ".H", hData),
                makeMatrix(name + ".R", rData),
                makeMatrix(name + ".Y", y));
    }

    private static MatrixParameter makeMatrix(final String name, final double[][] values) {
        final MatrixParameter m = new MatrixParameter(name, values.length, values[0].length);
        for (int i = 0; i < values.length; ++i) {
            for (int j = 0; j < values[0].length; ++j) {
                m.setParameterValue(i, j, values[i][j]);
            }
        }
        return m;
    }

    // -----------------------------------------------------------------------
    // Result container
    // -----------------------------------------------------------------------

    private static final class BenchmarkResult {
        final double logL;
        final long[] logLikNanos;
        final long[] combinedNanos;

        BenchmarkResult(final double logL, final long[] logLikNanos, final long[] combinedNanos) {
            this.logL = logL;
            this.logLikNanos = logLikNanos;
            this.combinedNanos = combinedNanos;
        }
    }

    private static void printRow(final String label,
                                 final long[] uniformNanos,
                                 final long[] irregularNanos) {
        final double u = medianMs(uniformNanos);
        final double ir = medianMs(irregularNanos);
        System.out.printf("%-32s  %12.3f  %12.3f  %9.1fx%n", label, u, ir, ir / u);
    }

    // -----------------------------------------------------------------------
    // Timing helpers
    // -----------------------------------------------------------------------

    private static double medianMs(final long[] nanos) {
        final long[] copy = Arrays.copyOf(nanos, nanos.length);
        Arrays.sort(copy);
        final int n = copy.length;
        final double medNs = (n % 2 == 0)
                ? 0.5 * (copy[n / 2 - 1] + copy[n / 2])
                : copy[n / 2];
        return medNs / 1.0e6;
    }

    private static double meanMs(final long[] nanos) {
        long sum = 0L;
        for (final long v : nanos) sum += v;
        return sum / (double) nanos.length / 1.0e6;
    }

    private static double minMs(final long[] nanos) {
        long min = Long.MAX_VALUE;
        for (final long v : nanos) if (v < min) min = v;
        return min / 1.0e6;
    }

    private static double maxMs(final long[] nanos) {
        long max = Long.MIN_VALUE;
        for (final long v : nanos) if (v > max) max = v;
        return max / 1.0e6;
    }

    private static String repeat(final char value, final int count) {
        final char[] chars = new char[count];
        Arrays.fill(chars, value);
        return new String(chars);
    }
}
