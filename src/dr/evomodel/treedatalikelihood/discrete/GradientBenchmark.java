package dr.evomodel.treedatalikelihood.discrete;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.xml.Reportable;

import java.util.Arrays;
import java.util.logging.Logger;

/**
 * Benchmarks the wall-clock cost of repeated gradient evaluations and reports
 * mean / std-dev / min / max per-call time.
 *
 * Usage in a BEAST XML file:
 * <pre>{@code
 * <report>
 *   <benchmarkGradient warmup="50" iterations="200">
 *     <exactLogCtmcRateGradient .../>
 *   </benchmarkGradient>
 * </report>
 * }</pre>
 *
 * The benchmark is lazy: it runs only when {@link #getReport()} is first called,
 * so it executes after all model setup is complete (same call site as a normal
 * BEAST {@code <report>} block).
 *
 * @author Filippo Monti
 */
public class GradientBenchmark implements Reportable {

    private final GradientWrtParameterProvider gradient;
    private final int warmup;
    private final int iterations;

    public GradientBenchmark(GradientWrtParameterProvider gradient, int warmup, int iterations) {
        this.gradient   = gradient;
        this.warmup     = warmup;
        this.iterations = iterations;
    }

    @Override
    public String getReport() {
        Logger.getLogger("dr.evomodel").info(
                "GradientBenchmark: running " + warmup + " warmup + " + iterations + " timed calls…");

        for (int i = 0; i < warmup; i++) {
            gradient.getGradientLogDensity();
        }

        final long[] nanos = new long[iterations];
        for (int i = 0; i < iterations; i++) {
            final long t0 = System.nanoTime();
            gradient.getGradientLogDensity();
            nanos[i] = System.nanoTime() - t0;
        }

        final double[] ms = new double[iterations];
        for (int i = 0; i < iterations; i++) ms[i] = nanos[i] * 1e-6;

        Arrays.sort(ms);

        double sum = 0.0, sum2 = 0.0;
        for (double v : ms) { sum += v; sum2 += v * v; }
        final double mean   = sum / iterations;
        final double stddev = Math.sqrt(Math.max(0.0, sum2 / iterations - mean * mean));
        final double min    = ms[0];
        final double max    = ms[iterations - 1];
        final double p50    = ms[iterations / 2];
        final double p95    = ms[(int) (iterations * 0.95)];

        return String.format(
                "GradientBenchmark (%s):%n" +
                "  warmup=%d  iterations=%d%n" +
                "  mean=%.4f ms  sd=%.4f ms%n" +
                "  min=%.4f ms  p50=%.4f ms  p95=%.4f ms  max=%.4f ms",
                gradient.getClass().getSimpleName(),
                warmup, iterations,
                mean, stddev,
                min, p50, p95, max);
    }
}
