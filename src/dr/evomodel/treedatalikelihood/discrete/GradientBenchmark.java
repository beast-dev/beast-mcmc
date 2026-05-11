package dr.evomodel.treedatalikelihood.discrete;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.xml.Reportable;

import java.util.Arrays;
import java.util.logging.Logger;

/**
 * Benchmarks the wall-clock cost of repeated gradient evaluations and reports
 * mean / std-dev / min / max per-call time.
 *
 * When {@code dirtyTarget} is provided and {@code dirty=true}, calls
 * {@link Likelihood#makeDirty()} before each iteration to simulate the realistic
 * HMC cost: post-order traversal + pre-order traversal + Fréchet kernel.
 * Without dirty, only the Fréchet kernel is timed (caches are warm).
 *
 * @author Filippo Monti
 */
public class GradientBenchmark implements Reportable {

    private final GradientWrtParameterProvider gradient;
    private final Likelihood dirtyTarget;
    private final int warmup;
    private final int iterations;
    private final boolean dirty;

    public GradientBenchmark(GradientWrtParameterProvider gradient, int warmup, int iterations,
                             boolean dirty, Likelihood dirtyTarget) {
        this.gradient    = gradient;
        this.warmup      = warmup;
        this.iterations  = iterations;
        this.dirty       = dirty;
        this.dirtyTarget = dirtyTarget;
    }

    @Override
    public String getReport() {
        Logger.getLogger("dr.evomodel").info(
                "GradientBenchmark: running " + warmup + " warmup + " + iterations +
                " timed calls" + (dirty ? " [dirty=true: includes pre/post-order traversal]" : "") + "…");

        for (int i = 0; i < warmup; i++) {
            if (dirty && dirtyTarget != null) dirtyTarget.makeDirty();
            gradient.getGradientLogDensity();
        }

        final long[] nanos = new long[iterations];
        for (int i = 0; i < iterations; i++) {
            if (dirty && dirtyTarget != null) dirtyTarget.makeDirty();
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
                "GradientBenchmark (%s, dirty=%b):%n" +
                "  warmup=%d  iterations=%d%n" +
                "  mean=%.4f ms  sd=%.4f ms%n" +
                "  min=%.4f ms  p50=%.4f ms  p95=%.4f ms  max=%.4f ms",
                gradient.getClass().getSimpleName(), dirty,
                warmup, iterations,
                mean, stddev,
                min, p50, p95, max);
    }
}
