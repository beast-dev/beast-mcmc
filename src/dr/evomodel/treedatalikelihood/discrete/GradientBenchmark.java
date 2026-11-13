package dr.evomodel.treedatalikelihood.discrete;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.xml.Reportable;

import java.io.File;
import java.lang.reflect.Method;
import java.time.Duration;
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

    private static final Logger LOGGER = Logger.getLogger("dr.evomodel");
    private static final String JFR_FILE_PROPERTY = "beast.gradient.benchmark.jfr";
    private static final String JFR_DIRTY_ONLY_PROPERTY = "beast.gradient.benchmark.jfr.dirtyOnly";
    private static final String JFR_PERIOD_MILLIS_PROPERTY = "beast.gradient.benchmark.jfr.periodMillis";
    private static final long DEFAULT_JFR_PERIOD_MILLIS = 10L;

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
        LOGGER.info(
                "GradientBenchmark: running " + warmup + " warmup + " + iterations +
                " timed calls" + (dirty ? " [dirty=true: includes pre/post-order traversal]" : "") + "…");

        for (int i = 0; i < warmup; i++) {
            if (dirty && dirtyTarget != null) dirtyTarget.makeDirty();
            gradient.getGradientLogDensity();
        }

        final long[] nanos = new long[iterations];
        final TimedLoopJfrRecording jfrRecording = startTimedLoopJfrRecording();
        try {
            for (int i = 0; i < iterations; i++) {
                if (dirty && dirtyTarget != null) dirtyTarget.makeDirty();
                final long t0 = System.nanoTime();
                gradient.getGradientLogDensity();
                nanos[i] = System.nanoTime() - t0;
            }
        } finally {
            stopTimedLoopJfrRecording(jfrRecording);
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

    private TimedLoopJfrRecording startTimedLoopJfrRecording() {
        final String fileName = System.getProperty(JFR_FILE_PROPERTY);
        if (fileName == null || fileName.trim().isEmpty()) {
            return null;
        }
        final boolean dirtyOnly = Boolean.parseBoolean(
                System.getProperty(JFR_DIRTY_ONLY_PROPERTY, "true"));
        if (dirtyOnly && !dirty) {
            return null;
        }

        try {
            final Class<?> recordingClass = Class.forName("jdk.jfr.Recording");
            final Object recording = recordingClass.getConstructor().newInstance();
            final File file = resolveJfrFile(fileName);
            final File parent = file.getParentFile();
            if (parent != null && !parent.exists() && !parent.mkdirs()) {
                LOGGER.warning("Could not create JFR output directory: " + parent);
                return null;
            }

            recordingClass.getMethod("setName", String.class).invoke(
                    recording, "GradientBenchmark dirty=" + dirty);
            recordingClass.getMethod("setToDisk", boolean.class).invoke(recording, true);

            final long periodMillis = Long.parseLong(
                    System.getProperty(JFR_PERIOD_MILLIS_PROPERTY,
                            Long.toString(DEFAULT_JFR_PERIOD_MILLIS)));
            enableJfrEvent(recordingClass, recording, "jdk.ExecutionSample", periodMillis);
            enableJfrEvent(recordingClass, recording, "jdk.ObjectAllocationSample", periodMillis);

            recordingClass.getMethod("start").invoke(recording);
            LOGGER.info("GradientBenchmark: JFR recording started after warmup: " + file);
            return new TimedLoopJfrRecording(recordingClass, recording, file);
        } catch (Throwable t) {
            LOGGER.warning("GradientBenchmark: could not start timed-loop JFR recording: "
                    + t.getMessage());
            return null;
        }
    }

    private File resolveJfrFile(String fileName) {
        String resolved = fileName
                .replace("{dirty}", Boolean.toString(dirty))
                .replace("{gradient}", gradient.getClass().getSimpleName());
        if (resolved.equals(fileName)
                && !Boolean.parseBoolean(System.getProperty(JFR_DIRTY_ONLY_PROPERTY, "true"))) {
            final int dot = resolved.lastIndexOf('.');
            final String suffix = "-dirty-" + dirty;
            resolved = dot > 0
                    ? resolved.substring(0, dot) + suffix + resolved.substring(dot)
                    : resolved + suffix;
        }
        return new File(resolved);
    }

    private static void enableJfrEvent(Class<?> recordingClass, Object recording,
                                       String eventName, long periodMillis) {
        try {
            final Object settings = recordingClass.getMethod("enable", String.class)
                    .invoke(recording, eventName);
            final Method withPeriod = settings.getClass().getMethod("withPeriod", Duration.class);
            withPeriod.invoke(settings, Duration.ofMillis(periodMillis));
            final Method withStackTrace = settings.getClass().getMethod("withStackTrace", boolean.class);
            withStackTrace.invoke(settings, true);
        } catch (Throwable ignored) {
            // Some JVMs do not expose every event or setting. Keep the benchmark runnable.
        }
    }

    private static void stopTimedLoopJfrRecording(TimedLoopJfrRecording recording) {
        if (recording == null) {
            return;
        }
        try {
            recording.recordingClass.getMethod("stop").invoke(recording.recording);
            recording.recordingClass.getMethod("dump", java.nio.file.Path.class)
                    .invoke(recording.recording, recording.file.toPath());
            LOGGER.info("GradientBenchmark: JFR recording written to " + recording.file);
        } catch (Throwable t) {
            LOGGER.warning("GradientBenchmark: could not write timed-loop JFR recording: "
                    + t.getMessage());
        } finally {
            try {
                recording.recordingClass.getMethod("close").invoke(recording.recording);
            } catch (Throwable ignored) {
                // Nothing useful to do at shutdown.
            }
        }
    }

    private static final class TimedLoopJfrRecording {
        private final Class<?> recordingClass;
        private final Object recording;
        private final File file;

        private TimedLoopJfrRecording(Class<?> recordingClass, Object recording, File file) {
            this.recordingClass = recordingClass;
            this.recording = recording;
            this.file = file;
        }
    }
}
