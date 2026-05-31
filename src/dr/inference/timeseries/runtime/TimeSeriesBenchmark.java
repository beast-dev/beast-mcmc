package dr.inference.timeseries.runtime;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

import java.util.Arrays;
import java.util.Locale;

/**
 * Small runnable benchmark for time-series likelihood and gradient evaluation.
 */
public final class TimeSeriesBenchmark implements Runnable {

    public enum Mode {
        LOG_LIKELIHOOD("logLikelihood"),
        GRADIENT("gradient");

        private final String label;

        Mode(final String label) {
            this.label = label;
        }

        public String getLabel() {
            return label;
        }

        public static Mode parse(final String value) {
            if (value == null) {
                throw new IllegalArgumentException("mode must not be null");
            }
            if ("logLikelihood".equalsIgnoreCase(value)
                    || "logL".equalsIgnoreCase(value)
                    || "likelihood".equalsIgnoreCase(value)) {
                return LOG_LIKELIHOOD;
            }
            if ("gradient".equalsIgnoreCase(value)
                    || "grad".equalsIgnoreCase(value)) {
                return GRADIENT;
            }
            throw new IllegalArgumentException("Unknown time-series benchmark mode: " + value);
        }
    }

    private final String id;
    private final Mode mode;
    private final Likelihood likelihood;
    private final GradientWrtParameterProvider gradient;
    private final int warmupIterations;
    private final int timedIterations;
    private final boolean perturbEachIteration;
    private final Parameter perturbParameter;
    private final int perturbIndex;
    private final double perturbation;
    private final String description;
    private double perturbBaseValue;

    public TimeSeriesBenchmark(final String id,
                               final Mode mode,
                               final Likelihood likelihood,
                               final GradientWrtParameterProvider gradient,
                               final int warmupIterations,
                               final int timedIterations,
                               final boolean perturbEachIteration,
                               final Parameter perturbParameter,
                               final int perturbIndex,
                               final double perturbation,
                               final String description) {
        if (mode == null) {
            throw new IllegalArgumentException("mode must not be null");
        }
        if (likelihood == null) {
            throw new IllegalArgumentException("likelihood must not be null");
        }
        if (mode == Mode.GRADIENT && gradient == null) {
            throw new IllegalArgumentException("gradient mode requires a gradient provider");
        }
        if (warmupIterations < 0) {
            throw new IllegalArgumentException("warmupIterations must be non-negative");
        }
        if (timedIterations < 1) {
            throw new IllegalArgumentException("timedIterations must be positive");
        }
        if (perturbEachIteration) {
            if (perturbParameter == null) {
                throw new IllegalArgumentException("perturbEachIteration requires a perturb parameter");
            }
            if (perturbIndex < 0 || perturbIndex >= perturbParameter.getDimension()) {
                throw new IllegalArgumentException("perturbIndex out of range for perturb parameter");
            }
            if (!(perturbation > 0.0)) {
                throw new IllegalArgumentException("perturbation must be positive");
            }
        }
        this.id = id == null ? "timeSeriesBenchmark" : id;
        this.mode = mode;
        this.likelihood = likelihood;
        this.gradient = gradient;
        this.warmupIterations = warmupIterations;
        this.timedIterations = timedIterations;
        this.perturbEachIteration = perturbEachIteration;
        this.perturbParameter = perturbParameter;
        this.perturbIndex = perturbIndex;
        this.perturbation = perturbation;
        this.description = description == null ? "" : description;
    }

    @Override
    public void run() {
        final double originalPerturbValue = perturbEachIteration
                ? perturbParameter.getParameterValue(perturbIndex)
                : Double.NaN;
        perturbBaseValue = originalPerturbValue;
        try {
            printHeader();
            if (mode == Mode.LOG_LIKELIHOOD) {
                runLogLikelihoodCase();
            } else {
                runGradientCase();
            }
        } finally {
            if (perturbEachIteration) {
                perturbParameter.setParameterValue(perturbIndex, originalPerturbValue);
                likelihood.makeDirty();
                if (gradient != null && gradient.getLikelihood() != likelihood) {
                    gradient.getLikelihood().makeDirty();
                }
            }
        }
    }

    private void printHeader() {
        System.out.println();
        System.out.println("=== Time-series benchmark: " + id + " ===");
        if (!description.isEmpty()) {
            System.out.println(description);
        }
        System.out.printf(Locale.US,
                "mode=%s  warmup=%d  timedIterations=%d  perturbEachIteration=%s%n",
                mode.getLabel(), warmupIterations, timedIterations, perturbEachIteration);
        if (perturbEachIteration) {
            System.out.printf(Locale.US,
                    "perturbParameter=%s[%d]  perturbation=%.4g%n",
                    perturbParameter.getParameterName(), perturbIndex, perturbation);
        }
    }

    private void runLogLikelihoodCase() {
        double checksum = 0.0;
        for (int i = 0; i < warmupIterations; ++i) {
            prepareIteration(i);
            final double logLikelihood = likelihood.getLogLikelihood();
            requireFinite(logLikelihood, "warmup log-likelihood");
            checksum += logLikelihood;
        }

        final long[] nanos = new long[timedIterations];
        double lastLogLikelihood = Double.NaN;
        for (int i = 0; i < timedIterations; ++i) {
            prepareIteration(warmupIterations + i);
            final long start = System.nanoTime();
            lastLogLikelihood = likelihood.getLogLikelihood();
            nanos[i] = System.nanoTime() - start;
            requireFinite(lastLogLikelihood, "timed log-likelihood");
            checksum += lastLogLikelihood;
        }

        System.out.printf(Locale.US, "lastLogLikelihood=%.10f  checksum=%.10f%n",
                lastLogLikelihood, checksum);
        printTiming("logLikelihood", nanos);
    }

    private void runGradientCase() {
        double checksum = 0.0;
        for (int i = 0; i < warmupIterations; ++i) {
            prepareIteration(i);
            final double[] values = gradient.getGradientLogDensity();
            checksum += checksum(values);
        }

        final long[] nanos = new long[timedIterations];
        double[] lastGradient = null;
        for (int i = 0; i < timedIterations; ++i) {
            prepareIteration(warmupIterations + i);
            final long start = System.nanoTime();
            lastGradient = gradient.getGradientLogDensity();
            nanos[i] = System.nanoTime() - start;
            checksum += checksum(lastGradient);
        }

        System.out.printf(Locale.US, "gradientDimension=%d  gradientChecksum=%.10f%n",
                lastGradient == null ? gradient.getDimension() : lastGradient.length, checksum);
        printTiming("gradient", nanos);
    }

    private void prepareIteration(final int iteration) {
        if (perturbEachIteration) {
            final double direction = (iteration & 1) == 0 ? 1.0 : -1.0;
            perturbParameter.setParameterValue(perturbIndex, perturbBaseValue + direction * perturbation);
        }
        likelihood.makeDirty();
        if (gradient != null && gradient.getLikelihood() != likelihood) {
            gradient.getLikelihood().makeDirty();
        }
    }

    private static void printTiming(final String label, final long[] nanos) {
        System.out.printf(Locale.US,
                "%s medianMs=%.6f  meanMs=%.6f  minMs=%.6f  maxMs=%.6f%n",
                label, medianMs(nanos), meanMs(nanos), minMs(nanos), maxMs(nanos));
    }

    private static double checksum(final double[] values) {
        double sum = 0.0;
        for (int i = 0; i < values.length; ++i) {
            requireFinite(values[i], "gradient[" + i + "]");
            sum += (i + 1.0) * values[i];
        }
        return sum;
    }

    private static void requireFinite(final double value, final String name) {
        if (!Double.isFinite(value)) {
            throw new IllegalStateException(name + " is not finite: " + value);
        }
    }

    private static double medianMs(final long[] nanos) {
        final long[] copy = Arrays.copyOf(nanos, nanos.length);
        Arrays.sort(copy);
        final int n = copy.length;
        final double median = (n & 1) == 0
                ? 0.5 * (copy[n / 2 - 1] + copy[n / 2])
                : copy[n / 2];
        return median / 1.0e6;
    }

    private static double meanMs(final long[] nanos) {
        long sum = 0L;
        for (final long value : nanos) {
            sum += value;
        }
        return sum / (double) nanos.length / 1.0e6;
    }

    private static double minMs(final long[] nanos) {
        long min = Long.MAX_VALUE;
        for (final long value : nanos) {
            if (value < min) {
                min = value;
            }
        }
        return min / 1.0e6;
    }

    private static double maxMs(final long[] nanos) {
        long max = Long.MIN_VALUE;
        for (final long value : nanos) {
            if (value > max) {
                max = value;
            }
        }
        return max / 1.0e6;
    }
}
