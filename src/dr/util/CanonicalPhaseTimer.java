package dr.util;

import java.util.Locale;

/**
 * Debug-only aggregate timing for canonical OU gradient phases.
 */
public final class CanonicalPhaseTimer {

    public static final String PROPERTY = "beast.debug.canonicalPhaseTiming";
    public static final String REPORT_EVERY_PROPERTY = "beast.debug.canonicalPhaseTimingReportEvery";

    private static final boolean ENABLED = Boolean.getBoolean(PROPERTY);
    private static final int REPORT_EVERY = Math.max(1, Integer.getInteger(REPORT_EVERY_PROPERTY, 50));

    private static long postorderNanos;
    private static long preorderNanos;
    private static long branchPrepNanos;
    private static long branchGradientNanos;
    private static long reductionNanos;
    private static long finalizeNanos;
    private static long totalNanos;
    private static long evaluations;
    private static int lastActiveBranches;
    private static int lastWorkerCount;
    private static boolean lastParallel;

    private CanonicalPhaseTimer() { }

    public static boolean isEnabled() {
        return ENABLED;
    }

    public static long start() {
        return ENABLED ? System.nanoTime() : 0L;
    }

    public static void recordPostorder(final long startNanos) {
        if (ENABLED) {
            postorderNanos += elapsedSince(startNanos);
        }
    }

    public static void recordPreorder(final long startNanos) {
        if (ENABLED) {
            preorderNanos += elapsedSince(startNanos);
        }
    }

    public static void recordBranchPrep(final long startNanos) {
        if (ENABLED) {
            branchPrepNanos += elapsedSince(startNanos);
        }
    }

    public static void recordBranchGradient(final long startNanos) {
        if (ENABLED) {
            branchGradientNanos += elapsedSince(startNanos);
        }
    }

    public static void recordReduction(final long startNanos) {
        if (ENABLED) {
            reductionNanos += elapsedSince(startNanos);
        }
    }

    public static void recordFinalize(final long startNanos) {
        if (ENABLED) {
            finalizeNanos += elapsedSince(startNanos);
        }
    }

    public static void recordGradientShape(final int activeBranches,
                                           final int workerCount,
                                           final boolean parallel) {
        if (ENABLED) {
            lastActiveBranches = activeBranches;
            lastWorkerCount = workerCount;
            lastParallel = parallel;
        }
    }

    public static void finishJointGradient(final long startNanos) {
        if (!ENABLED) {
            return;
        }
        totalNanos += elapsedSince(startNanos);
        evaluations++;
        if (evaluations % REPORT_EVERY == 0) {
            reportAndReset();
        }
    }

    private static long elapsedSince(final long startNanos) {
        return System.nanoTime() - startNanos;
    }

    private static void reportAndReset() {
        final double totalMs = toMillis(totalNanos);
        final String message = String.format(Locale.US,
                "[canonical-phase-timing] count=%d activeBranches=%d workers=%d parallel=%s "
                        + "total=%.3fms postorder=%.3fms preorder=%.3fms branchPrep=%.3fms "
                        + "branchGradient=%.3fms reduction=%.3fms finalize=%.3fms "
                        + "branchGradientShare=%.1f%%",
                evaluations,
                lastActiveBranches,
                lastWorkerCount,
                Boolean.toString(lastParallel),
                totalMs,
                toMillis(postorderNanos),
                toMillis(preorderNanos),
                toMillis(branchPrepNanos),
                toMillis(branchGradientNanos),
                toMillis(reductionNanos),
                toMillis(finalizeNanos),
                totalNanos == 0L ? 0.0 : 100.0 * branchGradientNanos / totalNanos);
        System.err.println(message);
        reset();
    }

    private static double toMillis(final long nanos) {
        return nanos / 1.0e6;
    }

    private static void reset() {
        postorderNanos = 0L;
        preorderNanos = 0L;
        branchPrepNanos = 0L;
        branchGradientNanos = 0L;
        reductionNanos = 0L;
        finalizeNanos = 0L;
        totalNanos = 0L;
        evaluations = 0L;
    }
}
