package dr.inference.operators;

import dr.math.MathUtils;

final class UnivariateSliceSampler {

    interface LogDensity {
        double logDensity(double x);
    }

    private UnivariateSliceSampler() {
        // no instances
    }

    static double sample(final double current,
                         final double lower,
                         final double upper,
                         final double windowSize,
                         final int maxSteppingOut,
                         final int maxShrinkIterations,
                         final LogDensity logDensity) {
        final double currentLogDensity = logDensity.logDensity(current);
        if (!Double.isFinite(currentLogDensity)) {
            return Double.NaN;
        }

        final double cutoff = currentLogDensity + MathUtils.randomLogDouble();
        final SliceInterval interval = constructInterval(
                current,
                cutoff,
                lower,
                upper,
                windowSize,
                maxSteppingOut,
                logDensity);
        return drawFromInterval(
                current,
                cutoff,
                interval,
                maxShrinkIterations,
                logDensity);
    }

    private static SliceInterval constructInterval(final double current,
                                                   final double cutoff,
                                                   final double lower,
                                                   final double upper,
                                                   final double windowSize,
                                                   final int maxSteppingOut,
                                                   final LogDensity logDensity) {
        double left = current - windowSize * MathUtils.nextDouble();
        double right = left + windowSize;
        if (left < lower) {
            left = lower;
        }
        if (right > upper) {
            right = upper;
        }

        int leftSteps = MathUtils.nextInt(maxSteppingOut);
        int rightSteps = (maxSteppingOut - 1) - leftSteps;

        while (leftSteps > 0 && left > lower) {
            final double nextLeft = Math.max(lower, left - windowSize);
            if (!(logDensity.logDensity(nextLeft) > cutoff)) {
                break;
            }
            left = nextLeft;
            leftSteps--;
        }

        while (rightSteps > 0 && right < upper) {
            final double nextRight = Math.min(upper, right + windowSize);
            if (!(logDensity.logDensity(nextRight) > cutoff)) {
                break;
            }
            right = nextRight;
            rightSteps--;
        }

        return new SliceInterval(left, right);
    }

    private static double drawFromInterval(final double current,
                                           final double cutoff,
                                           final SliceInterval interval,
                                           final int maxShrinkIterations,
                                           final LogDensity logDensity) {
        double left = interval.left;
        double right = interval.right;
        for (int i = 0; i < maxShrinkIterations; i++) {
            final double proposed = MathUtils.uniform(left, right);
            if (logDensity.logDensity(proposed) >= cutoff) {
                return proposed;
            }
            if (proposed < current) {
                left = proposed;
            } else {
                right = proposed;
            }
        }
        return Double.NaN;
    }

    private static final class SliceInterval {
        private final double left;
        private final double right;

        private SliceInterval(final double left, final double right) {
            this.left = left;
            this.right = right;
        }
    }
}
