/*
 * SericolaSeriesWeights.java
 *
 * Stable probability weights for the Sericola reward-density series.
 */

package dr.inference.markovjumps;

import dr.math.GammaFunction;

import java.util.Arrays;

final class SericolaSeriesWeights {

    private SericolaSeriesWeights() {
        // Utility class.
    }

    static int determinePoissonUpperTruncation(double mean, double epsilon) {
        return determinePoissonUpperTruncation(mean, epsilon, true);
    }

    static int determinePoissonUpperTruncation(double mean, double epsilon, boolean rescale) {
        validateMean(mean);
        if (!(epsilon > 0.0 && epsilon < 1.0)) {
            throw new IllegalArgumentException("epsilon must be in (0,1)");
        }
        if (!rescale) {
            return determinePoissonUpperTruncationUnscaled(mean, epsilon);
        }
        if (mean == 0.0) {
            return 0;
        }

        final int mode = poissonMode(mean);
        double lowerAndModeScaledMass = 1.0;
        double term = 1.0;
        for (int n = mode; n >= 1; --n) {
            term *= n / mean;
            if (term == 0.0) {
                break;
            }
            lowerAndModeScaledMass += term;
        }

        int n = mode;
        term = 1.0;
        final int hardCap = poissonHardCap(mean);

        while (true) {
            final double next = n + 1.0;
            final double ratio = mean / next;
            if (ratio < 1.0) {
                final double nextTerm = term * ratio;
                final double tailBound = nextTerm / (1.0 - ratio);
                final double boundedTotal = lowerAndModeScaledMass + tailBound;
                if (tailBound <= epsilon * boundedTotal) {
                    return n;
                }
            }

            n++;
            if (n > hardCap) {
                throw new IllegalStateException(
                        "Poisson truncation did not converge before hard cap; mean=" +
                                mean + " epsilon=" + epsilon + " hardCap=" + hardCap);
            }

            term *= mean / n;
            if (term == 0.0) {
                return n;
            }
            lowerAndModeScaledMass += term;
        }
    }

    static void fillPoissonWeights(double mean, int maxN, double[] out, int offset) {
        fillPoissonWeights(mean, maxN, out, offset, true);
    }

    static void fillPoissonWeights(double mean, int maxN, double[] out, int offset, boolean rescale) {
        validateMean(mean);
        if (maxN < 0) {
            throw new IllegalArgumentException("maxN must be >= 0");
        }
        if (out == null || offset < 0 || offset + maxN >= out.length) {
            throw new IllegalArgumentException("Output array too small for Poisson weights");
        }
        if (!rescale) {
            fillPoissonWeightsUnscaled(mean, maxN, out, offset);
            return;
        }

        Arrays.fill(out, offset, offset + maxN + 1, 0.0);

        if (mean == 0.0) {
            out[offset] = 1.0;
            return;
        }

        final int mode = Math.min(poissonMode(mean), maxN);
        final double modeWeight = Math.exp(logPoissonProbability(mean, mode));
        if (!(modeWeight > 0.0) || Double.isNaN(modeWeight) || Double.isInfinite(modeWeight)) {
            throw new IllegalStateException(
                    "Could not compute finite Poisson mode weight; mean=" + mean + " mode=" + mode);
        }

        out[offset + mode] = modeWeight;

        double weight = modeWeight;
        for (int n = mode; n >= 1; --n) {
            weight *= n / mean;
            if (weight == 0.0) {
                break;
            }
            out[offset + n - 1] = weight;
        }

        weight = modeWeight;
        for (int n = mode; n < maxN; ++n) {
            weight *= mean / (n + 1.0);
            if (weight == 0.0) {
                break;
            }
            out[offset + n + 1] = weight;
        }
    }

    static void fillBernsteinWeights(int degree, double x, double[] out) {
        fillBernsteinWeights(degree, x, out, true);
    }

    // Tolerance for ordinary floating-point noise around the [0,1] boundary
    // (e.g. an unconstrained continuous HMC coordinate landing at
    // 1.0000000000000002). Values outside [0,1] by more than this, or NaN,
    // indicate a genuinely invalid reward-mixture state -- see
    // RewardDensityDomainException.
    private static final double BOUNDARY_TOLERANCE = 1.0e-9;

    static void fillBernsteinWeights(int degree, double x, double[] out, boolean rescale) {
        if (degree < 0) {
            throw new IllegalArgumentException("degree must be >= 0");
        }
        if (out == null || out.length < degree + 1) {
            throw new IllegalArgumentException("Output array too small for Bernstein weights");
        }
        if (!rescale) {
            fillBernsteinWeightsUnscaled(degree, x, out);
            return;
        }

        Arrays.fill(out, 0, degree + 1, 0.0);

        if (degree == 0) {
            out[0] = 1.0;
            return;
        }
        if (Double.isNaN(x) || x < -BOUNDARY_TOLERANCE || x > 1.0 + BOUNDARY_TOLERANCE) {
            throw new RewardDensityDomainException(
                    "Reward value " + x + " is not a valid reward-mixture state (must be in [0,1])");
        }
        if (x < 0.0) {
            x = 0.0;
        } else if (x > 1.0) {
            x = 1.0;
        }
        if (x == 0.0) {
            out[0] = 1.0;
            return;
        }
        if (x == 1.0) {
            out[degree] = 1.0;
            return;
        }

        int mode = (int) Math.floor((degree + 1.0) * x);
        if (mode > degree) {
            mode = degree;
        }

        final double logModeWeight =
                GammaFunction.lnGamma(degree + 1.0)
                        - GammaFunction.lnGamma(mode + 1.0)
                        - GammaFunction.lnGamma(degree - mode + 1.0)
                        + mode * Math.log(x)
                        + (degree - mode) * Math.log1p(-x);

        final double modeWeight = Math.exp(logModeWeight);
        if (!(modeWeight > 0.0) || Double.isNaN(modeWeight) || Double.isInfinite(modeWeight)) {
            throw new IllegalStateException(
                    "Could not compute finite Bernstein mode weight; degree=" +
                            degree + " x=" + x + " mode=" + mode);
        }

        out[mode] = modeWeight;

        double weight = modeWeight;
        final double leftRatio = (1.0 - x) / x;
        for (int k = mode; k >= 1; --k) {
            weight *= (k / (degree - k + 1.0)) * leftRatio;
            if (weight == 0.0) {
                break;
            }
            out[k - 1] = weight;
        }

        weight = modeWeight;
        final double rightRatio = x / (1.0 - x);
        for (int k = mode; k < degree; ++k) {
            weight *= ((degree - k) / (k + 1.0)) * rightRatio;
            if (weight == 0.0) {
                break;
            }
            out[k + 1] = weight;
        }

        double sum = 0.0;
        for (int k = 0; k <= degree; ++k) {
            sum += out[k];
        }
        if (!(sum > 0.0) || Double.isNaN(sum) || Double.isInfinite(sum)) {
            throw new IllegalStateException(
                    "Could not normalize Bernstein weights; degree=" + degree + " x=" + x);
        }
        final double invSum = 1.0 / sum;
        for (int k = 0; k <= degree; ++k) {
            out[k] *= invSum;
        }
    }

    private static int determinePoissonUpperTruncationUnscaled(double mean, double epsilon) {
        final double target = 1.0 - epsilon;
        double term = Math.exp(-mean);
        double sum = term;

        int i = 0;
        final int hardCap = poissonLegacyHardCap(mean);

        while (sum < target && i < hardCap) {
            i++;
            term *= mean / i;
            sum += term;
        }
        return i;
    }

    private static void fillPoissonWeightsUnscaled(double mean, int maxN, double[] out, int offset) {
        Arrays.fill(out, offset, offset + maxN + 1, 0.0);
        double term = Math.exp(-mean);
        out[offset] = term;
        for (int n = 0; n < maxN; ++n) {
            term *= mean / (n + 1.0);
            out[offset + n + 1] = term;
        }
    }

    private static void fillBernsteinWeightsUnscaled(int degree, double x, double[] out) {
        Arrays.fill(out, 0, degree + 1, 0.0);

        if (degree == 0) {
            out[0] = 1.0;
            return;
        }
        if (!(x >= 0.0 && x <= 1.0) || Double.isNaN(x)) {
            throw new IllegalArgumentException("x must be in [0,1]");
        }
        if (x == 0.0) {
            out[0] = 1.0;
            return;
        }
        if (x == 1.0) {
            out[degree] = 1.0;
            return;
        }

        final double oneMinus = 1.0 - x;
        final double ratio = x / oneMinus;
        double weight = Math.pow(oneMinus, degree);
        out[0] = weight;

        for (int k = 0; k < degree; ++k) {
            weight *= ((degree - k) / (k + 1.0)) * ratio;
            out[k + 1] = weight;
        }
    }

    private static double logPoissonProbability(double mean, int n) {
        return -mean + n * Math.log(mean) - GammaFunction.lnGamma(n + 1.0);
    }

    private static int poissonMode(double mean) {
        if (mean < 1.0) {
            return 0;
        }
        if (mean > Integer.MAX_VALUE - 1.0) {
            throw new IllegalStateException(
                    "Poisson mean is too large for an indexed Sericola series: " + mean);
        }
        return (int) Math.floor(mean);
    }

    private static int poissonHardCap(double mean) {
        final double cap = Math.max(5000.0, mean + 12.0 * Math.sqrt(mean + 1.0) + 100.0);
        if (cap > Integer.MAX_VALUE - 1.0) {
            return Integer.MAX_VALUE - 1;
        }
        return (int) Math.ceil(cap);
    }

    private static int poissonLegacyHardCap(double mean) {
        final double cap = Math.max(5000.0, mean + 10.0 * Math.sqrt(mean + 1.0));
        if (cap > Integer.MAX_VALUE - 1.0) {
            return Integer.MAX_VALUE - 1;
        }
        return (int) cap;
    }

    private static void validateMean(double mean) {
        if (!(mean >= 0.0) || Double.isNaN(mean) || Double.isInfinite(mean)) {
            throw new IllegalArgumentException("Poisson mean must be finite and nonnegative: " + mean);
        }
    }
}
