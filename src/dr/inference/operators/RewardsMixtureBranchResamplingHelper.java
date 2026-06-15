package dr.inference.operators;

import dr.math.MathUtils;
/*
 * @author Filippo Monti
 */
public final class RewardsMixtureBranchResamplingHelper {

    private static final double SUPPORT_BOUNDARY_TOLERANCE = 1.0e-12;

    private RewardsMixtureBranchResamplingHelper() {
        // no instances
    }

    public static final class BranchWeights {
        public final double[] logAtomicWeights;
        public final double logAtomicTotalWeight;
        public final double logCtsWeight;

        public BranchWeights(final double[] logAtomicWeights,
                             final double logAtomicTotalWeight,
                             final double logCtsWeight) {
            this.logAtomicWeights = logAtomicWeights;
            this.logAtomicTotalWeight = logAtomicTotalWeight;
            this.logCtsWeight = logCtsWeight;
        }
    }

    public static int sampleIndicatorFromLogs(final double logAtomicWeight,
                                              final double logCtsWeight) {

        final boolean atomicFinite = Double.isFinite(logAtomicWeight);
        final boolean ctsFinite = Double.isFinite(logCtsWeight);

        if (!atomicFinite && !ctsFinite) {
            throw new IllegalArgumentException(
                    "Invalid total weight: atomic=" + logAtomicWeight + ", cts=" + logCtsWeight
            );
        }
        if (atomicFinite && !ctsFinite) {
            return 1;
        }
        if (!atomicFinite) {
            return 0;
        }

        // pAtomic = exp(logA) / (exp(logA) + exp(logC))
        //         = 1 / (1 + exp(logC - logA))
        final double logRatio = logCtsWeight - logAtomicWeight;
        final double pAtomic;

        if (logRatio >= 0.0) {
            final double e = Math.exp(-logRatio);
            pAtomic = e / (1.0 + e);
        } else {
            final double e = Math.exp(logRatio);
            pAtomic = 1.0 / (1.0 + e);
        }

        return MathUtils.nextDouble() < pAtomic ? 1 : 0;
    }

    public static int sampleAtomFromLogs(final double[] logAtomicWeights,
                                         final int nstates) {

        double logTotal = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < nstates; i++) {
            final double lw = logAtomicWeights[i];
            if (Double.isNaN(lw)) {
                throw new IllegalArgumentException("Invalid atomic log-weight at state " + i + ": " + lw);
            }
            logTotal = logAdd(logTotal, lw);
        }

        if (!Double.isFinite(logTotal)) {
            throw new IllegalArgumentException("All atomic log-weights are -Infinity.");
        }

        final double u = MathUtils.nextDouble();
        double cumulative = 0.0;

        for (int i = 0; i < nstates; i++) {
            final double lw = logAtomicWeights[i];
            if (Double.isFinite(lw)) {
                cumulative += Math.exp(lw - logTotal);
            }
            if (u < cumulative) {
                return i;
            }
        }

        return nstates - 1;
    }

    public static double bilinearFormStable(final double[] pre,
                                            final double[] D,
                                            final double[] post,
                                            final int n) {
        return bilinearFormStable(pre, 0, D, post, 0, n);
    }

    public static double bilinearFormStable(final double[] pre,
                                            final int preOffset,
                                            final double[] D,
                                            final double[] post,
                                            final int postOffset,
                                            final int n) {
        double acc = 0.0;
        double cAcc = 0.0;

        for (int i = 0; i < n; i++) {
            final double prei = pre[preOffset + i];
            if (prei == 0.0) {
                continue;
            }

            final int rowBase = i * n;
            double rowDot = 0.0;
            double cRow = 0.0;

            for (int j = 0; j < n; j++) {
                final double y = D[rowBase + j] * post[postOffset + j] - cRow;
                final double t = rowDot + y;
                cRow = (t - rowDot) - y;
                rowDot = t;
            }

            final double y = prei * rowDot - cAcc;
            final double t = acc + y;
            cAcc = (t - acc) - y;
            acc = t;
        }

        return acc;
    }

    public static double logAtomicWeight(final double pre,
                                         final double post,
                                         final double logLocalFactor) {
        if (!(pre > 0.0) || !(post > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        return Math.log(pre) + Math.log(post) + logLocalFactor;
    }

    public static double logAtomicWeight(final double pre,
                                         final double post,
                                         final double logLocalFactor,
                                         final double preScale,
                                         final double postScale) {
        if (!(pre > 0.0) || !(post > 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        return Math.log(pre) + Math.log(post) + logLocalFactor + safeScale(preScale) + safeScale(postScale);
    }

    public static double logContinuousWeight(final double[] pre,
                                             final double[] D,
                                             final double[] post,
                                             final int nstates) {
        final double inner = bilinearFormStable(pre, D, post, nstates);
        if (!(inner > 0.0) || Double.isNaN(inner)) {
            return Double.NEGATIVE_INFINITY;
        }

        return Math.log(inner);
    }

    public static double logContinuousWeight(final double[] pre,
                                             final double[] D,
                                             final double[] post,
                                             final int nstates,
                                             final double preScale,
                                             final double postScale) {
        final double inner = bilinearFormStable(pre, D, post, nstates);
        if (!(inner > 0.0) || Double.isNaN(inner)) {
            return Double.NEGATIVE_INFINITY;
        }

        return Math.log(inner) + safeScale(preScale) + safeScale(postScale);
    }

    public static boolean isContinuousRewardOutsideOpenSupport(final double reward,
                                                               final double minReward,
                                                               final double maxReward) {
        if (!Double.isFinite(reward) || !Double.isFinite(minReward) || !Double.isFinite(maxReward)) {
            return true;
        }
        if (!(minReward < maxReward)) {
            return true;
        }

        final double scale = Math.max(1.0, Math.max(Math.abs(minReward), Math.abs(maxReward)));
        final double tolerance = SUPPORT_BOUNDARY_TOLERANCE * scale;
        return reward <= minReward + tolerance || reward >= maxReward - tolerance;
    }

    private static double safeScale(final double scale) {
        return Double.isFinite(scale) ? scale : 0.0;
    }

    public static double logSum(final double[] x, final int n) {
        double acc = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            acc = logAdd(acc, x[i]);
        }
        return acc;
    }

    public static double logAdd(final double a, final double b) {
        if (a == Double.NEGATIVE_INFINITY) return b;
        if (b == Double.NEGATIVE_INFINITY) return a;
        if (a < b) {
            return b + Math.log1p(Math.exp(a - b));
        } else {
            return a + Math.log1p(Math.exp(b - a));
        }
    }
}
