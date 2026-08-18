package dr.inference.operators;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

import java.util.Arrays;

/**
 * Full-conditional Gibbs refresh for the state-to-reward-rate-slot mapping
 * ("rewardRatesMapping" / RewardRates.getStateIndices()).
 *
 * {@link OneZeroOneShuffleOperator} proposes a single extreme<->internal
 * transposition and accepts it via Metropolis-Hastings. That move set can
 * leave some target permutations unreachable in one step (any permutation
 * that differs from the current mapping by an extreme<->extreme or an
 * internal<->internal transposition, or a composition of such moves with no
 * monotonically-improving intermediate state) with no accepted single-step
 * path to climb, so it can get stuck arbitrarily far from a much better
 * mapping.
 *
 * This operator instead enumerates every one of the K! valid permutations of
 * the mapping, evaluates the full joint log-likelihood under each (primary
 * treeDataLikelihood plus any dependent CTMC/continuous likelihoods that also
 * depend on the reward rates), and Gibbs-samples the new mapping proportional
 * to relative likelihood. Because the mapping change is a global
 * invalidation (it can change the continuous-reward transition matrix used
 * by every branch, not just branches touching the swapped states -- see
 * RewardsAwareBranchModel.handleModelChangedEvent), there is no branch-local
 * sufficient-statistics shortcut available here the way there is for
 * RewardMixtureCategoricalGibbsOperator; each candidate requires a full
 * "set, invalidate, re-evaluate" pass, mirroring the cluster-move path in
 * RewardsMixtureIndicatorAndAtomIndicesOperator.
 *
 * Exhaustive enumeration is an exact draw from the true full conditional (the
 * candidate set is the entire state space for every current state, so the
 * normalizing constant is the same regardless of where the chain currently
 * is), which is why this can be accepted unconditionally as a GibbsOperator
 * with no Metropolis-Hastings correction. A cheaper alternative -- enumerate
 * only the K(K-1)/2 pairwise-swap neighbors of the current state -- was
 * considered and rejected: that neighborhood is current-state-dependent, so
 * its normalizing constant differs from state to state and a naive
 * proportional-to-likelihood draw among only those neighbors would not
 * satisfy detailed balance without an additional Multiple-Try-Metropolis-style
 * correction.
 *
 * Cost: K! full-model log-likelihood evaluations per call. This is cheap for
 * the K in {3,4,5} this project mostly uses (<=120 evaluations); K=7 (5040
 * evaluations) can be expensive if this operator is given a weight that makes
 * it fire often. No fallback for large K is implemented here.
 *
 * @author Filippo Monti
 */
public final class OneZeroOneShuffleGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private final Parameter rewardRatesValues;
    private final Parameter rewardRatesMapping;
    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeDataLikelihood[] dependentCtmcLikelihoods;
    private final TreeDataLikelihood[] dependentContinuousLikelihoods;
    private final double tol;

    private final int stateCount;
    private final int permutationCount;
    private final int[][] permutations;
    private final double[] logWeights;
    private final int[] storedMapping;

    public OneZeroOneShuffleGibbsOperator(final Parameter rewardRatesValues,
                                          final Parameter rewardRatesMapping,
                                          final TreeDataLikelihood treeDataLikelihood,
                                          final TreeDataLikelihood[] dependentCtmcLikelihoods,
                                          final TreeDataLikelihood[] dependentContinuousLikelihoods,
                                          final double weight,
                                          final double tol) {
        if (rewardRatesValues == null) {
            throw new IllegalArgumentException("rewardRatesValues must be non-null");
        }
        if (rewardRatesMapping == null) {
            throw new IllegalArgumentException("rewardRatesMapping must be non-null");
        }
        if (treeDataLikelihood == null) {
            throw new IllegalArgumentException("treeDataLikelihood must be non-null");
        }
        if (rewardRatesValues.getDimension() < 2) {
            throw new IllegalArgumentException("rewardRatesValues dimension must be >= 2");
        }
        if (tol < 0.0) {
            throw new IllegalArgumentException("tol must be >= 0");
        }

        this.rewardRatesValues = rewardRatesValues;
        this.rewardRatesMapping = rewardRatesMapping;
        this.treeDataLikelihood = treeDataLikelihood;
        this.dependentCtmcLikelihoods = dependentCtmcLikelihoods == null
                ? new TreeDataLikelihood[0]
                : Arrays.copyOf(dependentCtmcLikelihoods, dependentCtmcLikelihoods.length);
        this.dependentContinuousLikelihoods = dependentContinuousLikelihoods == null
                ? new TreeDataLikelihood[0]
                : Arrays.copyOf(dependentContinuousLikelihoods, dependentContinuousLikelihoods.length);
        this.tol = tol;

        this.stateCount = rewardRatesValues.getDimension();
        if (rewardRatesMapping.getDimension() != stateCount) {
            throw new IllegalArgumentException(
                    "rewardRatesMapping dimension must equal rewardRatesValues dimension. Found " +
                            rewardRatesMapping.getDimension() + " but expected " + stateCount);
        }
        for (final TreeDataLikelihood dependent : this.dependentCtmcLikelihoods) {
            if (dependent == null) {
                throw new IllegalArgumentException("dependentCtmcLikelihoods contains a null entry");
            }
        }
        for (final TreeDataLikelihood dependent : this.dependentContinuousLikelihoods) {
            if (dependent == null) {
                throw new IllegalArgumentException("dependentContinuousLikelihoods contains a null entry");
            }
        }

        validateRewardRatesValues();

        this.permutations = enumeratePermutations(stateCount);
        this.permutationCount = permutations.length;
        this.logWeights = new double[permutationCount];
        this.storedMapping = new int[stateCount];

        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return "oneZeroOneShuffleGibbsOperator(" + rewardRatesMapping.getParameterName() + ")";
    }

    @Override
    public double doOperation() {

        storeCurrentMapping();

        for (int p = 0; p < permutationCount; p++) {
            applyMapping(permutations[p]);
            logWeights[p] = computeLogTarget();
        }

        final double logTotal = RewardsMixtureBranchResamplingHelper.logSum(logWeights, permutationCount);

        if (!Double.isFinite(logTotal)) {
            applyMapping(storedMapping);
            return Double.NEGATIVE_INFINITY;
        }

        final int chosen = sampleIndex(logTotal);
        applyMapping(permutations[chosen]);

        return 0.0;
    }

    private int sampleIndex(final double logTotal) {
        final double u = MathUtils.nextDouble();
        double cumulative = 0.0;
        for (int p = 0; p < permutationCount; p++) {
            if (Double.isFinite(logWeights[p])) {
                cumulative += Math.exp(logWeights[p] - logTotal);
            }
            if (u < cumulative) {
                return p;
            }
        }
        return permutationCount - 1;
    }

    private double computeLogTarget() {
        treeDataLikelihood.makeDirty();
        double logTarget = treeDataLikelihood.getLogLikelihood();
        for (final TreeDataLikelihood dependent : dependentCtmcLikelihoods) {
            dependent.makeDirty();
            logTarget += dependent.getLogLikelihood();
        }
        for (final TreeDataLikelihood dependent : dependentContinuousLikelihoods) {
            dependent.makeDirty();
            logTarget += dependent.getLogLikelihood();
        }
        return logTarget;
    }

    private void storeCurrentMapping() {
        for (int s = 0; s < stateCount; s++) {
            storedMapping[s] = toIndexStrict(rewardRatesMapping.getParameterValue(s), "rewardRatesMapping[" + s + "]");
        }
    }

    private void applyMapping(final int[] mapping) {
        for (int s = 0; s < stateCount; s++) {
            rewardRatesMapping.setParameterValueQuietly(s, mapping[s]);
        }
        // A quiet set alone does not invalidate RewardsAwareBranchModel's cached
        // no-jump rates / atomic scales / continuous transition matrices (those
        // are only invalidated by handleModelChangedEvent, which only fires in
        // response to an actual change event) -- treeDataLikelihood.makeDirty()
        // forces a recompute, but with stale cached matrices unless this fires.
        rewardRatesMapping.fireParameterChangedEvent();
    }

    private void validateRewardRatesValues() {
        final double v0 = rewardRatesValues.getParameterValue(0);
        final double v1 = rewardRatesValues.getParameterValue(1);

        if (Math.abs(v0 - 0.0) > tol) {
            throw new IllegalArgumentException("rewardRatesValues[0] must be 0.0 but found " + v0);
        }
        if (Math.abs(v1 - 1.0) > tol) {
            throw new IllegalArgumentException("rewardRatesValues[1] must be 1.0 but found " + v1);
        }
        for (int j = 0; j < stateCount; j++) {
            final double a = rewardRatesValues.getParameterValue(j);
            if (Double.isNaN(a) || a < -tol || a > 1.0 + tol) {
                throw new IllegalArgumentException("rewardRatesValues out of [0,1] at index " + j + ": " + a);
            }
        }
    }

    private static int toIndexStrict(final double x, final String name) {
        final long r = Math.round(x);
        if (Math.abs(x - r) > 1e-9) {
            throw new IllegalArgumentException(name + " must be integer-like but found " + x);
        }
        if (r < Integer.MIN_VALUE || r > Integer.MAX_VALUE) {
            throw new IllegalArgumentException(name + " out of int range: " + x);
        }
        return (int) r;
    }

    /**
     * Enumerate all n! permutations of {0,...,n-1} via Heap's algorithm.
     * Returns an array of freshly-allocated int[] arrays (no shared/mutated
     * buffer), one per permutation, in the order Heap's algorithm produces
     * them (not lexicographic, but every permutation appears exactly once).
     */
    static int[][] enumeratePermutations(final int n) {
        if (n <= 0) {
            throw new IllegalArgumentException("n must be positive: " + n);
        }
        long count = 1;
        for (int i = 2; i <= n; i++) {
            count *= i;
        }
        if (count > Integer.MAX_VALUE / 2) {
            throw new IllegalArgumentException("n! is too large to enumerate: n=" + n);
        }

        final java.util.List<int[]> results = new java.util.ArrayList<int[]>((int) count);
        final int[] current = new int[n];
        for (int i = 0; i < n; i++) {
            current[i] = i;
        }
        final int[] c = new int[n];

        results.add(Arrays.copyOf(current, n));
        int i = 0;
        while (i < n) {
            if (c[i] < i) {
                if (i % 2 == 0) {
                    swap(current, 0, i);
                } else {
                    swap(current, c[i], i);
                }
                results.add(Arrays.copyOf(current, n));
                c[i]++;
                i = 0;
            } else {
                c[i] = 0;
                i++;
            }
        }

        return results.toArray(new int[0][]);
    }

    private static void swap(final int[] a, final int i, final int j) {
        final int tmp = a[i];
        a[i] = a[j];
        a[j] = tmp;
    }
}
