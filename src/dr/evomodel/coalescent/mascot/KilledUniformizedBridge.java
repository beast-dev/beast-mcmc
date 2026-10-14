/*
 * KilledUniformizedBridge.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.inference.markovjumps.StateChange;
import dr.inference.markovjumps.StateHistory;
import dr.math.MathUtils;

import java.util.Arrays;

/**
 * Endpoint-conditioned simulation of a homogeneous continuous-time Markov chain over one time cell by uniformization
 * (Hobolth and Stone 2009, Ann Appl Stat 3:1204-1231), for a generator G that may be substochastic: rows may sum to
 * less than zero, the deficit being a killing rate. With mu = max_s(-G[s][s]) and R = I + G / mu (entrywise
 * nonnegative, rows summing to at most one), the cell's transition matrix is T = exp(h G) = sum_n Pois(n; mu h) R^n.
 * Killing only lowers the row sums of R, so the series and the bridge construction are unchanged.
 * <p>
 * Every quantity here -- row and column actions, T[s][r], endpoint draws and paths -- uses the same truncated
 * series, so a filter computed with {@link #transitionRow} and draws from {@link #sampleStartAndPath} are exactly
 * consistent with each other.
 */
public final class KilledUniformizedBridge {

    /** Callers should split a cell so that mu h stays at or below this, which bounds the number of series terms. */
    public static final double MAX_RECOMMENDED_RATE_TIMES_LENGTH = 20.0;

    private static final double TAIL_TOLERANCE = 1.0e-16;
    private static final double MAX_RATE_TIMES_LENGTH = 600.0;
    private static final double ROW_SUM_TOLERANCE = 1.0e-12;

    private final int stateCount;
    private final double[] uniformized;
    private final double[] scratchA;
    private final double[] scratchB;
    private double[] poisson = new double[16];
    private double[] columnPowers = new double[0];
    private double[] jumpTimes = new double[8];
    private double rate;
    private double length;
    private int truncation;
    private boolean generatorSet;

    public KilledUniformizedBridge(int stateCount) {
        if (stateCount < 1) {
            throw new IllegalArgumentException("stateCount must be positive: " + stateCount);
        }
        this.stateCount = stateCount;
        this.uniformized = new double[stateCount * stateCount];
        this.scratchA = new double[stateCount];
        this.scratchB = new double[stateCount];
    }

    /**
     * Sets the cell's generator (row-major, stateCount x stateCount) and length. Off-diagonal entries must be
     * nonnegative and every row must sum to at most zero.
     */
    public void setGenerator(double[] generator, double cellLength) {
        int K = stateCount;
        if (generator.length != K * K) {
            throw new IllegalArgumentException("generator length " + generator.length + " is not " + (K * K));
        }
        if (!(cellLength >= 0.0) || !Double.isFinite(cellLength)) {
            throw new IllegalArgumentException("invalid cell length: " + cellLength);
        }
        double mu = 0.0;
        for (int s = 0; s < K; s++) {
            double rowSum = 0.0;
            for (int t = 0; t < K; t++) {
                double g = generator[s * K + t];
                if (!Double.isFinite(g) || (t != s && g < 0.0)) {
                    throw new IllegalArgumentException("invalid generator entry [" + s + "][" + t + "]: " + g);
                }
                rowSum += g;
            }
            mu = Math.max(mu, -generator[s * K + s]);
            if (rowSum > ROW_SUM_TOLERANCE * Math.max(1.0, -generator[s * K + s])) {
                throw new IllegalArgumentException("generator row " + s + " sums to " + rowSum + " > 0");
            }
        }

        rate = mu;
        length = cellLength;
        generatorSet = true;
        if (mu == 0.0 || cellLength == 0.0) {
            Arrays.fill(uniformized, 0.0);
            for (int s = 0; s < K; s++) {
                uniformized[s * K + s] = 1.0;
            }
            poisson[0] = 1.0;
            truncation = 0;
            return;
        }

        double muh = mu * cellLength;
        if (muh > MAX_RATE_TIMES_LENGTH) {
            throw new IllegalArgumentException("mu h = " + muh + " is too large; split the cell");
        }
        for (int s = 0; s < K; s++) {
            for (int t = 0; t < K; t++) {
                uniformized[s * K + t] = generator[s * K + t] / mu;
            }
            uniformized[s * K + s] = Math.max(0.0, 1.0 + generator[s * K + s] / mu);
        }

        // Keep at least K - 1 terms so that every state reachable within the cell keeps a positive probability.
        int minimumTerms = Math.max(K - 1, (int) Math.ceil(muh));
        double p = Math.exp(-muh);
        poisson[0] = p;
        double accumulated = p;
        int n = 0;
        while (n < minimumTerms || p >= TAIL_TOLERANCE * accumulated) {
            n++;
            p *= muh / n;
            if (n >= poisson.length) {
                poisson = Arrays.copyOf(poisson, 2 * poisson.length);
            }
            poisson[n] = p;
            accumulated += p;
        }
        truncation = n;
    }

    public double getUniformizationRate() {
        requireGenerator();
        return rate;
    }

    public int getTruncation() {
        requireGenerator();
        return truncation;
    }

    /** out = in T. {@code in} and {@code out} may be the same array. */
    public void transitionRow(double[] in, int inOffset, double[] out, int outOffset) {
        requireGenerator();
        int K = stateCount;
        double[] u = scratchA;
        double[] next = scratchB;
        System.arraycopy(in, inOffset, u, 0, K);
        double p0 = poisson[0];
        for (int s = 0; s < K; s++) {
            out[outOffset + s] = p0 * u[s];
        }
        for (int n = 1; n <= truncation; n++) {
            for (int t = 0; t < K; t++) {
                double sum = 0.0;
                for (int s = 0; s < K; s++) {
                    sum += u[s] * uniformized[s * K + t];
                }
                next[t] = sum;
            }
            double[] swap = u;
            u = next;
            next = swap;
            double pn = poisson[n];
            for (int s = 0; s < K; s++) {
                out[outOffset + s] += pn * u[s];
            }
        }
    }

    /** out = T in. {@code in} and {@code out} may be the same array. */
    public void transitionColumn(double[] in, int inOffset, double[] out, int outOffset) {
        requireGenerator();
        int K = stateCount;
        double[] u = scratchA;
        double[] next = scratchB;
        System.arraycopy(in, inOffset, u, 0, K);
        double p0 = poisson[0];
        for (int s = 0; s < K; s++) {
            out[outOffset + s] = p0 * u[s];
        }
        for (int n = 1; n <= truncation; n++) {
            multiplyColumn(u, next);
            double[] swap = u;
            u = next;
            next = swap;
            double pn = poisson[n];
            for (int s = 0; s < K; s++) {
                out[outOffset + s] += pn * u[s];
            }
        }
    }

    /** T[from][to] from the same truncated series as every other method. */
    public double transitionProbability(int from, int to) {
        requireGenerator();
        computeColumnPowers(to);
        double sum = 0.0;
        for (int n = 0; n <= truncation; n++) {
            sum += poisson[n] * columnPowers[n * stateCount + from];
        }
        return sum;
    }

    /**
     * Draws the state at the start of the cell jointly with the path through the cell, given the state at its end.
     * The start state s is drawn with probability proportional to {@code startWeights[s] T[s][endState]} and the
     * path from the chain conditioned on both endpoints. Jumps are placed strictly inside
     * {@code (startTime, endTime)} and appended to {@code path} in increasing time; virtual (self) jumps are not
     * recorded. Returns the start state.
     *
     * @throws IllegalStateException if no start state has positive weight
     */
    public int sampleStartAndPath(double[] startWeights, int startWeightsOffset, int endState,
                                  double startTime, double endTime, Path path) {
        requireGenerator();
        int K = stateCount;
        computeColumnPowers(endState);

        double total = 0.0;
        for (int n = 0; n <= truncation; n++) {
            double pn = poisson[n];
            int offset = n * K;
            for (int s = 0; s < K; s++) {
                total += startWeights[startWeightsOffset + s] * pn * columnPowers[offset + s];
            }
        }
        if (!(total > 0.0) || !Double.isFinite(total)) {
            throw new IllegalStateException("no start state can reach end state " + endState +
                    " in the cell [" + startTime + ", " + endTime + "] (total weight " + total + ", start weights " +
                    Arrays.toString(Arrays.copyOfRange(startWeights, startWeightsOffset, startWeightsOffset + K)) +
                    ")");
        }

        double u = MathUtils.nextDouble() * total;
        int chosenState = -1;
        int chosenJumps = -1;
        double cumulative = 0.0;
        search:
        for (int n = 0; n <= truncation; n++) {
            double pn = poisson[n];
            int offset = n * K;
            for (int s = 0; s < K; s++) {
                double w = startWeights[startWeightsOffset + s] * pn * columnPowers[offset + s];
                if (w > 0.0) {
                    chosenState = s;
                    chosenJumps = n;
                    cumulative += w;
                    if (u < cumulative) {
                        break search;
                    }
                }
            }
        }

        int jumps = chosenJumps;
        if (jumps == 0) {
            return chosenState;
        }

        if (jumpTimes.length < jumps) {
            jumpTimes = new double[Math.max(jumps, 2 * jumpTimes.length)];
        }
        double span = endTime - startTime;
        for (int j = 0; j < jumps; j++) {
            double fraction;
            do {
                fraction = MathUtils.nextDouble();
            } while (fraction == 0.0);
            jumpTimes[j] = startTime + fraction * span;
        }
        Arrays.sort(jumpTimes, 0, jumps);

        int current = chosenState;
        for (int j = 1; j <= jumps; j++) {
            int next;
            if (j == jumps) {
                next = endState;
            } else {
                next = drawIntermediateState(current, jumps - j);
            }
            if (next != current) {
                path.add(jumpTimes[j - 1], next);
            }
            current = next;
        }
        return chosenState;
    }

    /** A complete endpoint-conditioned history from {@code startState} to {@code endState}, for testing. */
    public StateHistory sampleBridge(int startState, int endState, double startTime, double endTime) {
        double[] weights = new double[stateCount];
        weights[startState] = 1.0;
        Path path = new Path();
        sampleStartAndPath(weights, 0, endState, startTime, endTime, path);
        StateHistory history = new StateHistory(startTime, startState, stateCount);
        for (int i = 0; i < path.size(); i++) {
            history.addChange(new StateChange(path.getTime(i), path.getState(i)));
        }
        history.addEndingState(new StateChange(endTime, endState));
        return history;
    }

    /** Draws z with probability proportional to R[current][z] (R^remaining e_end)[z]. */
    private int drawIntermediateState(int current, int remaining) {
        int K = stateCount;
        int offset = remaining * K;
        double total = 0.0;
        for (int z = 0; z < K; z++) {
            total += uniformized[current * K + z] * columnPowers[offset + z];
        }
        double u = MathUtils.nextDouble() * total;
        double cumulative = 0.0;
        int chosen = -1;
        for (int z = 0; z < K; z++) {
            double w = uniformized[current * K + z] * columnPowers[offset + z];
            if (w > 0.0) {
                chosen = z;
                cumulative += w;
                if (u < cumulative) {
                    break;
                }
            }
        }
        if (chosen < 0) {
            throw new IllegalStateException("bridge reached a state with no path to its end state");
        }
        return chosen;
    }

    /** columnPowers[n K + s] = (R^n e_end)[s] for n = 0..truncation. */
    private void computeColumnPowers(int endState) {
        int K = stateCount;
        if (endState < 0 || endState >= K) {
            throw new IllegalArgumentException("end state out of range: " + endState);
        }
        int needed = (truncation + 1) * K;
        if (columnPowers.length < needed) {
            columnPowers = new double[Math.max(needed, 2 * columnPowers.length)];
        }
        Arrays.fill(columnPowers, 0, K, 0.0);
        columnPowers[endState] = 1.0;
        for (int n = 1; n <= truncation; n++) {
            int previous = (n - 1) * K;
            int offset = n * K;
            for (int s = 0; s < K; s++) {
                double sum = 0.0;
                for (int t = 0; t < K; t++) {
                    sum += uniformized[s * K + t] * columnPowers[previous + t];
                }
                columnPowers[offset + s] = sum;
            }
        }
    }

    private void multiplyColumn(double[] in, double[] out) {
        int K = stateCount;
        for (int s = 0; s < K; s++) {
            double sum = 0.0;
            for (int t = 0; t < K; t++) {
                sum += uniformized[s * K + t] * in[t];
            }
            out[s] = sum;
        }
    }

    private void requireGenerator() {
        if (!generatorSet) {
            throw new IllegalStateException("setGenerator must be called first");
        }
    }

    /** Growable list of (time, state) jumps. */
    public static final class Path {
        private double[] times = new double[16];
        private int[] states = new int[16];
        private int size;

        public void add(double time, int state) {
            if (size == times.length) {
                times = Arrays.copyOf(times, 2 * size);
                states = Arrays.copyOf(states, 2 * size);
            }
            times[size] = time;
            states[size] = state;
            size++;
        }

        public void clear() {
            size = 0;
        }

        public int size() {
            return size;
        }

        public double getTime(int index) {
            return times[index];
        }

        public int getState(int index) {
            return states[index];
        }
    }
}
