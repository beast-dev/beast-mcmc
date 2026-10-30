/*
 * AbstractBranchHistoryModel.java
 *
 * Copyright © 2026, the BEAST Development Team.
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.treedatalikelihood.MAT;

import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.math.MathUtils;

import java.util.Arrays;

/**
 * Shared machinery for branch history models: the rate matrix, transition probabilities, the
 * integral over the times of an ordered sequence of hits and an exact sampler of hit paths given
 * their end states (uniformization; Hobolth and Stone 2009).
 */
public abstract class AbstractBranchHistoryModel implements BranchHistoryModel {

    protected final int stateCount;
    protected SubstitutionModel substitutionModel;
    protected double[] infinitesimalMatrix;
    protected double[] exitRates;

    // for sampling paths by uniformization: the rate mu (the largest exit rate) and the transition matrix
    // U = I + Q / mu of the uniformized chain, for the current rates; and work arrays
    private double uniformizationRate;
    private final double[] uniformizedMatrix;
    private double[] jumpRow;
    private double[] nextJumpRow;
    private final double[] stateWeights;

    // cache of transition probability matrices (and their logs) by distance, for the current rates: an
    // open-addressing table, emptied when the rates change or it is half full
    // (large enough for a distance per branch for each of several site rate categories)
    private static final int CACHE_CAPACITY = 32768;
    private static final int CACHE_HASH_SHIFT = 64 - Integer.numberOfTrailingZeros(CACHE_CAPACITY);
    private final long[] cachedDistances = new long[CACHE_CAPACITY];
    private final double[][] cachedTransitionProbabilities = new double[CACHE_CAPACITY][];
    private final double[][] cachedLogTransitionProbabilities = new double[CACHE_CAPACITY][];
    private int cachedCount = 0;

    protected AbstractBranchHistoryModel(int stateCount) {
        this.stateCount = stateCount;
        this.uniformizedMatrix = new double[stateCount * stateCount];
        this.jumpRow = new double[stateCount];
        this.nextJumpRow = new double[stateCount];
        this.stateWeights = new double[stateCount];
    }

    @Override
    public void setRates(SubstitutionModel substitutionModel, double[] infinitesimalMatrix, double[] exitRates) {
        this.substitutionModel = substitutionModel;
        this.infinitesimalMatrix = infinitesimalMatrix.clone();
        this.exitRates = exitRates.clone();
        clearTransitionProbabilityCache();

        double mu = 0.0;
        for (int s = 0; s < stateCount; s++) {
            mu = Math.max(mu, exitRates[s]);
        }
        uniformizationRate = mu;
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                uniformizedMatrix[i * stateCount + j] =
                        i == j ? 1.0 - exitRates[i] / mu : infinitesimalMatrix[i * stateCount + j] / mu;
            }
        }
        setEigenDecomposition(substitutionModel);
    }

    @Override
    public int getStateCount() {
        return stateCount;
    }

    // ---------------------------------------------------------------------------------------------
    // Transition probabilities from the eigen decomposition (without the cache)
    // ---------------------------------------------------------------------------------------------

    // the rate matrix's eigen decomposition, when the substitution model computes its transition
    // probabilities from it as BaseSubstitutionModel does (otherwise null)
    private double[] eigenValues;
    private double[] eigenVectors;
    private double[] inverseEigenVectors;
    private double[] eigenExp;

    private void setEigenDecomposition(SubstitutionModel substitutionModel) {
        eigenValues = null;
        if (!(substitutionModel instanceof BaseSubstitutionModel)) {
            return;
        }
        try {
            if (substitutionModel.getClass().getMethod("getTransitionProbabilities", double.class, double[].class)
                    .getDeclaringClass() != BaseSubstitutionModel.class) {
                return;
            }
        } catch (NoSuchMethodException e) {
            return;
        }
        EigenDecomposition eigen = substitutionModel.getEigenDecomposition();
        if (eigen == null || eigen.getEigenValues().length != stateCount ||
                eigen.getEigenVectors().length != stateCount * stateCount) {
            return;
        }
        eigenValues = eigen.getEigenValues().clone();
        eigenVectors = eigen.getEigenVectors().clone();
        inverseEigenVectors = eigen.getInverseEigenVectors().clone();
        eigenExp = new double[stateCount];
    }

    /**
     * @return true if transition probabilities can be computed directly ({@link #getLogTransitionDiagonal},
     * {@link #getLogTransitionProbability})
     */
    protected boolean hasDirectTransitionProbabilities() {
        return eigenValues != null;
    }

    /**
     * Sets out[offset + a] to log P_aa(distance) for every state a: the same values as the diagonal of
     * {@link #getLogTransitionProbabilities} (computed with the same arithmetic as BaseSubstitutionModel),
     * for S exponentials and logarithms rather than a whole matrix.
     */
    protected void getLogTransitionDiagonal(double distance, double[] out, int offset) {
        for (int i = 0; i < stateCount; i++) {
            eigenExp[i] = Math.exp(distance * eigenValues[i]);
        }
        for (int a = 0; a < stateCount; a++) {
            double p = 0.0;
            for (int k = 0; k < stateCount; k++) {
                p += eigenVectors[a * stateCount + k] * (inverseEigenVectors[k * stateCount + a] * eigenExp[k]);
            }
            out[offset + a] = Math.log(Math.max(Math.abs(p), 0.0));
        }
    }

    /**
     * @return log P_ab(distance), as {@link #getLogTransitionProbabilities} gives it
     */
    protected double getLogTransitionProbability(double distance, int a, int b) {
        double p = 0.0;
        for (int k = 0; k < stateCount; k++) {
            p += eigenVectors[a * stateCount + k] * (inverseEigenVectors[k * stateCount + b] * Math.exp(distance * eigenValues[k]));
        }
        return Math.log(Math.max(Math.abs(p), 0.0));
    }

    private void clearTransitionProbabilityCache() {
        if (cachedCount > 0) {
            MATProfile.count(MATProfile.Counter.TRANSITION_CACHE_CLEARS, 1);
            Arrays.fill(cachedTransitionProbabilities, null);
            Arrays.fill(cachedLogTransitionProbabilities, null);
            cachedCount = 0;
        }
    }

    @Override
    public void checkCompatibility(MATDataLikelihoodDelegate delegate) {
    }

    @Override
    public boolean needsComposition() {
        return false;
    }

    @Override
    public void initialise(DefaultMutationStore store) {
    }

    /**
     * @return the transition probabilities over the given distance (expected substitutions), row-major
     */
    protected double[] getTransitionProbabilities(double distance) {
        return cachedTransitionProbabilities[getCacheSlot(distance)];
    }

    protected double[] getLogTransitionProbabilities(double distance) {
        return cachedLogTransitionProbabilities[getCacheSlot(distance)];
    }

    /**
     * @return the slot of the cache holding the matrices for the distance (computed if not already there)
     */
    private int getCacheSlot(double distance) {
        final long key = Double.doubleToLongBits(distance);
        final int mask = CACHE_CAPACITY - 1;
        long h = key * 0x9E3779B97F4A7C15L;
        int slot = (int) (h >>> CACHE_HASH_SHIFT) & mask;
        while (cachedTransitionProbabilities[slot] != null) {
            if (cachedDistances[slot] == key) {
                return slot;
            }
            slot = (slot + 1) & mask;
        }
        if (cachedCount >= CACHE_CAPACITY / 2) {
            clearTransitionProbabilityCache();
            slot = (int) (h >>> CACHE_HASH_SHIFT) & mask;
        }
        MATProfile.count(MATProfile.Counter.TRANSITION_CACHE_MISSES, 1);
        double[] p = new double[stateCount * stateCount];
        substitutionModel.getTransitionProbabilities(distance, p);
        double[] logP = new double[p.length];
        for (int i = 0; i < p.length; i++) {
            logP[i] = Math.log(Math.max(p[i], 0.0));
        }
        cachedDistances[slot] = key;
        cachedTransitionProbabilities[slot] = p;
        cachedLogTransitionProbabilities[slot] = logP;
        cachedCount++;
        return slot;
    }

    /**
     * @return log(expm1(x) / x), i.e. log of (1/R) int_0^R exp(x s / R) ds, stable for all x
     */
    protected static double logExpm1OverX(double x) {
        if (Math.abs(x) < 1E-8) {
            return 0.5 * x;
        }
        if (x > 30.0) {
            return x - Math.log(x) + Math.log1p(-Math.exp(-x));
        }
        return Math.log(Math.expm1(x) / x);
    }

    // --------------------------------------------------------------------------------------------
    // Hit paths
    // --------------------------------------------------------------------------------------------

    /**
     * @return log of the integral over 0 < t_1 < ... < t_k < length of
     * exp(-sum_i c_i (t_{i+1} - t_i)), with t_0 = 0 and t_{k+1} = length: the probability that a path
     * through k+1 states with exit rates c_0..c_k makes its k jumps (and no others) in the branch,
     * excluding the jump rates themselves
     */
    static double logOrderedTimesIntegral(double[] c, int k, double length) {
        if (k == 0) {
            return -c[0] * length;
        }
        if (k == 1) {
            return -c[1] * length + Math.log(length) + logExpm1OverX((c[1] - c[0]) * length);
        }
        // Van Loan: the (0, k) element of exp(A) for the bidiagonal A with A_ii = -c_i L and A_i,i+1 = L,
        // shifted by the smallest rate so the exponential stays well scaled
        double m = Double.POSITIVE_INFINITY;
        for (int i = 0; i <= k; i++) {
            m = Math.min(m, c[i]);
        }
        int n = k + 1;
        double[][] a = new double[n][n];
        for (int i = 0; i < n; i++) {
            a[i][i] = -(c[i] - m) * length;
            if (i < k) {
                a[i][i + 1] = length;
            }
        }
        double[][] e = expm(a);
        return -m * length + Math.log(e[0][k]);
    }

    /**
     * Matrix exponential of a small matrix by scaling and squaring with a Taylor series.
     */
    static double[][] expm(double[][] a) {
        int n = a.length;
        double norm = 0.0;
        for (double[] row : a) {
            double sum = 0.0;
            for (double v : row) {
                sum += Math.abs(v);
            }
            norm = Math.max(norm, sum);
        }
        int squarings = norm > 0.5 ? (int) Math.ceil(Math.log(norm / 0.5) / Math.log(2.0)) : 0;
        double scale = Math.pow(2.0, -squarings);
        double[][] b = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                b[i][j] = a[i][j] * scale;
            }
        }
        double[][] result = identity(n);
        double[][] term = identity(n);
        for (int order = 1; order <= 20; order++) {
            term = multiply(term, b);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    term[i][j] /= order;
                    result[i][j] += term[i][j];
                }
            }
        }
        for (int s = 0; s < squarings; s++) {
            result = multiply(result, result);
        }
        return result;
    }

    private static double[][] identity(int n) {
        double[][] m = new double[n][n];
        for (int i = 0; i < n; i++) {
            m[i][i] = 1.0;
        }
        return m;
    }

    private static double[][] multiply(double[][] x, double[][] y) {
        int n = x.length;
        double[][] z = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < n; k++) {
                if (x[i][k] != 0.0) {
                    for (int j = 0; j < n; j++) {
                        z[i][j] += x[i][k] * y[k][j];
                    }
                }
            }
        }
        return z;
    }

    /**
     * Samples a path of hits from state a to state b over a branch of length L at rate nu, exactly from
     * the conditional distribution of a continuous-time Markov chain given its end states, by
     * uniformization.
     *
     * @param states filled with the states visited, states[0] = a ... states[k] = b
     * @param times  filled with the positions of the hits as fractions of the branch, times[1..k]
     * @return k, the number of hits
     */
    protected int samplePath(int a, int b, double length, double nu, int[][] states, double[][] times) {
        final double distance = nu * length;
        final double[] u = uniformizedMatrix;
        final int n = drawJumpCount(a, b, distance);

        // columns (U^m)_.b for m = 0..n, to draw the intermediate states
        double[][] toEnd = new double[n + 1][stateCount];
        toEnd[0][b] = 1.0;
        for (int m = 1; m <= n; m++) {
            for (int i = 0; i < stateCount; i++) {
                double sum = 0.0;
                for (int j = 0; j < stateCount; j++) {
                    sum += u[i * stateCount + j] * toEnd[m - 1][j];
                }
                toEnd[m][i] = sum;
            }
        }
        double[] jumpTimes = new double[n];
        for (int i = 0; i < n; i++) {
            jumpTimes[i] = MathUtils.nextDouble();
        }
        Arrays.sort(jumpTimes);

        // the state after each jump, given the state before it and that the path ends in b: real jumps
        // (to a different state) go in the path, virtual ones do not
        if (states[0] == null || states[0].length < n + 1) {
            states[0] = new int[n + 1];
        }
        if (times[0] == null || times[0].length < n + 1) {
            times[0] = new double[n + 1];
        }
        final double[] weights = stateWeights;
        int k = 0;
        int current = a;
        states[0][0] = a;
        for (int i = 1; i <= n; i++) {
            double total = 0.0;
            for (int x = 0; x < stateCount; x++) {
                weights[x] = u[current * stateCount + x] * toEnd[n - i][x];
                total += weights[x];
            }
            double draw = MathUtils.nextDouble() * total;
            int next = stateCount - 1;
            for (int x = 0; x < stateCount; x++) {
                if (draw < weights[x]) {
                    next = x;
                    break;
                }
                draw -= weights[x];
            }
            if (next != current) {
                k++;
                states[0][k] = next;
                times[0][k] = jumpTimes[i - 1];
                current = next;
            }
        }
        if (current != b) {
            throw new IllegalStateException("Path sampler failed to reach the end state");
        }
        return k;
    }

    /**
     * @return the number of jumps (real and virtual) of the uniformized chain over the distance, given that
     * it goes from a to b: P(n) is proportional to Poisson(n; mu d) (U^n)_ab
     */
    private int drawJumpCount(int a, int b, double distance) {
        final double mu = uniformizationRate;
        final double[] u = uniformizedMatrix;
        final double pab = getTransitionProbabilities(distance)[a * stateCount + b];
        final double target = MathUtils.nextDouble() * pab;
        double poisson = Math.exp(-mu * distance);
        double[] row = jumpRow;      // (U^n)_a.
        double[] next = nextJumpRow;
        Arrays.fill(row, 0.0);
        row[a] = 1.0;
        double cumulative = 0.0;
        int n = 0;
        while (true) {
            cumulative += poisson * row[b];
            if (cumulative >= target || n > 10000) {
                break;
            }
            Arrays.fill(next, 0.0);
            for (int i = 0; i < stateCount; i++) {
                if (row[i] != 0.0) {
                    for (int j = 0; j < stateCount; j++) {
                        next[j] += row[i] * u[i * stateCount + j];
                    }
                }
            }
            double[] swap = row;
            row = next;
            next = swap;
            n++;
            poisson *= mu * distance / n;
        }
        return n;
    }
}
