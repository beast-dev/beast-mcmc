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

    // cache of the most recent transition probability matrix
    private double cachedDistance = Double.NaN;
    private final double[] cachedTransitionProbabilities;
    private final double[] cachedLogTransitionProbabilities;

    protected AbstractBranchHistoryModel(int stateCount) {
        this.stateCount = stateCount;
        this.cachedTransitionProbabilities = new double[stateCount * stateCount];
        this.cachedLogTransitionProbabilities = new double[stateCount * stateCount];
    }

    @Override
    public void setRates(SubstitutionModel substitutionModel, double[] infinitesimalMatrix, double[] exitRates) {
        this.substitutionModel = substitutionModel;
        this.infinitesimalMatrix = infinitesimalMatrix.clone();
        this.exitRates = exitRates.clone();
        this.cachedDistance = Double.NaN;
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
        if (distance != cachedDistance) {
            substitutionModel.getTransitionProbabilities(distance, cachedTransitionProbabilities);
            for (int i = 0; i < cachedTransitionProbabilities.length; i++) {
                cachedLogTransitionProbabilities[i] = Math.log(Math.max(cachedTransitionProbabilities[i], 0.0));
            }
            cachedDistance = distance;
        }
        return cachedTransitionProbabilities;
    }

    protected double[] getLogTransitionProbabilities(double distance) {
        getTransitionProbabilities(distance);
        return cachedLogTransitionProbabilities;
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
            return -c[1] * length + Math.log(length) + MATDataLikelihoodDelegate.logExpm1OverX((c[1] - c[0]) * length);
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
        final double pab = getTransitionProbabilities(distance)[a * stateCount + b];
        double mu = 0.0;
        for (int s = 0; s < stateCount; s++) {
            mu = Math.max(mu, exitRates[s]);
        }
        // uniformized transition matrix U = I + Q / mu
        final double[] u = new double[stateCount * stateCount];
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                u[i * stateCount + j] = (i == j ? 1.0 - exitRates[i] / mu : infinitesimalMatrix[i * stateCount + j] / mu);
            }
        }

        // number of (real and virtual) jumps: P(n) proportional to Poisson(n; mu d) (U^n)_ab
        double target = MathUtils.nextDouble() * pab;
        double poisson = Math.exp(-mu * distance);
        double[] row = new double[stateCount];
        row[a] = 1.0;
        double cumulative = 0.0;
        int n = 0;
        while (true) {
            cumulative += poisson * row[b];
            if (cumulative >= target || n > 10000) {
                break;
            }
            double[] next = new double[stateCount];
            for (int i = 0; i < stateCount; i++) {
                if (row[i] != 0.0) {
                    for (int j = 0; j < stateCount; j++) {
                        next[j] += row[i] * u[i * stateCount + j];
                    }
                }
            }
            row = next;
            n++;
            poisson *= mu * distance / n;
        }

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

        if (states[0] == null || states[0].length < n + 1) {
            states[0] = new int[n + 1];
        }
        if (times[0] == null || times[0].length < n + 1) {
            times[0] = new double[n + 1];
        }
        int k = 0;
        int current = a;
        states[0][0] = a;
        for (int i = 1; i <= n; i++) {
            double total = 0.0;
            double[] weights = new double[stateCount];
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
}
