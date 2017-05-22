/*
 * SericolaSeriesMarkovReward.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.markovjumps;

import dr.evomodel.substmodel.DefaultEigenSystem;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.EigenSystem;
import dr.math.Binomial;
import dr.math.GammaFunction;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;

/**
 * @author Marc Suchard
 * @author Forrest Crawford
 */
public class SericolaSeriesMarkovReward implements MarkovReward {

    // Following Bladt, Meini, Neuts and Sericola (2002).
    // Assuming each state has a distinct reward, i.e. \phi + 1 = stateCount,
    // and states are sorted in increasing reward order

    private static final boolean DEBUG = false;

    public SericolaSeriesMarkovReward(double[] Q, double[] r, int dim) {
        this(Q, r, dim, 1E-10);
    }

    public SericolaSeriesMarkovReward(double[] Q, double[] r, int dim, double epsilon) {
        this.Q = Q;
        this.r = r;
        this.maxTime = 0;
        this.epsilon = epsilon;

        this.dim = dim;
        lambda = determineLambda();

        phi = dim - 1;

        if (DEBUG) {
            System.err.println("lambda = " + lambda);
        }

        P = initializeP(Q, lambda);
        eigenSystem = new DefaultEigenSystem(dim);
    }

    private double[][] initializeW(int times, int dim) {
        return new double[times][dim * dim];
    }

    private int getHfromX(double x, double time) {
        // TODO assert x > h[0] * time;
        int h = 1;
        while (x >= r[h] * time) {
            h++;
        }
        return h;
    }

    private void growC(double time, int extraN) {
        int newN = getNfromC();
        if (time > maxTime) {
            newN = determineNumberOfSteps(time, lambda) + extraN;
            maxTime = time;
        }

        // Grow C if necessary
        if (newN > getNfromC()) {
            if (DEBUG) {
                System.err.println("Growing C to N = " + newN + " with " + maxTime);
            }
            if (newN > 500) {
                System.err.println("Warning: > 500 recursion depth in SericolaSeriesMarkovReward");
            }
            initializeSpace(phi, newN);
            computeChnk();
        }
    }

    // START: internal structure of C

    private double[][][][] internalC; // TODO Linearize for store/restore; TODO reduce to minimal storage

    private void initializeSpace(int phi, int N) {
        internalC = new double[phi + 1][N + 1][N + 1][dim * dim];
        // indices [h][n][k][B_u][B_v]
    }

    private double[] C(int h, int n, int k) {
        return internalC[h][n][k];
    }

    private int getNfromC() {
        return (internalC == null) ? -1 : internalC[0].length - 1;
    }

    private int idx(int i, int j) {
        return i * dim + j; // row-major
    }

    // END: internal structure of C, TODO Change to expandable list

    private int[] getHfromX(double[] X, double time) {
        int[] H = new int[X.length];
        for (int i = 0; i < X.length; ++i) {
            H[i] = getHfromX(X[i], time);
        }
        return H;
//        return new int[] { 1 };      // AR nasty hack - revert shortly
    }

    public double computePdf(double x, double time, int i, int j) {
        if (x == time) return 0.0;
        else return computePdf(x, time)[i * dim + j];
    }

    public double[] computePdf(double x, double time) {
        return computePdf(new double[]{x}, time)[0];
    }

    public double[][] computePdf(double[] X, double time) {
        int[] H = getHfromX(X, time);

        growC(time, 1);

        double[][] W = initializeW(X.length, dim); // initialize with zeros

        final int N = getNfromC() - 1; // TODO N should be branch-length-specific to save computation
        for (int n = 0; n <= N; ++n) {
            accumulatePdf(W, X, H, n, time); // TODO This can be sped up when only a single entry is wanted
        }

        if (DEBUG) {
            for (int i = 0; i < W.length; ++i) {
                System.err.println("W'[" + i + "]:\n" + new Matrix(squareMatrix(W[i])));
            }
            System.err.println("");
        }

        return W;
    }

    public double computeCdf(double x, double time, int i, int j) {
        return computeCdf(x, time)[i * dim + j];
    }

    public double[] computeCdf(double x, double time) {
        return computeCdf(new double[]{x}, time)[0];
    }

    public double[][] computeCdf(double[] X, double time) {

        int[] H = getHfromX(X, time);

        growC(time, 0);

        double[][] W = initializeW(X.length, dim); // initialize with zeros

        final int N = getNfromC();
        for (int n = 0; n <= N; ++n) {
            accumulateCdf(W, X, H, n, time);
        }

        if (DEBUG) {
            for (int i = 0; i < W.length; ++i) {
                System.err.println("W[" + i + "]:\n" + new Matrix(squareMatrix(W[i])));
            }
            System.err.println("");
        }

        return W;
    }

    private double[] initializeP(double[] Q, double lambda) {
        double[] P = new double[dim * dim];
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                double identity = (i == j) ? 1.0 : 0.0;
                P[idx(i, j)] = identity + Q[idx(i, j)] / lambda;
            }
        }
        return P;
    }

    private void accumulateCdf(double[][] W, double[] X, int[] H, int n, double time) {

        final double premult = Math.exp(
                -lambda * time + n * (Math.log(lambda) + Math.log(time)) - GammaFunction.lnGamma(n + 1.0)
        );

        // TODO Make factorial/choose static look-up tables

        for (int t = 0; t < X.length; ++t) { // For each time point
            double x = X[t];
            int h = H[t];

            double xh = (x - r[h - 1] * time) / ((r[h] - r[h - 1]) * time);

            final int dim2 = dim * dim;
            double[] inc = new double[dim2]; // W^{\epsilon}(x(i),t,n)
            for (int k = 0; k <= n; k++) {
                final double binomialCoef = Binomial.choose(n, k) * Math.pow(xh, k) * Math.pow(1.0 - xh, n - k);
                for (int uv = 0; uv < dim2; ++uv) {
                    inc[uv] += binomialCoef * C(h, n, k)[uv];
                }
            }

            for (int uv = 0; uv < dim2; ++uv) {
                W[t][uv] += premult * inc[uv];
            }
        }
    }

    private void accumulatePdf(double[][] W, double[] X, int[] H, int n, double time) {

        final double premult = Math.exp(
                -lambda * time + n * (Math.log(lambda) + Math.log(time)) - GammaFunction.lnGamma(n + 1.0)
        );

        // TODO Make factorial/choose static look-up tables
        // AR - Binomial has a look-up-table built in for k=2.


        for (int t = 0; t < X.length; ++t) { // For each time point
            double x = X[t];
            int h = H[t];

            final double factor = lambda / (r[h] - r[h - 1]);

            double xh = (x - r[h - 1] * time) / ((r[h] - r[h - 1]) * time);

            final int dim2 = dim * dim;
            double[] inc = new double[dim2]; // W^{\epsilon}(x(i),t,n)
            for (int k = 0; k <= n; k++) {
                final double binomialCoef = Binomial.choose(n, k) * Math.pow(xh, k) * Math.pow(1.0 - xh, n - k);
                for (int uv = 0; uv < dim2; ++uv) {
                    inc[uv] += binomialCoef * (C(h, n + 1, k + 1)[uv] - C(h, n + 1, k)[uv]);
                }
            }

            for (int uv = 0; uv < dim2; ++uv) {
                W[t][uv] += factor * premult * inc[uv];
            }
        }
    }

    private double relationTwelve(int h, int n, int k, int u, int v) {
        // TODO ratios are independent of u,v,w
        double c = (r[u] - r[h]) / (r[u] - r[h - 1]) * C(h, n, k - 1)[idx(u, v)];

        double d = 0;
        for (int w = 0; w <= phi; ++w) {
            d += P[idx(u, w)] * C(h, n - 1, k - 1)[idx(w, v)];
        }
        d *= (r[h] - r[h - 1]) / (r[u] - r[h - 1]);

        return c + d;
    }

    private double relationThirteen(int h, int n, int k, int u, int v) {
        // TODO ratios of are independent of u,v,w
        double c = (r[h - 1] - r[u]) / (r[h] - r[u]) * C(h, n, k + 1)[idx(u, v)];

        double d = 0;
        for (int w = 0; w <= phi; ++w) {
            d += P[idx(u, w)] * C(h, n - 1, k)[idx(w, v)];
        }
        d *= (r[h] - r[h - 1]) / (r[h] - r[u]);

        return c + d;
    }

    private double[] product(double[] a, double[] b) {
        double[] c = new double[dim * dim];
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                final int ij = idx(i, j);
                for (int k = 0; k < dim; ++k) {
                    c[ij] += a[idx(i, k)] * b[idx(k, j)];
                }
            }
        }
        return c;
    }

    private void computeChnk() {

        double[] Pn = new double[dim * dim];
        for (int u = 0; u < dim; ++u) {
            Pn[idx(u, u)] = 1.0;
        }

        for (int h = 1; h <= phi; ++h) {
            // zero corner cases automatically
            for (int u = 0; u <= h - 1; ++u) {
                C(h, 0, 0)[idx(u, u)] = 1.0;
            }
        }

//        accumulate(0);
        final int N = getNfromC();
        for (int n = 1; n <= N; ++n) {
            // zero corner cases automatically
            for (int h = 1; h <= phi; ++h) {
                for (int k = 1; k <= n; ++k) {

                    for (int u = h; u <= phi; ++u) {
                        for (int v = 0; v <= phi; ++v) {
                            C(h, n, k)[idx(u, v)] = relationTwelve(h, n, k, u, v);
                        }
                    }

                }
                for (int u = h + 1; u <= phi; ++u) {
                    for (int v = 0; v <= phi; ++v) {
                        C(h + 1, n, 0)[idx(u, v)] = C(h, n, n)[idx(u, v)];
                    }
                }
            }

            Pn = product(Pn, P);
            for (int u = 0; u <= phi - 1; ++u) {
                for (int v = 0; v <= phi; ++v) {
                    C(phi, n, n)[idx(u, v)] = Pn[idx(u, v)];
                }
            }

            for (int h = phi; h >= 1; --h) {
                for (int k = n - 1; k >= 0; --k) {

                    for (int u = 0; u <= h - 1; u++) {
                        for (int v = 0; v <= phi; ++v) {
                            C(h, n, k)[idx(u, v)] = relationThirteen(h, n, k, u, v);
                        }
                    }

                    for (int u = 0; u <= h - 2; ++u) {
                        for (int v = 0; v <= phi; ++v) {
                            C(h - 1, n, n)[idx(u, v)] = C(h, n, 0)[idx(u, v)];
                        }
                    }
                }
            }

//            accumulate(n);
        }
    }

    private double determineLambda() {
        double lambda = Q[0]; // Q[idx(0,0)]
        for (int i = 1; i < dim; ++i) {
            int ii = idx(i, i);
            if (Q[ii] < lambda) {
                lambda = Q[ii];
            }
        }
        return -lambda;
    }

    private double[][] squareMatrix(final double[] mat) {
        double[][] rtn = new double[dim][dim];
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                rtn[i][j] = mat[idx(i, j)];
            }
        }
        return rtn;
    }

    private int determineNumberOfSteps(double time, double lambda) {
//        final double tolerance = (1.0 - epsilon) / Math.exp(-lambda * time);
//        final double logTolerance = Math.log(1.0 - epsilon);

        int i = -1;
//        double sum = 0.0;
//        int factorialI = 1;
//
//        while (sum < tolerance) {
//            i++;
//            sum += Math.pow(lambda * time, i) / factorialI;
//            factorialI *= (i + 1); // + 1 because used on next iterate
//        }
//
//        int firstN = i;
//
//                            LogTricks.logSum()
//        i = -1;
        final double tolerance2 =
//                -epsilon;
                (1.0 - epsilon);
        double sum2 = 0.0;
//        double sum2 = Double.NEGATIVE_INFINITY;
        while (Math.abs(sum2 - tolerance2) > epsilon && sum2 < 1.0) {
//        while (sum2 < tolerance2) {
            i++;
            double logDensity = -lambda * time + i * (Math.log(lambda) + Math.log(time)) - GammaFunction.lnGamma(i + 1);
            sum2 += Math.exp(logDensity);
//            sum2 = LogTricks.logSum(sum2, logDensity);
            if (DEBUG) {
                System.err.println(sum2 + " " + tolerance2 + " " + Math.abs(sum2 - tolerance2) + " " + epsilon * 0.01);
//            if (i > 500) System.exit(-1);
            }
        }

//        System.err.println("First: " + firstN);
//        System.err.println("Second:" + i);
//        System.exit(-1);

        return i;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Q: " + new Vector(Q) + "\n");
        sb.append("r: " + new Vector(r) + "\n");
        sb.append("lambda: " + lambda + "\n");
        sb.append("N: " + getNfromC() + "\n");
        sb.append("maxTime: " + maxTime + "\n");
        sb.append("cprob at maxTime: " + new Vector(computeConditionalProbabilities(maxTime)) + "\n");
        return sb.toString();
    }

    private EigenDecomposition getEigenDecomposition() {
        if (eigenDecomposition == null) {
            eigenDecomposition = eigenSystem.decomposeMatrix(squareMatrix(Q));
        }
        return eigenDecomposition;
    }

    private EigenDecomposition eigenDecomposition;

    public double[] computeConditionalProbabilities(double distance) {

        double[] matrix = new double[dim * dim];
        eigenSystem.computeExponential(getEigenDecomposition(), distance, matrix);

        return matrix;
    }

    public double computeConditionalProbability(double distance, int i, int j) {
        return eigenSystem.computeExponential(getEigenDecomposition(), distance, i, j);
    }

    private final double[] Q;
    private final double[] r;

    private final double lambda;
    private final double[] P;
    private final int phi;
    private final int dim;
    private final double epsilon;

    private final EigenSystem eigenSystem;

    private double maxTime;
}
