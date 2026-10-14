/*
 * KilledUniformizedBridgeTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent.mascot;

import dr.evomodel.coalescent.mascot.KilledUniformizedBridge;
import dr.inference.markovjumps.StateChange;
import dr.inference.markovjumps.StateHistory;
import dr.math.MathUtils;
import junit.framework.TestCase;

import java.util.List;

/** One-cell killed bridges against matrix exponentials and closed-form (Van Loan) endpoint-conditioned moments. */
public class KilledUniformizedBridgeTest extends TestCase {

    private static final double Z_LIMIT = 5.0;

    public void testTransitionMatchesPadeExponential() {
        MathUtils.setSeed(11);
        for (int K : new int[]{2, 3, 5}) {
            for (double target : new double[]{0.01, 0.5, 5.0, 19.0}) {
                double[] G = randomKilledGenerator(K, true);
                double mu = maxExitRate(G, K);
                double h = target / mu;
                KilledUniformizedBridge bridge = new KilledUniformizedBridge(K);
                bridge.setGenerator(G, h);
                double[][] T = MatrixExp.expm(toMatrix(G, K, h));

                double scale = 0.0;
                for (double[] row : T) {
                    for (double v : row) {
                        scale = Math.max(scale, Math.abs(v));
                    }
                }
                double tolerance = 1.0e-13 * Math.max(1.0, scale);
                for (int s = 0; s < K; s++) {
                    double[] unit = new double[K];
                    unit[s] = 1.0;
                    double[] row = new double[K];
                    double[] column = new double[K];
                    bridge.transitionRow(unit, 0, row, 0);
                    bridge.transitionColumn(unit, 0, column, 0);
                    for (int t = 0; t < K; t++) {
                        assertEquals("row K=" + K + " mu h=" + target, T[s][t], row[t], tolerance);
                        assertEquals("column K=" + K + " mu h=" + target, T[t][s], column[t], tolerance);
                        assertEquals("entry", T[s][t], bridge.transitionProbability(s, t), tolerance);
                    }
                }
            }
        }
    }

    /** States several jumps away keep a positive transition probability even in very short cells. */
    public void testSupportIsKeptInShortCells() {
        int K = 6;
        double[] G = new double[K * K];
        for (int s = 0; s + 1 < K; s++) {
            G[s * K + s + 1] = 1.0;
            G[s * K + s] = -1.0 - 0.5;
        }
        G[(K - 1) * K + K - 1] = -0.5;
        KilledUniformizedBridge bridge = new KilledUniformizedBridge(K);
        bridge.setGenerator(G, 1.0e-3);
        double far = bridge.transitionProbability(0, K - 1);
        double leading = Math.pow(1.0e-3, K - 1) / 120.0;   // h^5 / 5! to leading order
        assertTrue("T[0][K-1] = " + far, far > 0.0);
        assertEquals(1.0, far / leading, 1.0e-2);
    }

    public void testZeroGeneratorIsIdentity() {
        MathUtils.setSeed(3);
        int K = 3;
        KilledUniformizedBridge bridge = new KilledUniformizedBridge(K);
        bridge.setGenerator(new double[K * K], 2.0);
        for (int s = 0; s < K; s++) {
            for (int t = 0; t < K; t++) {
                assertEquals(s == t ? 1.0 : 0.0, bridge.transitionProbability(s, t), 0.0);
            }
            assertEquals(0, bridge.sampleBridge(s, s, 0.0, 2.0).getNumberOfJumps());
        }
        KilledUniformizedBridge.Path path = new KilledUniformizedBridge.Path();
        assertEquals(2, bridge.sampleStartAndPath(new double[]{0.2, 0.3, 0.5}, 0, 2, 0.0, 2.0, path));
        assertEquals(0, path.size());
    }

    public void testRejectsInvalidGenerators() {
        KilledUniformizedBridge bridge = new KilledUniformizedBridge(2);
        try {
            bridge.setGenerator(new double[]{-1.0, -0.5, 0.5, -0.5}, 1.0);
            fail("negative off-diagonal accepted");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        try {
            bridge.setGenerator(new double[]{-1.0, 2.0, 0.5, -0.5}, 1.0);
            fail("positive row sum accepted");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        KilledUniformizedBridge fresh = new KilledUniformizedBridge(2);
        try {
            fresh.transitionProbability(0, 0);
            fail("used before setGenerator");
        } catch (IllegalStateException expected) {
            // expected
        }
    }

    /** The start state is drawn in proportion to weight[s] T[s][end]. */
    public void testStartStateLaw() {
        MathUtils.setSeed(17);
        int K = 3;
        double[] G = randomKilledGenerator(K, false);
        double h = 1.5 / maxExitRate(G, K);
        KilledUniformizedBridge bridge = new KilledUniformizedBridge(K);
        bridge.setGenerator(G, h);
        double[] weights = {0.5, 0.2, 0.3};
        int end = 1;
        double[] expected = new double[K];
        double total = 0.0;
        for (int s = 0; s < K; s++) {
            expected[s] = weights[s] * bridge.transitionProbability(s, end);
            total += expected[s];
        }
        int n = 100000;
        int[] counts = new int[K];
        KilledUniformizedBridge.Path path = new KilledUniformizedBridge.Path();
        for (int i = 0; i < n; i++) {
            path.clear();
            counts[bridge.sampleStartAndPath(weights, 0, end, 0.0, h, path)]++;
        }
        for (int s = 0; s < K; s++) {
            double p = expected[s] / total;
            assertProportion("start state " + s, p, counts[s], n);
        }
    }

    /**
     * Endpoint-conditioned moments of the bridge against closed forms: expected labelled jumps and dwell times
     * (Van Loan integrals), P(no jump), the count generating function at 1/2, and jumps in the first third of the
     * cell (which checks where the jumps are placed, not only how many there are).
     */
    public void testBridgeMomentsMatchClosedForms() {
        MathUtils.setSeed(23);
        int K = 3;
        double[] G = randomKilledGenerator(K, false);
        double h = 2.0 / maxExitRate(G, K);
        KilledUniformizedBridge bridge = new KilledUniformizedBridge(K);
        bridge.setGenerator(G, h);
        double[][] generator = toMatrix(G, K, 1.0);
        double[][] T = MatrixExp.expm(toMatrix(G, K, h));
        double[][] offDiagonal = new double[K][K];
        for (int a = 0; a < K; a++) {
            for (int b = 0; b < K; b++) {
                if (a != b) {
                    offDiagonal[a][b] = generator[a][b];
                }
            }
        }
        double[][] earlyJumps = MatrixExp.mul(vanLoan(generator, offDiagonal, h / 3.0, K),
                MatrixExp.expm(toMatrix(G, K, 2.0 * h / 3.0)));
        double[][] halfGenerator = new double[K][K];
        for (int a = 0; a < K; a++) {
            for (int b = 0; b < K; b++) {
                halfGenerator[a][b] = (a == b ? 1.0 : 0.5) * generator[a][b] * h;
            }
        }
        double[][] halfPgf = MatrixExp.expm(halfGenerator);

        int n = 40000;
        for (int[] ends : new int[][]{{0, 0}, {0, 2}, {2, 1}}) {
            int s = ends[0];
            int r = ends[1];
            double normalizer = T[s][r];
            double[] sumCounts = new double[K * K];
            double[] sumSqCounts = new double[K * K];
            double[] sumDwell = new double[K];
            double[] sumSqDwell = new double[K];
            int noJump = 0;
            double sumEarly = 0.0;
            double sumSqEarly = 0.0;
            double sumPgf = 0.0;
            double sumSqPgf = 0.0;
            for (int i = 0; i < n; i++) {
                StateHistory history = bridge.sampleBridge(s, r, 0.0, h);
                int[] counts = history.getJumpCounts();
                double[] dwell = history.getWaitingTimes();
                for (int k = 0; k < K * K; k++) {
                    sumCounts[k] += counts[k];
                    sumSqCounts[k] += counts[k] * counts[k];
                }
                for (int a = 0; a < K; a++) {
                    sumDwell[a] += dwell[a];
                    sumSqDwell[a] += dwell[a] * dwell[a];
                }
                int jumps = history.getNumberOfJumps();
                if (jumps == 0) {
                    noJump++;
                }
                int early = 0;
                List<StateChange> changes = history.getStateChanges();
                for (int j = 1; j < changes.size() - 1; j++) {
                    if (changes.get(j).getTime() < h / 3.0) {
                        early++;
                    }
                }
                sumEarly += early;
                sumSqEarly += early * early;
                double pgf = Math.pow(0.5, jumps);
                sumPgf += pgf;
                sumSqPgf += pgf * pgf;
            }
            String label = "bridge " + s + "->" + r;
            for (int a = 0; a < K; a++) {
                for (int b = 0; b < K; b++) {
                    if (a == b) {
                        continue;
                    }
                    double[][] register = new double[K][K];
                    register[a][b] = generator[a][b];
                    double expected = vanLoan(generator, register, h, K)[s][r] / normalizer;
                    assertMean(label + " jumps " + a + "->" + b, expected, sumCounts[a * K + b],
                            sumSqCounts[a * K + b], n);
                }
                double[][] register = new double[K][K];
                register[a][a] = 1.0;
                assertMean(label + " dwell " + a, vanLoan(generator, register, h, K)[s][r] / normalizer,
                        sumDwell[a], sumSqDwell[a], n);
            }
            double noJumpProbability = s == r ? Math.exp(generator[s][s] * h) / normalizer : 0.0;
            assertProportion(label + " no jump", noJumpProbability, noJump, n);
            assertMean(label + " early jumps", earlyJumps[s][r] / normalizer, sumEarly, sumSqEarly, n);
            assertMean(label + " E[0.5^N]", halfPgf[s][r] / normalizer, sumPgf, sumSqPgf, n);
        }
    }

    static double[] randomKilledGenerator(int K, boolean withZeros) {
        double[] G = new double[K * K];
        for (int a = 0; a < K; a++) {
            double rowSum = 0.0;
            for (int b = 0; b < K; b++) {
                if (a == b) {
                    continue;
                }
                double rate = withZeros && MathUtils.nextDouble() < 0.25 ? 0.0 : 0.1 + 1.9 * MathUtils.nextDouble();
                G[a * K + b] = rate;
                rowSum += rate;
            }
            double killing = 3.0 * MathUtils.nextDouble();
            G[a * K + a] = -rowSum - killing;
        }
        return G;
    }

    static double maxExitRate(double[] G, int K) {
        double mu = 0.0;
        for (int s = 0; s < K; s++) {
            mu = Math.max(mu, -G[s * K + s]);
        }
        return mu;
    }

    static double[][] toMatrix(double[] G, int K, double scale) {
        double[][] out = new double[K][K];
        for (int a = 0; a < K; a++) {
            for (int b = 0; b < K; b++) {
                out[a][b] = G[a * K + b] * scale;
            }
        }
        return out;
    }

    /** int_0^t exp(G u) B exp(G (t - u)) du. */
    static double[][] vanLoan(double[][] G, double[][] B, double t, int K) {
        double[][] block = new double[2 * K][2 * K];
        for (int a = 0; a < K; a++) {
            for (int b = 0; b < K; b++) {
                block[a][b] = G[a][b] * t;
                block[K + a][K + b] = G[a][b] * t;
                block[a][K + b] = B[a][b] * t;
            }
        }
        double[][] exponential = MatrixExp.expm(block);
        double[][] out = new double[K][K];
        for (int a = 0; a < K; a++) {
            for (int b = 0; b < K; b++) {
                out[a][b] = exponential[a][K + b];
            }
        }
        return out;
    }

    static void assertMean(String label, double expected, double sum, double sumSquares, int n) {
        double mean = sum / n;
        double variance = Math.max(0.0, sumSquares / n - mean * mean);
        if (variance == 0.0) {
            // Every draw gave the same value, as for a rare count that never occurred: missing n |mean - expected|
            // events of expected rate has probability about exp(-n |mean - expected|), so require that it is not
            // implausibly small (below about 1e-6).
            assertTrue(label + ": every draw gave " + mean + ", expected mean " + expected,
                    n * Math.abs(mean - expected) < 14.0);
            return;
        }
        double z = (mean - expected) / Math.sqrt(variance / n);
        assertTrue(label + ": mean " + mean + ", expected " + expected + ", z " + z, Math.abs(z) < Z_LIMIT);
    }

    static void assertProportion(String label, double p, int count, int n) {
        if (p < 1.0e-12) {
            assertEquals(label + " (impossible event observed)", 0, count);
            return;
        }
        if (p > 1.0 - 1.0e-12) {
            assertEquals(label + " (certain event missed)", n, count);
            return;
        }
        double z = (count - n * p) / Math.sqrt(n * p * (1.0 - p));
        assertTrue(label + ": observed " + ((double) count / n) + ", expected " + p + ", z " + z,
                Math.abs(z) < Z_LIMIT);
    }
}
