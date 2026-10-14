/*
 * MascotSurrogateOracleTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent.mascot;

import junit.framework.TestCase;

/** Checks the independent Pi_h reference before it is used to judge the stochastic mapper. */
public class MascotSurrogateOracleTest extends TestCase {

    /** The oracle's re-implemented forward pass reproduces MASCOT's log likelihood. */
    public void testLogLikelihoodMatchesMascot() throws Exception {
        for (MascotMappingFixtures.Spec spec : MascotMappingFixtures.oracleCases()) {
            MascotMappingFixtures.Built built = spec.build(0.01);
            double mascot = built.likelihood.getLogLikelihood();
            double oracle = built.oracle().logLikelihood;
            assertEquals(spec.name, mascot, oracle, 1.0e-11 * (1.0 + Math.abs(mascot)));
        }
        MascotMappingFixtures.Built highPressure = MascotMappingFixtures.highPressure().build(0.001);
        assertEquals(-73.73126497724, highPressure.likelihood.getLogLikelihood(), 1.0e-7);
        assertEquals(highPressure.likelihood.getLogLikelihood(), highPressure.oracle().logLikelihood, 1.0e-9);
    }

    /** Marginals are distributions, dwell times fill each branch, and the joint's margins are the marginals. */
    public void testInternalConsistency() throws Exception {
        for (MascotMappingFixtures.Spec spec : MascotMappingFixtures.oracleCases()) {
            MascotMappingFixtures.Built built = spec.build(0.02);
            MascotSurrogateOracle oracle = built.oracle();
            double[][] marginals = oracle.nodeMarginals();
            double[][] dwell = oracle.expectedDwell();
            for (int v = 0; v < oracle.nodeCount; v++) {
                double total = 0.0;
                for (int s = 0; s < spec.K; s++) {
                    assertTrue(spec.name + " marginal in [0,1]", marginals[v][s] >= -1e-15 && marginals[v][s] <= 1 + 1e-15);
                    total += marginals[v][s];
                }
                assertEquals(spec.name + " marginal sum", 1.0, total, 1.0e-12);
                if (v != oracle.root) {
                    double length = oracle.getHeight(oracle.getParent(v)) - oracle.getHeight(v);
                    double sum = 0.0;
                    for (int s = 0; s < spec.K; s++) {
                        sum += dwell[v][s];
                    }
                    assertEquals(spec.name + " dwell sum", length, sum, 1.0e-9 * (1.0 + length));
                }
            }

            int[] order = internalOrder(built);
            double[] joint = oracle.internalJoint(order);
            for (int q = 0; q < order.length; q++) {
                double[] margin = new double[spec.K];
                for (int x = 0; x < joint.length; x++) {
                    margin[digit(x, q, spec.K)] += joint[x];
                }
                for (int s = 0; s < spec.K; s++) {
                    assertEquals(spec.name + " joint margin, node " + order[q], marginals[order[q]][s], margin[s],
                            1.0e-10);
                }
            }
        }
    }

    /**
     * With one mobile lineage and every other lineage frozen in a known state, Pi_h equals the exact structured
     * coalescent at any step size, and MASCOT's likelihood equals the exact one up to RK4 error.
     */
    public void testSingleMobileLineageEqualsExact() throws Exception {
        MascotMappingFixtures.Built built = MascotMappingFixtures.singleMobile().build(0.05);
        MascotSurrogateOracle oracle = built.oracle();
        ExactStructuredCoalescent exact = built.exact();

        assertEquals(exact.logLikelihood(), built.likelihood.getLogLikelihood(), 1.0e-6);

        double[][] marginals = oracle.nodeMarginals();
        double[][] exactMarginals = exact.marginals();
        double[][] jumps = oracle.expectedJumps();
        double[][] exactJumps = exact.expectedJumps();
        double[][] dwell = oracle.expectedDwell();
        double[][] exactDwell = exact.expectedDwell();
        for (int v = 0; v < oracle.nodeCount; v++) {
            for (int s = 0; s < 3; s++) {
                assertEquals("marginal " + v + "," + s, exactMarginals[v][s], marginals[v][s], 1.0e-10);
            }
            if (v == oracle.root) {
                continue;
            }
            for (int i = 0; i < 9; i++) {
                assertEquals("jumps " + v + "," + i, exactJumps[v][i], jumps[v][i], 1.0e-10);
            }
            for (int s = 0; s < 3; s++) {
                assertEquals("dwell " + v + "," + s, exactDwell[v][s], dwell[v][s], 1.0e-10);
            }
        }
    }

    /** Pi_h converges as maxStep shrinks: successive differences in the node marginals decrease. */
    public void testConvergesInStepSize() throws Exception {
        for (MascotMappingFixtures.Spec spec : MascotMappingFixtures.oracleCases()) {
            MascotMappingFixtures.Built built = spec.build(0.04);
            double[][] coarse = built.oracle(0.04).nodeMarginals();
            double[][] medium = built.oracle(0.02).nodeMarginals();
            double[][] fine = built.oracle(0.01).nodeMarginals();
            double first = maxDifference(coarse, medium);
            double second = maxDifference(medium, fine);
            assertTrue(spec.name + ": " + first + " then " + second, second <= 0.5 * first + 1.0e-12);
        }
    }

    static int[] internalOrder(MascotMappingFixtures.Built built) {
        int[] order = new int[built.tree.getInternalNodeCount()];
        for (int i = 0; i < order.length; i++) {
            order[i] = built.tree.getInternalNode(i).getNumber();
        }
        java.util.Arrays.sort(order);
        return order;
    }

    static int digit(int index, int position, int K) {
        for (int q = 0; q < position; q++) {
            index /= K;
        }
        return index % K;
    }

    private static double maxDifference(double[][] a, double[][] b) {
        double max = 0.0;
        for (int v = 0; v < a.length; v++) {
            for (int s = 0; s < a[v].length; s++) {
                max = Math.max(max, Math.abs(a[v][s] - b[v][s]));
            }
        }
        return max;
    }
}
