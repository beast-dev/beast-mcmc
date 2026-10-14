/*
 * MascotStochasticMapTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent.mascot;

import dr.evolution.tree.NodeRef;
import dr.evomodel.coalescent.mascot.MascotLikelihood;
import dr.evomodel.coalescent.mascot.MascotStochasticMap;
import dr.evomodel.tree.TreeModel;
import dr.inference.markovjumps.StateChange;
import dr.inference.markovjumps.StateHistory;
import dr.math.MathUtils;
import dr.math.distributions.ChiSquareDistribution;
import junit.framework.TestCase;

import java.util.List;

/**
 * MASCOT stochastic mapping: every branch gets a complete {@link StateHistory}, and node states and paths form one
 * joint draw from Pi_h. Draws are compared with the independent {@link MascotSurrogateOracle} and, where Pi_h is
 * exact, with the brute-force {@link ExactStructuredCoalescent}.
 */
public class MascotStochasticMapTest extends TestCase {

    private static final double Z_LIMIT = 5.0;
    private static final double EARLY_FRACTION = 1.0 / 3.0;

    public void testEveryHistoryIsCompleteAndConsistent() throws Exception {
        MathUtils.setSeed(101);
        for (MascotMappingFixtures.Spec spec : MascotMappingFixtures.oracleCases()) {
            checkStructure(spec.build(0.01), 200);
            checkStructure(spec.build(0.25), 200);
        }
        checkStructure(MascotMappingFixtures.singleMobile().build(0.05), 200);
        checkStructure(MascotMappingFixtures.highPressure().build(0.01), 50);
    }

    public void testC1MatchesOracle() throws Exception {
        checkAgainstOracle(MascotMappingFixtures.c1().build(0.01), 20000, 1);
    }

    public void testC2MatchesOracle() throws Exception {
        checkAgainstOracle(MascotMappingFixtures.c2().build(0.01), 20000, 2);
    }

    public void testC3AmbiguousTipMatchesOracle() throws Exception {
        checkAgainstOracle(MascotMappingFixtures.c3().build(0.02), 20000, 3);
    }

    public void testC4EpochsMatchOracle() throws Exception {
        checkAgainstOracle(MascotMappingFixtures.c4().build(0.01), 20000, 4);
    }

    public void testC5StrictClockMatchesOracle() throws Exception {
        checkAgainstOracle(MascotMappingFixtures.c5().build(0.01), 20000, 5);
    }

    public void testC6BranchRatesMatchOracle() throws Exception {
        checkAgainstOracle(MascotMappingFixtures.c6().build(0.02), 20000, 6);
    }

    /** At a coarse step the draws still follow Pi_h exactly (Pi_h itself is further from the continuum). */
    public void testCoarseStepMatchesOracle() throws Exception {
        checkAgainstOracle(MascotMappingFixtures.c2().build(0.3), 20000, 7);
    }

    /** One mobile lineage among frozen lineages in known states: Pi_h is the exact structured coalescent. */
    public void testSingleMobileLineageMatchesExact() throws Exception {
        MathUtils.setSeed(8);
        MascotMappingFixtures.Built built = MascotMappingFixtures.singleMobile().build(0.05);
        ExactStructuredCoalescent exact = built.exact();
        double[][] marginals = exact.marginals();
        double[][] jumps = exact.expectedJumps();
        double[][] dwell = exact.expectedDwell();
        Sample sample = Sample.draw(built, 20000, null, null);
        int K = built.spec.K;
        for (int v = 0; v < sample.nodeCount; v++) {
            for (int s = 0; s < K; s++) {
                KilledUniformizedBridgeTest.assertProportion("node " + v + " state " + s, marginals[v][s],
                        sample.stateCounts[v][s], sample.n);
            }
            if (v == sample.root) {
                continue;
            }
            for (int i = 0; i < K * K; i++) {
                if (i % (K + 1) != 0) {
                    KilledUniformizedBridgeTest.assertMean("jumps " + v + "," + i, jumps[v][i],
                            sample.sumJumps[v][i], sample.sumSqJumps[v][i], sample.n);
                }
            }
            for (int s = 0; s < K; s++) {
                KilledUniformizedBridgeTest.assertMean("dwell " + v + "," + s, dwell[v][s], sample.sumDwell[v][s],
                        sample.sumSqDwell[v][s], sample.n);
            }
        }
    }

    /** One draw per state: shared until the state is stored, the model changes, or a redraw is forced. */
    public void testRedrawPolicy() throws Exception {
        MathUtils.setSeed(9);
        MascotMappingFixtures.Built built = MascotMappingFixtures.c3().build(0.02);
        MascotLikelihood likelihood = built.likelihood;

        MascotStochasticMap first = likelihood.getStochasticMap();
        assertSame(first, likelihood.getStochasticMap());
        assertEquals(likelihood.getLogLikelihood(), first.getLogLikelihood(), 0.0);

        likelihood.forceRedraw();
        MascotStochasticMap redrawn = likelihood.getStochasticMap();
        assertNotSame(first, redrawn);

        likelihood.storeModelState();
        MascotStochasticMap afterStore = likelihood.getStochasticMap();
        assertNotSame("storing (once per MCMC iteration) starts a new log step", redrawn, afterStore);
        likelihood.acceptModelState();
        assertSame("accepting keeps the draw", afterStore, likelihood.getStochasticMap());

        double before = likelihood.getLogLikelihood();
        likelihood.storeModelState();
        built.popSizes.setParameterValue(0, 1.3 * built.popSizes.getParameterValue(0));
        MascotStochasticMap proposed = likelihood.getStochasticMap();
        assertEquals(likelihood.getLogLikelihood(), proposed.getLogLikelihood(), 0.0);
        assertTrue(Math.abs(proposed.getLogLikelihood() - before) > 1.0e-6);

        likelihood.restoreModelState();
        MascotStochasticMap restored = likelihood.getStochasticMap();
        assertNotSame(proposed, restored);
        assertEquals(before, restored.getLogLikelihood(), 0.0);
    }

    /** Histories are handed out as copies, so a caller mutating one cannot change what others read. */
    public void testHistoriesAreCopies() throws Exception {
        MathUtils.setSeed(10);
        MascotMappingFixtures.Built built = MascotMappingFixtures.c1().build(0.01);
        MascotStochasticMap map = built.likelihood.getStochasticMap();
        NodeRef tip = built.tree.getExternalNode(0);
        StateHistory a = map.getBranchHistory(tip.getNumber());
        String before = a.toStringChanges(1, built.dataType);
        a.rescaleTimesOfEvents(0.0, 100.0);
        assertEquals(before, map.getBranchHistory(tip.getNumber()).toStringChanges(1, built.dataType));
        assertEquals(before, built.likelihood.getBranchHistory(tip).toStringChanges(1, built.dataType));
    }

    // ---------------------------------------------------------------------------------------------------------

    static void checkStructure(MascotMappingFixtures.Built built, int draws) {
        MascotLikelihood likelihood = built.likelihood;
        TreeModel tree = built.tree;
        int K = built.spec.K;
        for (int d = 0; d < draws; d++) {
            likelihood.forceRedraw();
            MascotStochasticMap map = likelihood.getStochasticMap();
            assertEquals(likelihood.getLogLikelihood(), map.getLogLikelihood(), 0.0);
            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);
                int v = node.getNumber();
                int state = map.getNodeState(v);
                assertTrue(state >= 0 && state < K);
                if (tree.isExternal(node)) {
                    assertTrue(built.spec.name + ": tip " + v + " drawn outside its data",
                            built.tipPartials[v][state] > 0.0);
                }
                StateHistory history = map.getBranchHistory(v);
                if (tree.isRoot(node)) {
                    assertNull(history);
                    continue;
                }
                NodeRef parent = tree.getParent(node);
                double bottom = tree.getNodeHeight(node);
                double top = tree.getNodeHeight(parent);
                int parentState = map.getNodeState(parent.getNumber());

                List<StateChange> changes = history.getStateChanges();   // also requires a finalized history
                assertEquals(bottom, history.getStartingTime(), 0.0);
                assertEquals(top, history.getEndingTime(), 0.0);
                assertEquals(state, history.getStartingState());
                int current = state;
                double time = bottom;
                for (int j = 1; j < changes.size() - 1; j++) {
                    StateChange change = changes.get(j);
                    assertTrue(built.spec.name + ": change time not strictly increasing inside the branch",
                            change.getTime() > time && change.getTime() < top);
                    assertTrue(built.spec.name + ": self-change recorded", change.getState() != current);
                    current = change.getState();
                    time = change.getTime();
                }
                assertEquals(built.spec.name + ": replayed path does not end in the parent's state",
                        parentState, current);
                assertEquals(changes.size() - 2, history.getNumberOfJumps());

                double dwell = 0.0;
                for (double w : history.getWaitingTimes()) {
                    dwell += w;
                }
                assertEquals(top - bottom, dwell, 1.0e-12 * (1.0 + top));
                int counted = 0;
                for (int c : history.getJumpCounts()) {
                    counted += c;
                }
                assertEquals(history.getNumberOfJumps(), counted);
            }
        }
    }

    private static void checkAgainstOracle(MascotMappingFixtures.Built built, int draws, long seed) {
        checkAgainstOracle(built, draws, seed, built.oracle());
    }

    /** Compares draws from {@code built.likelihood} (whatever its sampler) with {@code oracle}'s Pi_h. */
    static void checkAgainstOracle(MascotMappingFixtures.Built built, int draws, long seed,
                                   MascotSurrogateOracle oracle) {
        MathUtils.setSeed(seed);
        int K = built.spec.K;
        int[] order = MascotSurrogateOracleTest.internalOrder(built);
        Sample sample = Sample.draw(built, draws, order, oracle);
        String name = built.spec.name + " (maxStep " + built.maxStep + ")";

        double[][] marginals = oracle.nodeMarginals();
        for (int v = 0; v < sample.nodeCount; v++) {
            for (int s = 0; s < K; s++) {
                KilledUniformizedBridgeTest.assertProportion(name + " node " + v + " state " + s, marginals[v][s],
                        sample.stateCounts[v][s], sample.n);
            }
        }
        assertJointMatches(name, oracle.internalJoint(order), sample.jointCounts, sample.n);

        double[][] jumps = oracle.expectedJumps();
        double[][] dwell = oracle.expectedDwell();
        double[] noJump = oracle.jumpCountGeneratingFunction(0.0);
        double[] half = oracle.jumpCountGeneratingFunction(0.5);
        double[] early = oracle.expectedEarlyJumps(EARLY_FRACTION);
        for (int v = 0; v < sample.nodeCount; v++) {
            if (v == sample.root) {
                continue;
            }
            String branch = name + " branch above " + v;
            for (int a = 0; a < K; a++) {
                for (int b = 0; b < K; b++) {
                    if (a != b) {
                        KilledUniformizedBridgeTest.assertMean(branch + " jumps " + a + "->" + b,
                                jumps[v][a * K + b], sample.sumJumps[v][a * K + b],
                                sample.sumSqJumps[v][a * K + b], sample.n);
                    }
                }
                KilledUniformizedBridgeTest.assertMean(branch + " dwell " + a, dwell[v][a], sample.sumDwell[v][a],
                        sample.sumSqDwell[v][a], sample.n);
            }
            KilledUniformizedBridgeTest.assertProportion(branch + " no jump", noJump[v], sample.noJumpCounts[v],
                    sample.n);
            KilledUniformizedBridgeTest.assertMean(branch + " E[0.5^N]", half[v], sample.sumHalfPower[v],
                    sample.sumSqHalfPower[v], sample.n);
            KilledUniformizedBridgeTest.assertMean(branch + " early jumps", early[v], sample.sumEarly[v],
                    sample.sumSqEarly[v], sample.n);
        }
    }

    /** Pearson chi-square over joint configurations, pooling those expected fewer than five times. */
    private static void assertJointMatches(String label, double[] expected, int[] counts, int n) {
        double chiSquare = 0.0;
        int cells = 0;
        double pooledExpected = 0.0;
        int pooledObserved = 0;
        for (int x = 0; x < expected.length; x++) {
            double e = n * expected[x];
            if (e < 5.0) {
                pooledExpected += e;
                pooledObserved += counts[x];
                if (expected[x] < 1.0e-12) {
                    assertEquals(label + " impossible joint configuration " + x + " observed", 0, counts[x]);
                }
            } else {
                chiSquare += (counts[x] - e) * (counts[x] - e) / e;
                cells++;
            }
        }
        if (pooledExpected >= 5.0) {
            chiSquare += (pooledObserved - pooledExpected) * (pooledObserved - pooledExpected) / pooledExpected;
            cells++;
        }
        if (cells < 2) {
            return;
        }
        double p = 1.0 - ChiSquareDistribution.cdf(chiSquare, cells - 1);
        assertTrue(label + " joint chi-square " + chiSquare + " on " + (cells - 1) + " df, p = " + p, p > 1.0e-5);
    }

    /** Sufficient statistics of repeated draws, all computed from the StateHistory objects. */
    private static final class Sample {
        int n;
        int nodeCount;
        int root;
        int[][] stateCounts;
        int[] jointCounts;
        double[][] sumJumps;
        double[][] sumSqJumps;
        double[][] sumDwell;
        double[][] sumSqDwell;
        int[] noJumpCounts;
        double[] sumHalfPower;
        double[] sumSqHalfPower;
        double[] sumEarly;
        double[] sumSqEarly;

        static Sample draw(MascotMappingFixtures.Built built, int draws, int[] jointOrder,
                           MascotSurrogateOracle oracle) {
            MascotLikelihood likelihood = built.likelihood;
            TreeModel tree = built.tree;
            int K = built.spec.K;
            Sample sample = new Sample();
            sample.n = draws;
            sample.nodeCount = tree.getNodeCount();
            sample.root = tree.getRoot().getNumber();
            int nodes = sample.nodeCount;
            sample.stateCounts = new int[nodes][K];
            sample.sumJumps = new double[nodes][K * K];
            sample.sumSqJumps = new double[nodes][K * K];
            sample.sumDwell = new double[nodes][K];
            sample.sumSqDwell = new double[nodes][K];
            sample.noJumpCounts = new int[nodes];
            sample.sumHalfPower = new double[nodes];
            sample.sumSqHalfPower = new double[nodes];
            sample.sumEarly = new double[nodes];
            sample.sumSqEarly = new double[nodes];
            double[][][] cells = new double[nodes][][];
            if (oracle != null) {
                int size = 1;
                for (int q = 0; q < jointOrder.length; q++) {
                    size *= K;
                }
                sample.jointCounts = new int[size];
                for (int v = 0; v < nodes; v++) {
                    if (v != sample.root) {
                        cells[v] = oracle.cellBoundaries(v);
                    }
                }
            }

            for (int d = 0; d < draws; d++) {
                likelihood.forceRedraw();
                MascotStochasticMap map = likelihood.getStochasticMap();
                if (jointOrder != null) {
                    int index = 0;
                    for (int q = jointOrder.length - 1; q >= 0; q--) {
                        index = index * K + map.getNodeState(jointOrder[q]);
                    }
                    sample.jointCounts[index]++;
                }
                for (int v = 0; v < nodes; v++) {
                    sample.stateCounts[v][map.getNodeState(v)]++;
                    if (v == sample.root) {
                        continue;
                    }
                    StateHistory history = map.getBranchHistory(v);
                    int[] counts = history.getJumpCounts();
                    for (int k = 0; k < K * K; k++) {
                        sample.sumJumps[v][k] += counts[k];
                        sample.sumSqJumps[v][k] += counts[k] * counts[k];
                    }
                    double[] dwell = history.getWaitingTimes();
                    for (int a = 0; a < K; a++) {
                        sample.sumDwell[v][a] += dwell[a];
                        sample.sumSqDwell[v][a] += dwell[a] * dwell[a];
                    }
                    int jumps = history.getNumberOfJumps();
                    if (jumps == 0) {
                        sample.noJumpCounts[v]++;
                    }
                    double halfPower = Math.pow(0.5, jumps);
                    sample.sumHalfPower[v] += halfPower;
                    sample.sumSqHalfPower[v] += halfPower * halfPower;
                    if (cells[v] != null) {
                        int earlyJumps = countEarlyJumps(history, cells[v]);
                        sample.sumEarly[v] += earlyJumps;
                        sample.sumSqEarly[v] += earlyJumps * earlyJumps;
                    }
                }
            }
            return sample;
        }

        private static int countEarlyJumps(StateHistory history, double[][] cells) {
            double[] starts = cells[0];
            double[] ends = cells[1];
            List<StateChange> changes = history.getStateChanges();
            int early = 0;
            int cell = 0;
            for (int j = 1; j < changes.size() - 1; j++) {
                double t = changes.get(j).getTime();
                while (cell + 1 < starts.length && t >= ends[cell]) {
                    cell++;
                }
                if ((t - starts[cell]) / (ends[cell] - starts[cell]) < EARLY_FRACTION) {
                    early++;
                }
            }
            return early;
        }
    }
}
