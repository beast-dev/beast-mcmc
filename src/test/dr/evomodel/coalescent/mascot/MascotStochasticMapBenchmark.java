/*
 * MascotStochasticMapBenchmark.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent.mascot;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.io.NewickImporter;
import dr.evomodel.coalescent.basta.ConstantPopulationSizeModel;
import dr.evomodel.coalescent.mascot.MascotLikelihood;
import dr.evomodel.coalescent.mascot.MascotStochasticMap;
import dr.evomodel.coalescent.mascot.StochasticMapSampler;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.markovjumps.StateHistory;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;

/**
 * Speed and accuracy of the two MASCOT stochastic-mapping samplers (plan M4). Not a unit test.
 * <pre>
 * java -Xmx16g -cp "build:lib/*" test.dr.evomodel.coalescent.mascot.MascotStochasticMapBenchmark [speed] [accuracy]
 * </pre>
 */
public class MascotStochasticMapBenchmark {

    private static final StochasticMapSampler[] SAMPLERS = {StochasticMapSampler.GRID, StochasticMapSampler.H_TRANSFORM};

    public static void main(String[] args) throws Exception {
        List<String> parts = Arrays.asList(args);
        boolean all = parts.isEmpty();
        if (all || parts.contains("speed")) {
            speed();
        }
        if (all || parts.contains("accuracy")) {
            accuracy();
        }
    }

    // ---------------------------------------------------------------------------------------------------------
    // Speed
    // ---------------------------------------------------------------------------------------------------------

    private static void speed() throws Exception {
        MathUtils.setSeed(1);
        run(200, 5, 0.01, true, 11, false);   // JIT warm-up
        System.out.println("| tips | K | maxStep | tips sampled | cells | logL ms | grid prepare+draw ms | grid draw ms | " +
                "h-transform prepare+draw ms | h-transform draw ms | grid / h-transform (draw) |");
        System.out.println("|---|---|---|---|---|---|---|---|---|---|---|");
        for (int n : new int[]{200, 500, 1000}) {
            for (int K : new int[]{3, 5, 10, 20}) {
                for (double maxStep : new double[]{0.01, 0.001}) {
                    for (boolean heterochronous : new boolean[]{false, true}) {
                        try {
                            run(n, K, maxStep, heterochronous, 1000L * n + K, true);
                        } catch (OutOfMemoryError e) {
                            System.out.printf(Locale.US, "| %d | %d | %.3f | %s | out of memory |%n", n, K, maxStep,
                                    heterochronous ? "serial" : "together");
                        }
                    }
                }
            }
        }
    }

    private static void run(int n, int K, double maxStep, boolean heterochronous, long seed, boolean print)
            throws Exception {
        MascotLikelihood likelihood = randomLikelihood(n, K, maxStep, heterochronous, seed);
        double logL = timeMillis(() -> {
            likelihood.makeDirty();
            likelihood.getLogLikelihood();
        }, 3);
        double[] prepare = new double[2];
        double[] draw = new double[2];
        int cells = 0;
        for (int i = 0; i < 2; i++) {
            likelihood.setStochasticMapSampler(SAMPLERS[i]);
            prepare[i] = timeMillis(() -> {
                likelihood.makeDirty();
                likelihood.getLogLikelihood();
                likelihood.getStochasticMap();
            }, 3) - logL;
            int repeats = Math.max(3, (int) Math.min(50, 2000.0 / Math.max(prepare[i], 1.0)));
            draw[i] = timeMillis(() -> {
                likelihood.forceRedraw();
                likelihood.getStochasticMap();
            }, repeats);
            cells = likelihood.getStochasticMap().getCellCount();
        }
        if (print) {
            System.out.printf(Locale.US, "| %d | %d | %.3f | %s | %d | %.1f | %.1f | %.2f | %.1f | %.2f | %.1f |%n",
                    n, K, maxStep, heterochronous ? "serial" : "together", cells, logL, prepare[0], draw[0],
                    prepare[1], draw[1], draw[0] / draw[1]);
        }
    }

    private interface Task {
        void run() throws Exception;
    }

    /** Median wall time of {@code repeats} runs. */
    private static double timeMillis(Task task, int repeats) throws Exception {
        double[] times = new double[repeats];
        for (int r = 0; r < repeats; r++) {
            long start = System.nanoTime();
            task.run();
            times[r] = (System.nanoTime() - start) / 1.0e6;
        }
        Arrays.sort(times);
        return times[repeats / 2];
    }

    /** A random coalescent tree (pair rate 1/3) with random tip locations, migration rates and population sizes. */
    static MascotLikelihood randomLikelihood(int n, int K, double maxStep, boolean heterochronous, long seed)
            throws Exception {
        Random random = new Random(seed);
        double[] tipHeights = new double[n];
        for (int i = 0; i < n; i++) {
            tipHeights[i] = heterochronous ? 2.0 * random.nextDouble() : 0.0;
        }
        Integer[] order = new Integer[n];
        for (int i = 0; i < n; i++) {
            order[i] = i;
        }
        Arrays.sort(order, (a, b) -> Double.compare(tipHeights[a], tipHeights[b]));
        List<Integer> active = new ArrayList<>();
        Map<Integer, String> newick = new HashMap<>();
        Map<Integer, Double> height = new HashMap<>();
        int next = 0;
        int id = n;
        double t = 0.0;
        while (next < n || active.size() > 1) {
            int k = active.size();
            double rate = k * (k - 1) / 2.0 / 3.0;
            double wait = rate > 0 ? -Math.log(random.nextDouble()) / rate : Double.POSITIVE_INFINITY;
            double nextSample = next < n ? tipHeights[order[next]] : Double.POSITIVE_INFINITY;
            if (t + wait < nextSample) {
                t += wait;
                int a = active.remove(random.nextInt(active.size()));
                int b = active.remove(random.nextInt(active.size()));
                newick.put(id, "(" + newick.get(a) + ":" + (t - height.get(a)) + "," + newick.get(b) + ":" +
                        (t - height.get(b)) + ")");
                height.put(id, t);
                active.add(id++);
            } else {
                t = nextSample;
                int tip = order[next++];
                newick.put(tip, "T" + tip);
                height.put(tip, tipHeights[tip]);
                active.add(tip);
            }
        }
        String[] codes = new String[K];
        for (int i = 0; i < K; i++) {
            codes[i] = Integer.toString(i);
        }
        GeneralDataType dataType = new GeneralDataType(codes);
        TreeModel tree = new DefaultTreeModel("treeModel",
                new NewickImporter(newick.get(active.get(0)) + ";").importTree(null));
        Patterns patterns = new Patterns(dataType, tree);
        int[] pattern = new int[tree.getTaxonCount()];
        for (int i = 0; i < pattern.length; i++) {
            pattern[i] = random.nextInt(K);
        }
        patterns.addPattern(pattern);
        double[] rates = new double[K * (K - 1)];
        for (int i = 0; i < rates.length; i++) {
            rates[i] = (0.1 + 0.9 * random.nextDouble()) * 4.0 / K;
        }
        double[] sizes = new double[K];
        for (int i = 0; i < K; i++) {
            sizes[i] = (1.0 + 4.0 * random.nextDouble()) * 3.0 / K;
        }
        return new MascotLikelihood(null, tree, patterns, new Parameter.Default("migrationRates", rates), null, null,
                new ConstantPopulationSizeModel("constant", new Parameter.Default("popSizes", sizes), K, -1), K, null,
                maxStep, false, MascotLikelihood.DEFAULT_ANCESTRAL_STATE_TAG_NAME,
                MascotLikelihood.DEFAULT_MODE_STATE_TAG_NAME);
    }

    // ---------------------------------------------------------------------------------------------------------
    // Accuracy
    // ---------------------------------------------------------------------------------------------------------

    private static final int ACCURACY_DRAWS = 40000;
    private static final int REFERENCE_REFINEMENT = 10;

    private static void accuracy() throws Exception {
        System.out.println();
        System.out.println("Reference: Pi_h at maxStep / " + REFERENCE_REFINEMENT + " (a proxy for the continuous " +
                "surrogate); " + ACCURACY_DRAWS + " draws per row. 'own' = the sampler's exact target where it has " +
                "one (grid: Pi_h at maxStep).");
        System.out.println("| case | maxStep | max mu h | Pi_h vs reference (max abs, deterministic) | sampler | " +
                "node max abs z | node max abs diff | jumps max abs z | jumps max abs diff | dwell max abs z | " +
                "draw ms |");
        System.out.println("|---|---|---|---|---|---|---|---|---|---|---|");
        List<Object[]> cases = new ArrayList<>();
        for (MascotMappingFixtures.Spec spec : MascotMappingFixtures.oracleCases()) {
            cases.add(new Object[]{spec, 0.1});
            cases.add(new Object[]{spec, 0.01});
        }
        cases.add(new Object[]{MascotMappingFixtures.c2(), 0.3});   // stiff: h * stiffness near RK4's limit
        long seed = 1;
        for (Object[] row : cases) {
            MascotMappingFixtures.Spec spec = (MascotMappingFixtures.Spec) row[0];
            double maxStep = (Double) row[1];
            MascotMappingFixtures.Built built = spec.build(maxStep);
            MascotSurrogateOracle reference = built.oracle(maxStep / REFERENCE_REFINEMENT);
            MascotSurrogateOracle own = built.oracle(maxStep);
            double[][] refMarginals = reference.nodeMarginals();
            double[][] refJumps = reference.expectedJumps();
            double[][] refDwell = reference.expectedDwell();
            double gap = Math.max(maxAbs(own.nodeMarginals(), refMarginals), maxAbs(own.expectedJumps(), refJumps));
            for (StochasticMapSampler sampler : SAMPLERS) {
                built.likelihood.setStochasticMapSampler(sampler);
                MathUtils.setSeed(seed++);
                Moments moments = Moments.draw(built, ACCURACY_DRAWS);
                System.out.printf(Locale.US, "| %s | %.2f | %.2f | %.1e | %s | %.1f | %.4f | %.1f | %.4f | %.1f | %.3f |%n",
                        spec.name, maxStep, moments.maxRateTimesLength, gap, sampler,
                        moments.maxZ(moments.states, refMarginals), moments.maxDiff(moments.states, refMarginals),
                        moments.maxZ(moments.jumps, refJumps), moments.maxDiff(moments.jumps, refJumps),
                        moments.maxZ(moments.dwell, refDwell), moments.drawMillis);
            }
        }
    }

    private static double maxAbs(double[][] a, double[][] b) {
        double max = 0.0;
        for (int v = 0; v < a.length; v++) {
            if (a[v] == null) {
                continue;
            }
            for (int i = 0; i < a[v].length; i++) {
                max = Math.max(max, Math.abs(a[v][i] - b[v][i]));
            }
        }
        return max;
    }

    /** Running sums and sums of squares of per-node indicators, per-branch jump counts and dwell times. */
    private static final class Moments {
        int n;
        double[][][] states;   // [node][K][sum, sumSq]
        double[][][] jumps;    // [node][K K][sum, sumSq]
        double[][][] dwell;    // [node][K][sum, sumSq]
        double maxRateTimesLength;
        double drawMillis;

        static Moments draw(MascotMappingFixtures.Built built, int draws) {
            MascotLikelihood likelihood = built.likelihood;
            int nodes = built.tree.getNodeCount();
            int K = built.spec.K;
            int root = built.tree.getRoot().getNumber();
            Moments m = new Moments();
            m.n = draws;
            m.states = new double[nodes][K][2];
            m.jumps = new double[nodes][K * K][2];
            m.dwell = new double[nodes][K][2];
            long start = System.nanoTime();
            for (int d = 0; d < draws; d++) {
                likelihood.forceRedraw();
                MascotStochasticMap map = likelihood.getStochasticMap();
                m.maxRateTimesLength = map.getMaxRateTimesLength();
                for (int v = 0; v < nodes; v++) {
                    int state = map.getNodeState(v);
                    for (int s = 0; s < K; s++) {
                        double indicator = s == state ? 1.0 : 0.0;
                        m.states[v][s][0] += indicator;
                        m.states[v][s][1] += indicator;
                    }
                    if (v == root) {
                        continue;
                    }
                    StateHistory history = map.getBranchHistory(v);
                    int[] counts = history.getJumpCounts();
                    for (int k = 0; k < K * K; k++) {
                        m.jumps[v][k][0] += counts[k];
                        m.jumps[v][k][1] += counts[k] * counts[k];
                    }
                    double[] times = history.getWaitingTimes();
                    for (int s = 0; s < K; s++) {
                        m.dwell[v][s][0] += times[s];
                        m.dwell[v][s][1] += times[s] * times[s];
                    }
                }
            }
            m.drawMillis = (System.nanoTime() - start) / 1.0e6 / draws;
            return m;
        }

        double maxZ(double[][][] sums, double[][] expected) {
            double max = 0.0;
            for (int v = 0; v < sums.length; v++) {
                if (expected[v] == null) {
                    continue;
                }
                for (int i = 0; i < expected[v].length; i++) {
                    double mean = sums[v][i][0] / n;
                    double variance = Math.max(sums[v][i][1] / n - mean * mean, 0.0);
                    if (variance > 0.0) {
                        max = Math.max(max, Math.abs(mean - expected[v][i]) / Math.sqrt(variance / n));
                    }
                }
            }
            return max;
        }

        double maxDiff(double[][][] sums, double[][] expected) {
            double max = 0.0;
            for (int v = 0; v < sums.length; v++) {
                if (expected[v] == null) {
                    continue;
                }
                for (int i = 0; i < expected[v].length; i++) {
                    max = Math.max(max, Math.abs(sums[v][i][0] / n - expected[v][i]));
                }
            }
            return max;
        }
    }
}
