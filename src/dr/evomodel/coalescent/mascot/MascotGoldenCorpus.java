/*
 * MascotGoldenCorpus.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evomodel.coalescent.StructuredCoalescentSchedule;
import dr.evomodel.coalescent.StructuredCoalescentTipData;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Golden test cases for the native MASCOT kernel: inputs exactly as the native API receives them (schedule, tip
 * partials, migration matrices and inverse population sizes from {@link MascotRuntime#updateEpochRates}) and the
 * outputs of {@link GenericMascotLikelihoodDelegate}, the reference. Every double is written as the 16 hex digits of
 * its raw bits, so a reader reproduces it exactly.
 * <p>
 * Usage: {@code MascotGoldenCorpus <outputDirectory>} writes the generated cases (seeded, so the same every time).
 */
public final class MascotGoldenCorpus {

    public static final String FORMAT = "mascot-golden 1";

    /** One case: the inputs of an evaluation. */
    public static final class Case {
        public final String name;
        public final int stateCount;
        public final double maxStep;
        public final double[] boundaries;
        public final double[] theta;
        public final double[] branchRates;
        final MascotPreparedInput prepared;

        Case(String name, int stateCount, double maxStep, double[] boundaries, double[] theta, double[] branchRates,
             MascotPreparedInput prepared) {
            this.name = name;
            this.stateCount = stateCount;
            this.maxStep = maxStep;
            this.boundaries = boundaries;
            this.theta = theta;
            this.branchRates = branchRates;
            this.prepared = prepared;
        }
    }

    private MascotGoldenCorpus() {
    }

    public static void main(String[] args) throws IOException {
        File directory = new File(args[0]);
        if (!directory.isDirectory() && !directory.mkdirs()) {
            throw new IOException("cannot create " + directory);
        }
        for (Case c : generatedCases()) {
            write(c, new File(directory, c.name + ".golden"));
        }
    }

    /**
     * The generated cases:
     * <ul>
     * <li>a K sweep (2 to 128, including K not a multiple of 2 or 4) on one 150-tip heterochronous tree with branch
     * rates, at the production step 0.02;</li>
     * <li>small trees (2 to 80 tips, K = 2, 3, 5, 7) cycling through 1-3 epochs with boundaries on and between
     * events, simultaneous samples (zero-length intervals), zero migration rates, ambiguous tips, branch rates
     * absent or present, coarse and fine steps, and stiff (small Ne) regimes;</li>
     * <li>a coalescent-rate error case (two lineages in different demes with no migration).</li>
     * </ul>
     */
    public static List<Case> generatedCases() {
        List<Case> cases = new ArrayList<>();
        int[] sweep = {2, 3, 4, 5, 6, 7, 8, 9, 12, 16, 17, 24, 32, 33, 48, 63, 64, 128};
        for (int k : sweep) {
            cases.add(randomCase("sweepK" + k, 1000 + k, 150, k, 1, 0.02, true, false, false, false, 0.0));
        }
        int[] tipCounts = {2, 3, 5, 10, 20, 40, 80};
        int[] stateCounts = {2, 3, 5, 7};
        double[] steps = {0.02, 0.005, 0.1};
        int variant = 0;
        for (int n : tipCounts) {
            for (int k : stateCounts) {
                int epochs = 1 + variant % 3;
                double maxStep = steps[variant % steps.length];
                boolean branchRates = variant % 2 == 0;
                boolean simultaneous = variant % 4 == 1;
                boolean zeroRates = variant % 5 == 2;
                boolean ambiguous = variant % 3 == 2;
                double logNeShift = variant % 6 == 3 ? -3.0 : 0.0;
                cases.add(randomCase(String.format("smallN%dK%dv%d", n, k, variant), 2000 + variant, n, k, epochs,
                        maxStep, branchRates, simultaneous, zeroRates, ambiguous, logNeShift));
                variant++;
            }
        }
        cases.add(coalescentRateErrorCase());
        return cases;
    }

    static Case randomCase(String name, long seed, int tipCount, int stateCount, int epochCount, double maxStep,
                           boolean withBranchRates, boolean simultaneousSamples, boolean zeroRates,
                           boolean ambiguousTips, double logNeShift) {
        Random random = new Random(seed);
        int nodeCount = 2 * tipCount - 1;

        double[] sampleTimes = new double[tipCount];
        for (int i = 1; i < tipCount; i++) {
            double t = 4.0 * random.nextDouble();
            // simultaneous samples: a coarse grid of sampling times gives zero-length intervals
            sampleTimes[i] = simultaneousSamples ? Math.floor(t) : t;
        }
        Arrays.sort(sampleTimes);
        int[] tipOrder = shuffledRange(tipCount, random);

        int intervalCount = nodeCount - 1;
        double[] lengths = new double[intervalCount];
        int[] types = new int[intervalCount];
        int[] samples = new int[intervalCount];
        int[] child1 = new int[intervalCount];
        int[] child2 = new int[intervalCount];
        int[] parents = new int[intervalCount];
        double[] eventTimes = new double[intervalCount];

        List<Integer> active = new ArrayList<>();
        active.add(tipOrder[0]);
        double time = sampleTimes[0];
        int nextSample = 1;
        int nextParent = tipCount;
        for (int i = 0; i < intervalCount; i++) {
            int a = active.size();
            double rate = a >= 2 ? 0.25 * a * (a - 1) : 0.0;
            double wait = rate > 0.0 ? -Math.log(1.0 - random.nextDouble()) / rate : Double.POSITIVE_INFINITY;
            child1[i] = child2[i] = parents[i] = samples[i] = StructuredCoalescentSchedule.NO_NODE;
            if (nextSample < tipCount && time + wait >= sampleTimes[nextSample]) {
                lengths[i] = sampleTimes[nextSample] - time;
                types[i] = StructuredCoalescentSchedule.SAMPLE;
                samples[i] = tipOrder[nextSample];
                active.add(tipOrder[nextSample]);
                time = sampleTimes[nextSample];
                nextSample++;
            } else {
                lengths[i] = wait;
                types[i] = StructuredCoalescentSchedule.COALESCENT;
                int first = active.remove(random.nextInt(active.size()));
                int second = active.remove(random.nextInt(active.size()));
                child1[i] = first;
                child2[i] = second;
                parents[i] = nextParent;
                active.add(nextParent);
                nextParent++;
                time += wait;
            }
            eventTimes[i] = time;
        }
        StructuredCoalescentSchedule schedule = new StructuredCoalescentSchedule(sampleTimes[0], tipOrder[0],
                nodeCount, lengths, types, samples, child1, child2, parents);

        double[][] tips = new double[nodeCount][];
        for (int tip = 0; tip < tipCount; tip++) {
            tips[tip] = new double[stateCount];
            tips[tip][random.nextInt(stateCount)] = 1.0;
            if (ambiguousTips && random.nextDouble() < 0.2) {
                tips[tip][random.nextInt(stateCount)] = 1.0;
            }
        }
        MascotPreparedInput prepared = MascotPreparedInput.prepare(schedule,
                StructuredCoalescentTipData.fromPartials(nodeCount, stateCount, tips));

        double[] boundaries = epochBoundaries(epochCount, eventTimes, types, random);
        int migrationRates = stateCount * (stateCount - 1);
        int perEpoch = migrationRates + stateCount;
        double[] theta = new double[epochCount * perEpoch];
        for (int e = 0; e < epochCount; e++) {
            for (int i = 0; i < migrationRates; i++) {
                boolean off = zeroRates && random.nextDouble() < 0.3;
                theta[e * perEpoch + i] = off ? 0.0 : (0.5 + random.nextDouble()) * 4.0 / (stateCount - 1);
            }
            for (int s = 0; s < stateCount; s++) {
                theta[e * perEpoch + migrationRates + s] = 0.672 + (random.nextDouble() - 0.5) + logNeShift;
            }
        }

        double[] branchRates = null;
        if (withBranchRates) {
            branchRates = new double[nodeCount];
            for (int node = 0; node < nodeCount; node++) {
                branchRates[node] = 0.5 + random.nextDouble();
            }
        }
        return new Case(name, stateCount, maxStep, boundaries, theta, branchRates, prepared);
    }

    /** Two lineages sampled in different demes with no migration: lambda = 0 at their coalescence. */
    static Case coalescentRateErrorCase() {
        int stateCount = 2;
        StructuredCoalescentSchedule schedule = new StructuredCoalescentSchedule(0.0, 0, 3,
                new double[]{0.5, 1.0},
                new int[]{StructuredCoalescentSchedule.SAMPLE, StructuredCoalescentSchedule.COALESCENT},
                new int[]{1, StructuredCoalescentSchedule.NO_NODE},
                new int[]{StructuredCoalescentSchedule.NO_NODE, 0},
                new int[]{StructuredCoalescentSchedule.NO_NODE, 1},
                new int[]{StructuredCoalescentSchedule.NO_NODE, 2});
        double[][] tips = {{1.0, 0.0}, {0.0, 1.0}, null};
        MascotPreparedInput prepared = MascotPreparedInput.prepare(schedule,
                StructuredCoalescentTipData.fromPartials(3, stateCount, tips));
        return new Case("errorCoalescentRate", stateCount, 0.02, new double[]{0.0, Double.POSITIVE_INFINITY},
                new double[]{0.0, 0.0, 0.0, 0.0}, null, prepared);
    }

    /** One boundary exactly at a coalescent event time, a second between two events. */
    private static double[] epochBoundaries(int epochCount, double[] eventTimes, int[] types, Random random) {
        double[] boundaries = new double[epochCount + 1];
        boundaries[epochCount] = Double.POSITIVE_INFINITY;
        if (epochCount >= 2) {
            List<Double> coalescentTimes = new ArrayList<>();
            for (int i = 0; i < types.length; i++) {
                if (types[i] == StructuredCoalescentSchedule.COALESCENT && eventTimes[i] > 0.0) {
                    coalescentTimes.add(eventTimes[i]);
                }
            }
            double onEvent = coalescentTimes.get(coalescentTimes.size() / 2);
            boundaries[1] = onEvent;
            if (epochCount == 3) {
                // between events, or above the root when the boundary on an event is already at the root
                double span = Math.max(eventTimes[eventTimes.length - 1] - onEvent, 1.0);
                boundaries[2] = onEvent + span * (0.25 + 0.5 * random.nextDouble());
            }
        }
        return boundaries;
    }

    private static int[] shuffledRange(int n, Random random) {
        int[] order = new int[n];
        for (int i = 0; i < n; i++) {
            order[i] = i;
        }
        for (int i = n - 1; i > 0; i--) {
            int j = random.nextInt(i + 1);
            int t = order[i];
            order[i] = order[j];
            order[j] = t;
        }
        return order;
    }

    public static void write(Case c, File file) throws IOException {
        int stateCount = c.stateCount;
        int epochCount = c.boundaries.length - 1;
        StructuredCoalescentSchedule s = c.prepared.schedule;
        int lineageCount = c.prepared.maxLineageId + 1;

        double[] migration = new double[epochCount * stateCount * stateCount];
        double[] inversePopulation = new double[epochCount * stateCount];
        MascotRuntime.EpochRates rates = new MascotRuntime.EpochRates(stateCount);
        String status = "ok";
        try {
            for (int e = 0; e < epochCount; e++) {
                MascotRuntime.updateEpochRates(c.theta, e, stateCount, rates);
                System.arraycopy(rates.migrationMatrix, 0, migration, e * stateCount * stateCount,
                        stateCount * stateCount);
                System.arraycopy(rates.inversePopulation, 0, inversePopulation, e * stateCount, stateCount);
            }
        } catch (MascotLikelihoodDelegate.NumericalException ex) {
            status = "invalidRates";
        }

        MascotLikelihoodDelegate.Result result = null;
        if (status.equals("ok")) {
            try {
                result = new GenericMascotLikelihoodDelegate(stateCount, c.boundaries, c.maxStep)
                        .evaluate(c.prepared, c.theta, c.branchRates, true, true, false, false);
            } catch (MascotLikelihoodDelegate.NumericalException ex) {
                status = "numerical";
            }
        }

        try (PrintWriter w = new PrintWriter(file)) {
            w.println(FORMAT);
            w.println("name " + c.name);
            w.println("dims " + stateCount + " " + lineageCount + " " + s.getIntervalCount() + " " +
                    s.initialSampleNode + " " + epochCount);
            w.println("maxStep " + hex(c.maxStep));
            w.println("initialTime " + hex(s.initialTime));
            w.println("boundaries " + hex(c.boundaries));
            w.println("migration " + hex(migration));
            w.println("inversePopulation " + hex(inversePopulation));
            w.println("theta " + hex(c.theta));
            w.println("branchRates " + (c.branchRates == null ? "none" : hex(c.branchRates)));
            w.println("intervals");
            for (int i = 0; i < s.getIntervalCount(); i++) {
                w.println(hex(s.intervalLengths[i]) + " " + s.eventTypes[i] + " " + s.sampleNodes[i] + " " +
                        s.coalescentChild1[i] + " " + s.coalescentChild2[i] + " " + s.coalescentParent[i]);
            }
            w.println("tips");
            double[] partials = new double[stateCount];
            for (int node = 0; node < lineageCount; node++) {
                boolean has = c.prepared.tipData.hasTipPartials(node) || c.prepared.tipData.hasObservedState(node);
                if (has) {
                    c.prepared.tipData.writeTipPartials(node, partials, 0);
                    w.println(node + " 1 " + hex(partials));
                } else {
                    w.println(node + " 0");
                }
            }
            w.println("status " + status);
            if (result != null) {
                w.println("logL " + hex(result.logLikelihood));
                w.println("thetaGradient " + hex(result.gradient));
                w.println("clockGradient " + (result.clockGradient == null ? "none" : hex(result.clockGradient)));
                w.println("ancestralScores " + hex(result.ancestralStateScores));
            }
            w.println("end");
        }
    }

    static String hex(double x) {
        String bits = Long.toHexString(Double.doubleToRawLongBits(x));
        return "0000000000000000".substring(bits.length()) + bits;
    }

    static String hex(double[] x) {
        StringBuilder b = new StringBuilder();
        for (int i = 0; i < x.length; i++) {
            if (i > 0) {
                b.append(' ');
            }
            b.append(hex(x[i]));
        }
        return b.toString();
    }
}
