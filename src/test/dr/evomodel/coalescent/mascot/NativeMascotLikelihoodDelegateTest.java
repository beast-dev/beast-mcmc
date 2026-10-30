/*
 * NativeMascotLikelihoodDelegateTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent.mascot;

import beagle.mascot.MascotJNIWrapper;
import dr.evomodel.coalescent.StructuredCoalescentSchedule;
import dr.evomodel.coalescent.mascot.GenericMascotLikelihoodDelegate;
import dr.evomodel.coalescent.mascot.MascotGoldenCorpus;
import dr.evomodel.coalescent.mascot.MascotLikelihoodDelegate;
import dr.evomodel.coalescent.mascot.MascotPreparedInput;
import dr.evomodel.coalescent.mascot.NativeMascotLikelihoodDelegate;
import dr.evomodel.coalescent.mascot.StochasticMapSampler;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.atomic.AtomicBoolean;

/**
 * The native MASCOT delegate (beagle-lib's hmsbeagle-jni-mascot) against the Java delegate on the generated golden
 * cases: strict is bitwise on gradients and ancestral scores (logL up to an ulp of log(lambda) per coalescence), FMA is
 * within 1e-11 relative. Skipped, with a message, when the native library is not available (set mascot.library.path
 * to its build directory).
 */
public class NativeMascotLikelihoodDelegateTest extends TestCase {

    private static final double FMA_TOLERANCE = 1.0e-11;

    private static boolean nativeAvailable() {
        String reason = MascotJNIWrapper.getUnavailableReason();
        if (reason != null) {
            System.err.println("NativeMascotLikelihoodDelegateTest skipped: " + reason);
        }
        return reason == null;
    }

    private static final class Outputs {
        double logL;
        double[] gradient;
        double[] clock;
        double[] ancestral;
        boolean numerical;
    }

    private static Outputs evaluate(MascotLikelihoodDelegate delegate, MascotGoldenCorpus.Case c) {
        MascotPreparedInput prepared = c.getPreparedInput();
        int lineages = prepared.getSchedule().nodeCount;
        Outputs o = new Outputs();
        o.gradient = new double[c.theta.length];
        o.clock = c.branchRates == null ? null : new double[lineages];
        o.ancestral = new double[lineages * c.stateCount];
        try {
            o.logL = delegate.calculateLikelihoodAndDerivatives(prepared, c.theta, c.branchRates, o.gradient, o.clock,
                    o.ancestral, false);
            assertEquals(c.name + ": forward and forward+reverse logL", o.logL,
                    delegate.calculateLikelihood(prepared, c.theta, c.branchRates, false), 0.0);
        } catch (MascotLikelihoodDelegate.NumericalException e) {
            o.numerical = true;
        }
        return o;
    }

    private static NativeMascotLikelihoodDelegate nativeDelegate(MascotGoldenCorpus.Case c, boolean strict) {
        return NativeMascotLikelihoodDelegate.create(c.stateCount, c.getPreparedInput().getSchedule().nodeCount,
                c.boundaries, c.stepRule, strict);
    }

    private static int coalescences(MascotGoldenCorpus.Case c) {
        int count = 0;
        for (int type : c.getPreparedInput().getSchedule().eventTypes) {
            count += type == StructuredCoalescentSchedule.COALESCENT ? 1 : 0;
        }
        return count;
    }

    private static void assertBitwise(String label, double[] expected, double[] actual) {
        if (expected == null) {
            assertNull(label, actual);
            return;
        }
        assertEquals(label + " length", expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            boolean same = Double.isNaN(expected[i]) ? Double.isNaN(actual[i]) : expected[i] == actual[i];
            assertTrue(label + "[" + i + "]: " + expected[i] + " vs " + actual[i], same);
        }
    }

    /** Worst |a - e| / max(|e|, 1e-3 max |expected|), as the native golden tests measure it. */
    private static double relativeDifference(double[] expected, double[] actual) {
        if (expected == null) {
            return 0.0;
        }
        double scale = 0.0;
        for (double e : expected) {
            if (!Double.isNaN(e)) {
                scale = Math.max(scale, Math.abs(e));
            }
        }
        double worst = 0.0;
        for (int i = 0; i < expected.length; i++) {
            if (Double.isNaN(expected[i]) || Double.isNaN(actual[i])) {
                if (Double.isNaN(expected[i]) != Double.isNaN(actual[i])) {
                    return Double.POSITIVE_INFINITY;
                }
                continue;
            }
            worst = Math.max(worst, Math.abs(actual[i] - expected[i]) / Math.max(Math.abs(expected[i]), 1.0e-3 * scale));
        }
        return worst;
    }

    public void testStrictMatchesJava() {
        if (!nativeAvailable()) {
            return;
        }
        for (MascotGoldenCorpus.Case c : MascotGoldenCorpus.generatedCases()) {
            GenericMascotLikelihoodDelegate javaDelegate =
                    new GenericMascotLikelihoodDelegate(c.stateCount, c.boundaries, c.stepRule);
            NativeMascotLikelihoodDelegate strictDelegate = nativeDelegate(c, true);
            Outputs java = evaluate(javaDelegate, c);
            Outputs strict = evaluate(strictDelegate, c);
            assertEquals(c.name + ": numerical failure", java.numerical, strict.numerical);
            if (java.numerical) {
                continue;
            }
            double ulp = Math.ulp(java.logL);
            assertEquals(c.name + ": logL", java.logL, strict.logL, 2.0 * ulp + coalescences(c) * 3.6e-15);
            assertBitwise(c.name + ": theta gradient", java.gradient, strict.gradient);
            assertBitwise(c.name + ": clock gradient", java.clock, strict.clock);
            assertBitwise(c.name + ": ancestral scores", java.ancestral, strict.ancestral);
            assertTrue(c.name + ": step counts", Arrays.equals(javaDelegate.getLastStepCounts(),
                    strictDelegate.getLastStepCounts()));
        }
    }

    public void testFmaWithinTolerance() {
        if (!nativeAvailable()) {
            return;
        }
        for (MascotGoldenCorpus.Case c : MascotGoldenCorpus.generatedCases()) {
            Outputs java = evaluate(new GenericMascotLikelihoodDelegate(c.stateCount, c.boundaries, c.stepRule), c);
            NativeMascotLikelihoodDelegate fmaDelegate = nativeDelegate(c, false);
            Outputs fma = evaluate(fmaDelegate, c);
            assertEquals(c.name + ": numerical failure", java.numerical, fma.numerical);
            if (java.numerical) {
                continue;
            }
            double worst = Math.abs(fma.logL - java.logL) / Math.abs(java.logL);
            worst = Math.max(worst, relativeDifference(java.gradient, fma.gradient));
            worst = Math.max(worst, relativeDifference(java.clock, fma.clock));
            worst = Math.max(worst, relativeDifference(java.ancestral, fma.ancestral));
            assertTrue(c.name + ": FMA differs from Java by " + worst, worst <= FMA_TOLERANCE);
        }
    }

    /** Two delegates with their own instances, evaluated concurrently, reproduce their serial results. */
    public void testConcurrentInstances() throws InterruptedException {
        if (!nativeAvailable()) {
            return;
        }
        List<MascotGoldenCorpus.Case> cases = new ArrayList<>();
        for (MascotGoldenCorpus.Case c : MascotGoldenCorpus.generatedCases()) {
            if (c.name.equals("sweepK5") || c.name.equals("sweepK8")) {
                cases.add(c);
            }
        }
        assertEquals(2, cases.size());
        final AtomicBoolean same = new AtomicBoolean(true);
        List<Thread> threads = new ArrayList<>();
        for (final MascotGoldenCorpus.Case c : cases) {
            final NativeMascotLikelihoodDelegate delegate = nativeDelegate(c, false);
            final Outputs serial = evaluate(delegate, c);
            threads.add(new Thread(() -> {
                for (int i = 0; i < 20; i++) {
                    Outputs o = evaluate(delegate, c);
                    if (o.logL != serial.logL || relativeDifference(serial.gradient, o.gradient) != 0.0) {
                        same.set(false);
                    }
                }
            }));
        }
        for (Thread t : threads) {
            t.start();
        }
        for (Thread t : threads) {
            t.join();
        }
        assertTrue("concurrent native instances reproduce their serial results", same.get());
    }

    /**
     * Single-interval changes accepted or rejected at random, as an MCMC would: an incremental delegate reproduces a
     * full pass bitwise and resumes at the changed interval; a rate change resumes from the start.
     */
    public void testIncrementalRestartMatchesFullPass() {
        if (!nativeAvailable()) {
            return;
        }
        if (!MascotJNIWrapper.supportsIncrementalRestart()) {
            System.err.println("testIncrementalRestartMatchesFullPass skipped: the library predates API 1.2");
            return;
        }
        Random random = new Random(7);
        int tested = 0;
        for (MascotGoldenCorpus.Case c : MascotGoldenCorpus.generatedCases()) {
            if (!c.name.equals("sweepK5") && !c.name.equals("sweepK8") && !c.name.equals("adaptivesweepK5") &&
                    !c.name.startsWith("smallN80K7")) {
                continue;
            }
            tested++;
            MascotPreparedInput prepared = c.getPreparedInput();
            double[] lengths = prepared.getSchedule().intervalLengths;
            double[] original = lengths.clone();
            try {
                for (boolean strict : new boolean[]{true, false}) {
                    String label = c.name + (strict ? " strict" : " fma");
                    NativeMascotLikelihoodDelegate full = NativeMascotLikelihoodDelegate.create(c.stateCount,
                            prepared.getSchedule().nodeCount, c.boundaries, c.stepRule, strict, false);
                    NativeMascotLikelihoodDelegate incremental = NativeMascotLikelihoodDelegate.create(c.stateCount,
                            prepared.getSchedule().nodeCount, c.boundaries, c.stepRule, strict, true);
                    double stored = incremental.calculateLikelihood(prepared, c.theta, c.branchRates, false);
                    assertEquals(label + ": first pass", full.calculateLikelihood(prepared, c.theta, c.branchRates,
                            false), stored, 0.0);
                    assertEquals(label + ": first pass restart", 0, incremental.getLastRestartInterval());
                    incremental.storeCheckpoints();
                    for (int step = 0; step < 40; step++) {
                        int j;
                        do {
                            j = random.nextInt(lengths.length);
                        } while (lengths[j] == 0.0);
                        double old = lengths[j];
                        lengths[j] = old * (0.8 + 0.4 * random.nextDouble());
                        double expected = full.calculateLikelihood(prepared, c.theta, c.branchRates, false);
                        double actual = incremental.calculateLikelihood(prepared, c.theta, c.branchRates, false);
                        assertEquals(label + " step " + step + ": logL", expected, actual, 0.0);
                        assertEquals(label + " step " + step + ": restart", j, incremental.getLastRestartInterval());
                        if (random.nextBoolean()) {
                            incremental.storeCheckpoints();
                            stored = actual;
                        } else {
                            lengths[j] = old;
                            incremental.restoreCheckpoints();
                            if (step % 4 == 0) {
                                assertEquals(label + " step " + step + ": restored logL", stored,
                                        incremental.calculateLikelihood(prepared, c.theta, c.branchRates, false), 0.0);
                                assertEquals(label + " step " + step + ": restored restart", lengths.length,
                                        incremental.getLastRestartInterval());
                            }
                        }
                    }
                    double[] theta = c.theta.clone();
                    theta[theta.length - 1] += 0.01;
                    assertEquals(label + ": rate change", full.calculateLikelihood(prepared, theta, c.branchRates,
                            false), incremental.calculateLikelihood(prepared, theta, c.branchRates, false), 0.0);
                    assertEquals(label + ": rate change restart", 0, incremental.getLastRestartInterval());
                }
            } finally {
                System.arraycopy(original, 0, lengths, 0, lengths.length);
            }
        }
        assertEquals("cases tested", 4, tested);
    }

    /** Stochastic mapping under the native delegate is the Java mapper at mappingMaxStep (here the same step). */
    public void testStochasticMappingRunsInJava() {
        if (!nativeAvailable()) {
            return;
        }
        for (MascotGoldenCorpus.Case c : MascotGoldenCorpus.generatedCases()) {
            if (!c.name.startsWith("smallN20")) {
                continue;
            }
            double[] heights = nodeHeights(c.getPreparedInput().getSchedule());
            GenericMascotLikelihoodDelegate java = new GenericMascotLikelihoodDelegate(c.stateCount, c.boundaries,
                    c.maxStep);
            NativeMascotLikelihoodDelegate nativeDelegate = nativeDelegate(c, false);
            double expected = java.prepareStochasticMap(c.getPreparedInput(), c.theta, c.branchRates, heights, false,
                    StochasticMapSampler.GRID);
            double actual = nativeDelegate.prepareStochasticMap(c.getPreparedInput(), c.theta, c.branchRates, heights,
                    false, StochasticMapSampler.GRID);
            assertEquals(c.name + ": mapping forward pass", expected, actual, 0.0);
            assertNotNull(c.name + ": a map is drawn", nativeDelegate.drawStochasticMap());
        }
    }

    private static double[] nodeHeights(StructuredCoalescentSchedule s) {
        double[] heights = new double[s.nodeCount];
        double time = s.initialTime;
        heights[s.initialSampleNode] = time;
        for (int i = 0; i < s.getIntervalCount(); i++) {
            time += s.intervalLengths[i];
            int node = s.eventTypes[i] == StructuredCoalescentSchedule.SAMPLE ? s.sampleNodes[i] : s.coalescentParent[i];
            heights[node] = time;
        }
        return heights;
    }
}
