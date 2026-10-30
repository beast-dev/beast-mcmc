/*
 * MascotStepRuleTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent.mascot;

import dr.evomodel.coalescent.StructuredCoalescentSchedule;
import dr.evomodel.coalescent.mascot.GenericMascotLikelihoodDelegate;
import dr.evomodel.coalescent.mascot.MascotGoldenCorpus;
import dr.evomodel.coalescent.mascot.MascotStepRule;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** The fixed and adaptive RK4 step rules (R3) of the Java MASCOT delegate, on the generated golden cases. */
public class MascotStepRuleTest extends TestCase {

    private static List<MascotGoldenCorpus.Case> cases(String... names) {
        List<MascotGoldenCorpus.Case> selected = new ArrayList<>();
        for (MascotGoldenCorpus.Case c : MascotGoldenCorpus.generatedCases()) {
            if (names.length == 0 || Arrays.asList(names).contains(c.name)) {
                selected.add(c);
            }
        }
        return selected;
    }

    private static double logL(MascotGoldenCorpus.Case c, MascotStepRule rule, double[] theta, double[] gradient) {
        GenericMascotLikelihoodDelegate d = new GenericMascotLikelihoodDelegate(c.stateCount, c.boundaries, rule);
        return gradient == null
                ? d.calculateLikelihood(c.getPreparedInput(), theta, c.branchRates, false)
                : d.calculateLikelihoodAndDerivatives(c.getPreparedInput(), theta, c.branchRates, gradient, null, null,
                false);
    }

    private static int[] stepCounts(MascotGoldenCorpus.Case c, MascotStepRule rule, double[] theta) {
        GenericMascotLikelihoodDelegate d = new GenericMascotLikelihoodDelegate(c.stateCount, c.boundaries, rule);
        d.calculateLikelihood(c.getPreparedInput(), theta, c.branchRates, false);
        return d.getLastStepCounts();
    }

    private static final double TIME_TOLERANCE = 1.0e-14;

    /** Lengths of the integration segments: intervals split at epoch boundaries, as MascotForwardEvaluator does. */
    private static double[] segmentLengths(MascotGoldenCorpus.Case c) {
        StructuredCoalescentSchedule s = c.getPreparedInput().getSchedule();
        List<Double> lengths = new ArrayList<>();
        double time = s.initialTime;
        for (double length : s.intervalLengths) {
            double end = time + length;
            while (end > time + TIME_TOLERANCE) {
                double segmentEnd = end;
                for (int b = 1; b < c.boundaries.length; b++) {
                    if (c.boundaries[b] > time + TIME_TOLERANCE) {
                        if (c.boundaries[b] <= end + TIME_TOLERANCE) {
                            segmentEnd = Math.min(c.boundaries[b], end);
                        }
                        break;
                    }
                }
                lengths.add(segmentEnd - time);
                time = segmentEnd;
            }
            time = end;
        }
        double[] out = new double[lengths.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = lengths.get(i);
        }
        return out;
    }

    public void testFixedRuleIsTheMaxStepRule() {
        for (MascotGoldenCorpus.Case c : cases()) {
            if (c.name.startsWith("error") || c.stepRule.getKind() != MascotStepRule.Kind.FIXED) {
                continue;
            }
            assertEquals(c.name, logL(c, MascotStepRule.fixed(c.maxStep), c.theta, null),
                    new GenericMascotLikelihoodDelegate(c.stateCount, c.boundaries, c.maxStep)
                            .calculateLikelihood(c.getPreparedInput(), c.theta, c.branchRates, false), 0.0);
            double[] lengths = segmentLengths(c);
            int[] counts = stepCounts(c, MascotStepRule.fixed(c.maxStep), c.theta);
            assertEquals(c.name + ": segments", lengths.length, counts.length);
            for (int i = 0; i < counts.length; i++) {
                assertEquals(c.name + ": segment " + i, Math.max(1, (int) Math.ceil(lengths[i] / c.maxStep)), counts[i]);
            }
        }
    }

    public void testAdaptiveClamps() {
        for (MascotGoldenCorpus.Case c : cases("sweepK5", "smallN40K5v22", "smallN80K7v27")) {
            double[] lengths = segmentLengths(c);
            int[] free = stepCounts(c, MascotStepRule.adaptive(0.0, 0.25, 0.0, 0), c.theta);
            int[] capped = stepCounts(c, MascotStepRule.adaptive(0.02, 0.25, 0.0, 0), c.theta);
            int[] floored = stepCounts(c, MascotStepRule.adaptive(0.0, 0.25, 0.05, 0), c.theta);
            int[] limited = stepCounts(c, MascotStepRule.adaptive(0.0, 0.25, 0.0, 3), c.theta);
            boolean anyAdaptive = false;
            for (int i = 0; i < lengths.length; i++) {
                assertTrue(free[i] >= 1);
                assertEquals(c.name + ": cap", Math.max(free[i], (int) Math.ceil(lengths[i] / 0.02)), capped[i]);
                assertEquals(c.name + ": floor", Math.min(free[i], Math.max(1, (int) Math.floor(lengths[i] / 0.05))),
                        floored[i]);
                assertEquals(c.name + ": limit", Math.min(free[i], 3), limited[i]);
                anyAdaptive |= free[i] > 1;
            }
            assertTrue(c.name + ": the rule never asked for more than one step", anyAdaptive);
        }
    }

    /** A smaller tau gives more steps and a logL closer to a fine fixed-step reference. */
    public void testAdaptiveConvergesToFineSteps() {
        for (MascotGoldenCorpus.Case c : cases("sweepK5", "smallN40K5v22")) {
            double reference = logL(c, MascotStepRule.fixed(1.0e-4), c.theta, null);
            double previous = Double.POSITIVE_INFINITY;
            int previousSteps = 0;
            for (double tau : new double[]{1.0, 0.25, 0.05, 0.01}) {
                MascotStepRule rule = MascotStepRule.adaptive(0.0, tau, 0.0, 0);
                double error = Math.abs(logL(c, rule, c.theta, null) - reference);
                int steps = 0;
                for (int n : stepCounts(c, rule, c.theta)) {
                    steps += n;
                }
                assertTrue(c.name + ": tau " + tau + " error " + error + " not below " + previous, error < previous);
                assertTrue(c.name + ": tau " + tau + " uses no more steps", steps > previousSteps);
                previous = error;
                previousSteps = steps;
            }
            assertTrue(c.name + ": tau 0.01 error " + previous, previous < 1.0e-8);
        }
    }

    /**
     * Gradients are those of the discretized likelihood: central differences match wherever the perturbation leaves
     * every step count unchanged (logL jumps where a count changes).
     */
    public void testAdaptiveGradientsMatchFiniteDifferences() {
        for (MascotGoldenCorpus.Case c : cases("smallN20K3v17", "smallN40K5v22", "sweepK5")) {
            MascotStepRule rule = MascotStepRule.adaptive(0.1, 0.25, 0.0, 0);
            double[] gradient = new double[c.theta.length];
            logL(c, rule, c.theta, gradient);
            int[] counts = stepCounts(c, rule, c.theta);
            double scale = 0.0;
            for (double g : gradient) {
                scale = Math.max(scale, Math.abs(g));
            }
            int checked = 0;
            for (int i = 0; i < c.theta.length; i++) {
                if (c.theta[i] == 0.0) {
                    continue;
                }
                double h = 1.0e-6 * Math.max(Math.abs(c.theta[i]), 1.0);
                double[] up = c.theta.clone();
                double[] down = c.theta.clone();
                up[i] += h;
                down[i] -= h;
                if (!Arrays.equals(counts, stepCounts(c, rule, up)) || !Arrays.equals(counts, stepCounts(c, rule, down))) {
                    continue;
                }
                double numeric = (logL(c, rule, up, null) - logL(c, rule, down, null)) / (2.0 * h);
                double error = Math.abs(numeric - gradient[i]) / Math.max(Math.abs(gradient[i]), 1.0e-2 * scale);
                assertTrue(c.name + ": theta[" + i + "] analytic " + gradient[i] + " numeric " + numeric,
                        error < 1.0e-5);
                checked++;
            }
            assertTrue(c.name + ": only " + checked + " of " + c.theta.length + " entries checked",
                    checked >= c.theta.length / 2);
        }
    }
}
