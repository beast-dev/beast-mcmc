/*
 * MascotStepRule.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

/**
 * The number of equal RK4 steps for each integration segment (an interval between tree events, split at epoch
 * boundaries). The native kernels (beagle-lib mascot.h, MascotStepRule) follow the same rule with the same operation
 * order, so all implementations take identical step counts.
 * <ul>
 * <li>FIXED: n = max(1, ceil(length / maxStep)).</li>
 * <li>ADAPTIVE: from the stiffness lambda of the state at the segment's start (see {@link #stiffness}),
 * n = max(1, ceil(length * lambda / tau)); then at least ceil(length / maxStep) if maxStep &gt; 0, so h never exceeds
 * maxStep; then at most max(1, floor(length / minStep)) if minStep &gt; 0, and at most maxStepsPerSegment if &gt; 0. The
 * last two guard against a pathological lambda and win over maxStep.</li>
 * </ul>
 * The step count is piecewise constant in the parameters, so gradients are those of the discretized likelihood.
 *
 * @author Marc A. Suchard
 */
public final class MascotStepRule {

    public enum Kind { FIXED, ADAPTIVE }

    public static final double DEFAULT_TAU = 0.25;

    private static final double LARGEST_STEP_COUNT = 1 << 30;

    final Kind kind;
    final double maxStep;
    final double tau;
    final double minStep;
    final int maxStepsPerSegment;

    private MascotStepRule(Kind kind, double maxStep, double tau, double minStep, int maxStepsPerSegment) {
        this.kind = kind;
        this.maxStep = maxStep;
        this.tau = tau;
        this.minStep = minStep;
        this.maxStepsPerSegment = maxStepsPerSegment;
    }

    public static MascotStepRule fixed(double maxStep) {
        if (!(maxStep > 0.0) || Double.isInfinite(maxStep)) {
            throw new IllegalArgumentException("maxStep must be positive and finite");
        }
        return new MascotStepRule(Kind.FIXED, maxStep, 0.0, 0.0, 0);
    }

    /**
     * @param maxStep            upper bound on the step, or &lt;= 0 for none
     * @param tau                target of length * lambda per step
     * @param minStep            lower bound on the step, or &lt;= 0 for none
     * @param maxStepsPerSegment upper bound on the step count, or &lt;= 0 for none
     */
    public static MascotStepRule adaptive(double maxStep, double tau, double minStep, int maxStepsPerSegment) {
        if (!(tau > 0.0) || Double.isInfinite(tau)) {
            throw new IllegalArgumentException("tau must be positive and finite");
        }
        if (Double.isNaN(maxStep) || Double.isNaN(minStep)) {
            throw new IllegalArgumentException("maxStep and minStep must not be NaN");
        }
        return new MascotStepRule(Kind.ADAPTIVE, maxStep, tau, minStep, maxStepsPerSegment);
    }

    public Kind getKind() {
        return kind;
    }

    public double getMaxStep() {
        return maxStep;
    }

    public double getTau() {
        return tau;
    }

    public double getMinStep() {
        return minStep;
    }

    public int getMaxStepsPerSegment() {
        return maxStepsPerSegment;
    }

    /**
     * lambda = clock * maxOut + maxLoad, in this order:
     * clock = 1 without branch rates, else the largest active-lineage rate (from 0.0, lineages in slot order);
     * maxOut = largest -M_kk (from 0.0, k = 0..K-1); maxLoad = largest q_k * load_k (from 0.0, k = 0..K-1), where
     * load_k sums p_lk over active lineages in slot order from 0.0.
     */
    static double stiffness(double[] probabilities, int activeCount, int stateCount, double[] migrationMatrix,
                            double[] inversePopulation, double[] activeClockRates) {
        double clock = 1.0;
        if (activeClockRates != null) {
            clock = 0.0;
            for (int i = 0; i < activeCount; i++) {
                clock = Math.max(clock, activeClockRates[i]);
            }
        }
        double maxOut = 0.0;
        for (int k = 0; k < stateCount; k++) {
            maxOut = Math.max(maxOut, -migrationMatrix[k * stateCount + k]);
        }
        double maxLoad = 0.0;
        for (int k = 0; k < stateCount; k++) {
            double load = 0.0;
            for (int i = 0; i < activeCount; i++) {
                load += probabilities[i * stateCount + k];
            }
            maxLoad = Math.max(maxLoad, inversePopulation[k] * load);
        }
        return clock * maxOut + maxLoad;
    }

    /** Step count for a segment of the given length; lambda is only used by ADAPTIVE. */
    int steps(double length, double lambda) {
        if (kind == Kind.FIXED) {
            return Math.max(1, (int) Math.ceil(length / maxStep));
        }
        double target = length * lambda / tau;
        int steps = target > 0.0 ? (int) Math.ceil(Math.min(target, LARGEST_STEP_COUNT)) : 1;
        steps = Math.max(1, steps);
        if (maxStep > 0.0) {
            steps = Math.max(steps, (int) Math.ceil(Math.min(length / maxStep, LARGEST_STEP_COUNT)));
        }
        if (minStep > 0.0) {
            steps = Math.min(steps, Math.max(1, (int) Math.floor(Math.min(length / minStep, LARGEST_STEP_COUNT))));
        }
        if (maxStepsPerSegment > 0) {
            steps = Math.min(steps, maxStepsPerSegment);
        }
        return steps;
    }
}
