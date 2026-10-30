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
 * The step budget (ADAPTIVE with maxStep &gt; 0 and stepBudget &gt; 0) bounds a whole pass: over its segments, in order,
 * the evaluator sums the step counts and the reference counts of a FIXED maxStep grid ({@link #referenceSteps}), and
 * fails with {@link MascotLikelihoodDelegate.StepBudgetException} before integrating a segment that takes the steps
 * above stepBudget times the reference. So a pass never costs more than stepBudget times the FIXED maxStep pass, in time
 * or tape memory, however stiff the parameters an HMC trajectory proposes; the likelihood there is -infinity.
 * The step count is piecewise constant in the parameters, so gradients are those of the discretized likelihood.
 *
 * @author Marc A. Suchard
 */
public final class MascotStepRule {

    public enum Kind { FIXED, ADAPTIVE }

    public static final double DEFAULT_TAU = 0.25;
    /** A pass may take at most this many times the steps of the FIXED maxStep grid (ADAPTIVE). */
    public static final int DEFAULT_STEP_BUDGET = 100;

    private static final double LARGEST_STEP_COUNT = 1 << 30;

    final Kind kind;
    final double maxStep;
    final double tau;
    final double minStep;
    final int maxStepsPerSegment;
    final int stepBudget;

    private MascotStepRule(Kind kind, double maxStep, double tau, double minStep, int maxStepsPerSegment,
                           int stepBudget) {
        this.kind = kind;
        this.maxStep = maxStep;
        this.tau = tau;
        this.minStep = minStep;
        this.maxStepsPerSegment = maxStepsPerSegment;
        this.stepBudget = stepBudget;
    }

    public static MascotStepRule fixed(double maxStep) {
        if (!(maxStep > 0.0) || Double.isInfinite(maxStep)) {
            throw new IllegalArgumentException("maxStep must be positive and finite");
        }
        return new MascotStepRule(Kind.FIXED, maxStep, 0.0, 0.0, 0, 0);
    }

    /** {@link #adaptive(double, double, double, int, int)} with the default step budget. */
    public static MascotStepRule adaptive(double maxStep, double tau, double minStep, int maxStepsPerSegment) {
        return adaptive(maxStep, tau, minStep, maxStepsPerSegment, DEFAULT_STEP_BUDGET);
    }

    /**
     * @param maxStep            upper bound on the step, or &lt;= 0 for none
     * @param tau                target of length * lambda per step
     * @param minStep            lower bound on the step, or &lt;= 0 for none
     * @param maxStepsPerSegment upper bound on the step count, or &lt;= 0 for none
     * @param stepBudget         a pass may take at most this many times the steps of the FIXED maxStep grid, or
     *                           &lt;= 0 for no bound; it applies only with maxStep &gt; 0
     */
    public static MascotStepRule adaptive(double maxStep, double tau, double minStep, int maxStepsPerSegment,
                                          int stepBudget) {
        if (!(tau > 0.0) || Double.isInfinite(tau)) {
            throw new IllegalArgumentException("tau must be positive and finite");
        }
        if (Double.isNaN(maxStep) || Double.isNaN(minStep)) {
            throw new IllegalArgumentException("maxStep and minStep must not be NaN");
        }
        return new MascotStepRule(Kind.ADAPTIVE, maxStep, tau, minStep, maxStepsPerSegment, stepBudget);
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

    public int getStepBudget() {
        return stepBudget;
    }

    /** Whether passes are held to the step budget: ADAPTIVE with maxStep &gt; 0 and stepBudget &gt; 0. */
    public boolean hasStepBudget() {
        return kind == Kind.ADAPTIVE && maxStep > 0.0 && stepBudget > 0;
    }

    /** The step count of the FIXED maxStep grid for a segment, the budget's reference (maxStep &gt; 0). */
    int referenceSteps(double length) {
        return Math.max(1, (int) Math.ceil(Math.min(length / maxStep, LARGEST_STEP_COUNT)));
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
