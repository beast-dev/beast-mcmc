/*
 * StepProfile.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.markovchain;

import dr.inference.operators.MCMCOperator;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.function.Supplier;

/**
 * Wall-clock time of the phases of every MCMC step, per operator: storing the model, the operator's
 * proposal, the evaluation of the posterior, accepting or restoring, and adaptation. Enabled with
 * -Dmcmc.profile=true (otherwise the timing calls are not made); printed after the operator analysis.
 * Only the chain after the adaptation delay is reported (the sums are reset when a chain starts, as the
 * operators' counts are).
 *
 * Other profiles can register a report to be printed with this one ({@link #addReporter}).
 *
 * @author Andrew Rambaut
 */
public final class StepProfile {

    public static final boolean ENABLED = Boolean.getBoolean("mcmc.profile");

    public enum Phase {
        STORE("store"),
        PROPOSE("propose"),
        EVALUATE("evaluate"),
        ACCEPT_RESTORE("accept/restore"),
        ADAPT("adapt");

        final String label;

        Phase(String label) {
            this.label = label;
        }
    }

    private static final List<Supplier<String>> reporters = new CopyOnWriteArrayList<Supplier<String>>();
    private static final List<Runnable> resetHooks = new CopyOnWriteArrayList<Runnable>();
    private static volatile StepProfile active = null;
    private static volatile MCMCOperator currentOperator = null;

    /**
     * Registers an action to run whenever a chain's step profile is reset (so that other profiles cover the
     * same steps).
     */
    public static void addResetHook(Runnable hook) {
        resetHooks.add(hook);
    }

    /**
     * @return the name of the operator whose step is in progress, or null outside a step
     */
    public static String getCurrentOperatorName() {
        final MCMCOperator operator = currentOperator;
        return operator == null ? null : operator.getOperatorName();
    }

    /**
     * @return the number of steps of the named operator profiled so far by the active chain (0 if none)
     */
    public static long getSteps(String operatorName) {
        final StepProfile profile = active;
        if (profile != null) {
            for (Map.Entry<MCMCOperator, long[]> entry : profile.times.entrySet()) {
                if (entry.getKey().getOperatorName().equals(operatorName)) {
                    return entry.getValue()[Phase.values().length];
                }
            }
        }
        return 0;
    }

    /**
     * Registers a report to be printed after this profile's (e.g. a likelihood's own phase timings).
     */
    public static void addReporter(Supplier<String> reporter) {
        reporters.add(reporter);
    }

    public static List<String> getReports() {
        List<String> reports = new ArrayList<String>();
        for (Supplier<String> reporter : reporters) {
            reports.add(reporter.get());
        }
        return reports;
    }

    // per operator: the nanoseconds in each phase, then the number of steps
    private final Map<MCMCOperator, long[]> times = new LinkedHashMap<MCMCOperator, long[]>();

    public StepProfile() {
        active = this;
    }

    /**
     * Marks the start of a step of the operator (so that other profiles can attribute their phases to it).
     */
    public void startStep(MCMCOperator operator) {
        currentOperator = operator;
    }

    /**
     * Adds the time since start to the operator's phase (a step is counted at its STORE phase).
     *
     * @return the current time, the start of the next phase
     */
    public long add(MCMCOperator operator, Phase phase, long start) {
        final long now = System.nanoTime();
        long[] t = times.get(operator);
        if (t == null) {
            t = new long[Phase.values().length + 1];
            times.put(operator, t);
        }
        t[phase.ordinal()] += now - start;
        if (phase == Phase.STORE) {
            t[Phase.values().length]++;
        }
        return now;
    }

    public void reset() {
        times.clear();
        for (Runnable hook : resetHooks) {
            hook.run();
        }
    }

    public String report() {
        final Phase[] phases = Phase.values();
        final StringBuilder sb = new StringBuilder();
        sb.append("\nStep profile (microseconds per step, wall clock, after the adaptation delay)\n");
        sb.append(String.format("%-50s %10s", "Operator", "Steps"));
        for (Phase phase : phases) {
            sb.append(String.format(" %14s", phase.label));
        }
        sb.append(String.format(" %14s%n", "total"));
        final long[] all = new long[phases.length + 1];
        for (Map.Entry<MCMCOperator, long[]> entry : times.entrySet()) {
            appendRow(sb, entry.getKey().getOperatorName(), entry.getValue());
            for (int i = 0; i < all.length; i++) {
                all[i] += entry.getValue()[i];
            }
        }
        appendRow(sb, "(all operators)", all);
        return sb.toString();
    }

    private static void appendRow(StringBuilder sb, String name, long[] t) {
        final int steps = (int) t[Phase.values().length];
        sb.append(String.format("%-50s %10d", name, steps));
        long total = 0;
        for (int i = 0; i < Phase.values().length; i++) {
            total += t[i];
            sb.append(String.format(" %14.2f", steps == 0 ? 0.0 : t[i] / 1.0e3 / steps));
        }
        sb.append(String.format(" %14.2f%n", steps == 0 ? 0.0 : total / 1.0e3 / steps));
    }
}
