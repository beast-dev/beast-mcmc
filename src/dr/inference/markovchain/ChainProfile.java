/*
 * ChainProfile.java
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

import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Profileable;
import dr.inference.operators.MCMCOperator;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.function.Supplier;

/**
 * Wall-clock profile of a Markov chain, printed after the operator analysis. Which aspects are profiled is
 * chosen with the system property {@code mcmc.profile}, a comma-separated list of:
 * <ul>
 * <li>{@code steps}: the phases of every step per operator: storing the model, the operator's proposal,
 * the evaluation of the posterior, accepting or restoring, and adaptation;</li>
 * <li>{@code likelihoods}: every evaluation of every likelihood component (the components of the compound
 * likelihoods and the chain's posterior itself), with the number of calculations reported by
 * {@link Profileable} likelihoods, in total and per operator;</li>
 * <li>any further name a component asks for with {@link #isEnabled} (e.g. {@code mat} for the MAT
 * likelihood's own phases), which registers its report with {@link #addReporter};</li>
 * <li>{@code all} (or {@code true}): everything.</li>
 * </ul>
 * When an aspect is off its timing calls are not made (the flags are compile-time constants). Only the
 * chain after the adaptation delay is reported: the sums are reset when a chain starts, as the operators'
 * counts are.
 *
 * @author Andrew Rambaut
 */
public final class ChainProfile {

    public static final String PROPERTY = "mcmc.profile";

    private static final Set<String> aspects = parseAspects(System.getProperty(PROPERTY));

    /**
     * @return whether the named aspect is to be profiled
     */
    public static boolean isEnabled(String aspect) {
        return aspects.contains("all") || aspects.contains(aspect.toLowerCase());
    }

    public static final boolean STEPS = isEnabled("steps");
    public static final boolean LIKELIHOODS = isEnabled("likelihoods");
    public static final boolean ENABLED = !aspects.isEmpty();

    private static Set<String> parseAspects(String value) {
        final Set<String> set = new HashSet<String>();
        if (value == null) {
            return set;
        }
        final String v = value.trim().toLowerCase();
        if (v.isEmpty() || v.equals("false") || v.equals("none")) {
            return set;
        }
        if (v.equals("true")) {
            set.add("all");
            return set;
        }
        for (String aspect : v.split(",")) {
            if (!aspect.trim().isEmpty()) {
                set.add(aspect.trim());
            }
        }
        return set;
    }

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
    private static volatile ChainProfile active = null;
    private static volatile MCMCOperator currentOperator = null;

    /**
     * Registers an action to run whenever a chain's profile is reset (so that other profiles cover the
     * same steps).
     */
    public static void addResetHook(Runnable hook) {
        resetHooks.add(hook);
    }

    /**
     * Registers a report to be printed after the chain's own (e.g. a likelihood's internal phases).
     */
    public static void addReporter(Supplier<String> reporter) {
        reporters.add(reporter);
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
        final ChainProfile profile = active;
        if (profile != null) {
            for (Map.Entry<MCMCOperator, long[]> entry : profile.stepCounts.entrySet()) {
                if (entry.getKey().getOperatorName().equals(operatorName)) {
                    return entry.getValue()[0];
                }
            }
        }
        return 0;
    }

    /**
     * Evaluates the likelihood, timing it for the active chain's profile. Callers use it in place of
     * {@code likelihood.getLogLikelihood()} when {@link #LIKELIHOODS} is on; it may be called from any
     * thread.
     */
    public static double evaluate(Likelihood likelihood) {
        final long calculationsBefore = getCalculationCount(likelihood);
        final long start = System.nanoTime();
        final double logLikelihood = likelihood.getLogLikelihood();
        final long elapsed = System.nanoTime() - start;
        final ChainProfile profile = active;
        if (profile != null) {
            profile.addEvaluation(likelihood, elapsed, getCalculationCount(likelihood) - calculationsBefore);
        }
        return logLikelihood;
    }

    private static long getCalculationCount(Likelihood likelihood) {
        return likelihood instanceof Profileable ? ((Profileable) likelihood).getTotalCalculationCount() : 0L;
    }

    // the chain's posterior, the root of the likelihood tree in the report
    private final Likelihood root;

    // per operator: the number of steps (whatever the aspects)
    private final Map<MCMCOperator, long[]> stepCounts = new LinkedHashMap<MCMCOperator, long[]>();

    // per operator: the nanoseconds in each phase (the steps aspect)
    private final Map<MCMCOperator, long[]> stepTimes = new LinkedHashMap<MCMCOperator, long[]>();

    // per likelihood: nanoseconds, evaluations, calculations, updated under the array's lock (the
    // likelihoods aspect)
    private final Map<Likelihood, long[]> likelihoodTimes = new ConcurrentHashMap<Likelihood, long[]>();

    // the same per operator
    private final Map<MCMCOperator, Map<Likelihood, long[]>> likelihoodTimesByOperator =
            new ConcurrentHashMap<MCMCOperator, Map<Likelihood, long[]>>();

    public ChainProfile(Likelihood root) {
        this.root = root;
        active = this;
    }

    /**
     * Marks the start of a step of the operator: counts it and lets evaluations and other profiles be
     * attributed to it.
     */
    public void startStep(MCMCOperator operator) {
        currentOperator = operator;
        long[] count = stepCounts.get(operator);
        if (count == null) {
            count = new long[1];
            stepCounts.put(operator, count);
        }
        count[0]++;
    }

    /**
     * Adds the time since start to the operator's phase.
     *
     * @return the current time, the start of the next phase
     */
    public long add(MCMCOperator operator, Phase phase, long start) {
        final long now = System.nanoTime();
        long[] t = stepTimes.get(operator);
        if (t == null) {
            t = new long[Phase.values().length];
            stepTimes.put(operator, t);
        }
        t[phase.ordinal()] += now - start;
        return now;
    }

    private void addEvaluation(Likelihood likelihood, long nanos, long calculations) {
        addEvaluation(likelihoodTimes, likelihood, nanos, calculations);
        final MCMCOperator operator = currentOperator;
        if (operator != null) {
            Map<Likelihood, long[]> times = likelihoodTimesByOperator.get(operator);
            if (times == null) {
                likelihoodTimesByOperator.putIfAbsent(operator, new ConcurrentHashMap<Likelihood, long[]>());
                times = likelihoodTimesByOperator.get(operator);
            }
            addEvaluation(times, likelihood, nanos, calculations);
        }
    }

    private static void addEvaluation(Map<Likelihood, long[]> times, Likelihood likelihood, long nanos, long calculations) {
        long[] t = times.get(likelihood);
        if (t == null) {
            times.putIfAbsent(likelihood, new long[3]);
            t = times.get(likelihood);
        }
        synchronized (t) {
            t[0] += nanos;
            t[1]++;
            t[2] += calculations;
        }
    }

    public void reset() {
        stepCounts.clear();
        stepTimes.clear();
        likelihoodTimes.clear();
        likelihoodTimesByOperator.clear();
        for (Runnable hook : resetHooks) {
            hook.run();
        }
    }

    /**
     * @return the report of every enabled aspect, including the registered reports
     */
    public String report() {
        final StringBuilder sb = new StringBuilder();
        if (STEPS) {
            sb.append(stepReport());
        }
        if (LIKELIHOODS) {
            sb.append(likelihoodReport());
        }
        for (Supplier<String> reporter : reporters) {
            sb.append(reporter.get());
        }
        return sb.toString();
    }

    private long getSteps(MCMCOperator operator) {
        final long[] count = stepCounts.get(operator);
        return count == null ? 0 : count[0];
    }

    private long getTotalSteps() {
        long steps = 0;
        for (long[] count : stepCounts.values()) {
            steps += count[0];
        }
        return steps;
    }

    // the steps aspect

    private String stepReport() {
        final Phase[] phases = Phase.values();
        final StringBuilder sb = new StringBuilder();
        sb.append("\nStep profile (microseconds per step, wall clock, after the adaptation delay)\n");
        sb.append(String.format("%-50s %10s", "Operator", "Steps"));
        for (Phase phase : phases) {
            sb.append(String.format(" %14s", phase.label));
        }
        sb.append(String.format(" %14s%n", "total"));
        final long[] all = new long[phases.length];
        for (Map.Entry<MCMCOperator, long[]> entry : stepTimes.entrySet()) {
            appendStepRow(sb, entry.getKey().getOperatorName(), getSteps(entry.getKey()), entry.getValue());
            for (int i = 0; i < all.length; i++) {
                all[i] += entry.getValue()[i];
            }
        }
        appendStepRow(sb, "(all operators)", getTotalSteps(), all);
        return sb.toString();
    }

    private static void appendStepRow(StringBuilder sb, String name, long steps, long[] t) {
        sb.append(String.format("%-50s %10d", name, steps));
        long total = 0;
        for (int i = 0; i < Phase.values().length; i++) {
            total += t[i];
            sb.append(String.format(" %14.2f", steps == 0 ? 0.0 : t[i] / 1.0e3 / steps));
        }
        sb.append(String.format(" %14.2f%n", steps == 0 ? 0.0 : total / 1.0e3 / steps));
    }

    // the likelihoods aspect

    private static final class Node {
        final Likelihood likelihood;
        final int depth;
        final List<Likelihood> children = new ArrayList<Likelihood>();

        Node(Likelihood likelihood, int depth) {
            this.likelihood = likelihood;
            this.depth = depth;
        }
    }

    /**
     * The likelihoods in the order of the report: a depth-first walk of the compound likelihoods from the
     * root, then any other likelihood that was evaluated (e.g. the ends of a path likelihood).
     */
    private List<Node> getLikelihoodTree() {
        final List<Node> tree = new ArrayList<Node>();
        final Set<Likelihood> seen = new HashSet<Likelihood>();
        addToTree(tree, seen, root, 0);
        for (Likelihood likelihood : likelihoodTimes.keySet()) {
            if (!seen.contains(likelihood)) {
                addToTree(tree, seen, likelihood, 0);
            }
        }
        return tree;
    }

    private static void addToTree(List<Node> tree, Set<Likelihood> seen, Likelihood likelihood, int depth) {
        if (!seen.add(likelihood)) {
            return;
        }
        final Node node = new Node(likelihood, depth);
        tree.add(node);
        if (likelihood instanceof CompoundLikelihood) {
            for (Likelihood child : ((CompoundLikelihood) likelihood).getLikelihoods()) {
                if (!seen.contains(child)) {
                    node.children.add(child);
                }
                addToTree(tree, seen, child, depth + 1);
            }
        }
    }

    private static String getLikelihoodName(Node node) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < node.depth; i++) {
            sb.append("  ");
        }
        final String id = node.likelihood.getId();
        sb.append(id != null ? id : node.likelihood.prettyName());
        return sb.toString();
    }

    /**
     * @return the likelihood's own nanoseconds: its evaluations less those of its children in the tree (a
     * compound's loop or thread pool), negative if threaded children overlapped
     */
    private static long getOwnNanos(Map<Likelihood, long[]> times, Node node) {
        final long[] t = times.get(node.likelihood);
        long own = t == null ? 0 : t[0];
        for (Likelihood child : node.children) {
            final long[] c = times.get(child);
            if (c != null) {
                own -= c[0];
            }
        }
        return own;
    }

    private String likelihoodReport() {
        final StringBuilder sb = new StringBuilder();
        final List<Node> tree = getLikelihoodTree();
        final long steps = getTotalSteps();
        final long chainNanos = STEPS ? getTotalNanos() : getInclusiveNanos(likelihoodTimes);
        sb.append("\nLikelihood profile (wall clock per evaluation, after the adaptation delay; a likelihood's calls include its components, own excludes them; % of ")
                .append(STEPS ? "the chain's step time" : "the posterior's evaluations").append(")\n");
        sb.append(String.format("%-50s %12s %12s %12s %14s %12s %12s %7s%n", "Likelihood", "Calls", "Total(s)", "us/call", "calcs/call", "us/step", "own us/step", "%"));
        for (Node node : tree) {
            final long[] t = likelihoodTimes.get(node.likelihood);
            if (t == null) {
                continue;
            }
            final long own = getOwnNanos(likelihoodTimes, node);
            sb.append(String.format("%-50s %12d %12.3f %12.2f %14.1f %12.2f %12.2f %7.2f%n", getLikelihoodName(node),
                    t[1], t[0] / 1.0e9, t[0] / 1.0e3 / t[1], (double) t[2] / t[1],
                    perStep(t[0], steps), perStep(own, steps), percent(own, chainNanos)));
        }
        if (STEPS) {
            // the phases summed over the operators, so that the own column adds up to the chain
            final long[] phases = new long[Phase.values().length];
            for (long[] t : stepTimes.values()) {
                for (int i = 0; i < phases.length; i++) {
                    phases[i] += t[i];
                }
            }
            appendOwnRow(sb, "(evaluate phase outside the likelihoods)", phases[Phase.EVALUATE.ordinal()] - getInclusiveNanos(likelihoodTimes), steps, chainNanos);
            for (Phase phase : Phase.values()) {
                if (phase != Phase.EVALUATE) {
                    appendOwnRow(sb, "(" + phase.label + " phase)", phases[phase.ordinal()], steps, chainNanos);
                }
            }
            appendOwnRow(sb, "(step total)", chainNanos, steps, chainNanos);
        }

        for (MCMCOperator operator : stepCounts.keySet()) {
            final Map<Likelihood, long[]> times = likelihoodTimesByOperator.get(operator);
            if (times == null) {
                continue;
            }
            final long operatorSteps = getSteps(operator);
            long largest = 0;
            for (long[] t : times.values()) {
                largest = Math.max(largest, t[0]);
            }
            final long[] phases = stepTimes.get(operator);
            long total = 0;
            if (phases != null) {
                for (long nanos : phases) {
                    total += nanos;
                }
            }
            final long operatorNanos = STEPS && phases != null ? total : getInclusiveNanos(times);
            sb.append(String.format("%nLikelihood profile of %s (%d steps; microseconds per step of the operator; rows under 1%% of the largest omitted; %% of %s)%n",
                    operator.getOperatorName(), operatorSteps, STEPS && phases != null ? "the operator's step time" : "the posterior's evaluations"));
            sb.append(String.format("%-50s %12s %12s %14s %12s %12s %7s%n", "Likelihood", "Calls/step", "us/call", "calcs/call", "us/step", "own us/step", "%"));
            long omittedOwn = 0;
            long allOwn = 0;
            for (Node node : tree) {
                final long[] t = times.get(node.likelihood);
                if (t == null) {
                    continue;
                }
                final long own = getOwnNanos(times, node);
                allOwn += own;
                if (t[0] * 100 < largest) {
                    omittedOwn += own;
                    continue;
                }
                sb.append(String.format("%-50s %12.2f %12.2f %14.1f %12.2f %12.2f %7.2f%n", getLikelihoodName(node),
                        operatorSteps == 0 ? 0.0 : (double) t[1] / operatorSteps, t[0] / 1.0e3 / t[1], (double) t[2] / t[1],
                        perStep(t[0], operatorSteps), perStep(own, operatorSteps), percent(own, operatorNanos)));
            }
            if (omittedOwn != 0) {
                appendOwnRow(sb, "(omitted rows)", omittedOwn, operatorSteps, operatorNanos);
            }
            if (STEPS && phases != null) {
                // the step's phases, so that the own column adds up to the step: the evaluate phase outside
                // the likelihoods is the chain's own work (the Hastings terms of the likelihoods, the checks
                // of the score; the full-evaluation test's extra evaluations count against it)
                final long evaluate = phases[Phase.EVALUATE.ordinal()];
                appendOwnRow(sb, "(evaluate phase outside the likelihoods)", evaluate - allOwn, operatorSteps, operatorNanos);
                for (Phase phase : Phase.values()) {
                    if (phase != Phase.EVALUATE) {
                        appendOwnRow(sb, "(" + phase.label + " phase)", phases[phase.ordinal()], operatorSteps, operatorNanos);
                    }
                }
                appendOwnRow(sb, "(step total)", total, operatorSteps, operatorNanos);
            }
        }
        return sb.toString();
    }

    private static void appendOwnRow(StringBuilder sb, String label, long nanos, long steps, long totalNanos) {
        sb.append(String.format("%-50s %12s %12s %14s %12s %12.2f %7.2f%n", label, "", "", "", "", perStep(nanos, steps), percent(nanos, totalNanos)));
    }

    private static double perStep(long nanos, long steps) {
        return steps == 0 ? 0.0 : nanos / 1.0e3 / steps;
    }

    private static double percent(long nanos, long totalNanos) {
        return totalNanos == 0 ? 0.0 : 100.0 * nanos / totalNanos;
    }

    /**
     * @return the nanoseconds of every phase of every step
     */
    private long getTotalNanos() {
        long total = 0;
        for (long[] phases : stepTimes.values()) {
            for (long nanos : phases) {
                total += nanos;
            }
        }
        return total;
    }

    /**
     * @return the nanoseconds of the root's evaluations in the times (the whole of the likelihood tree)
     */
    private long getInclusiveNanos(Map<Likelihood, long[]> times) {
        final long[] t = times.get(root);
        return t == null ? 0 : t[0];
    }
}
