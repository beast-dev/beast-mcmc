/*
 * MATProfile.java
 *
 * Copyright © 2026, the BEAST Development Team.
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
 */

package dr.evomodel.treedatalikelihood.MAT;

import dr.inference.markovchain.ChainProfile;

import java.util.function.Supplier;

/**
 * Wall-clock time of the phases of the MAT's work (the resampler's moves, the store's scans and the
 * likelihood's evaluation) and counts of what they handle, summed over the run. Enabled with
 * the {@code mat} aspect of the {@link ChainProfile} (-Dmcmc.profile=mat, or steps,mat for the per-step
 * figures); otherwise {@link #start} and {@link #stop} cost nothing. The report is printed with the chain's
 * profile at the end of the chain. Phases nest: the description
 * of each says which phase it is part of, so that times can be subtracted.
 *
 * The timings are static: with several MAT likelihoods in one run they are summed.
 *
 * @author Andrew Rambaut
 */
public final class MATProfile {

    public static final boolean ENABLED = ChainProfile.isEnabled("mat");

    public enum Phase {
        // a subtree move with the local block
        LOCAL_BEFORE("local.before", "old phase of a move with the local block (whole)"),
        LOCAL_RUNS("local.runs", "  ambiguous runs of a node recomputed (in local.before/after, gibbs)"),
        LOCAL_CHAINS("local.chains", "  old chain and candidate sites (in local.before)"),
        LOCAL_SCORE_OLD("local.scoreOld", "  scoring the old states of P and the old chain (in local.before)"),
        LOCAL_AFTER("local.after", "new phase of a move with the local block (whole)"),
        LOCAL_NEW_CHAIN("local.newChain", "  new chain and clearing the old chain's states (in local.after)"),
        LOCAL_SAMPLE_NEW("local.sampleNew", "  drawing the new states of P and the new chain (in local.after)"),
        // a subtree move with the path block
        PATH_BEFORE("path.before", "old phase of a move with the path block (whole)"),
        PATH_CANDIDATES("path.candidates", "  block and candidate sites (in path.before)"),
        PATH_SCORE_OLD("path.scoreOld", "  scoring the old states of the block (in path.before)"),
        PATH_AFTER("path.after", "new phase of a move with the path block (whole)"),
        PATH_SAMPLE_NEW("path.sampleNew", "  drawing the new states of the block (in path.after)"),
        PATH_HISTORIES("path.histories", "  scoring and drawing explicit hits (in path.before/after)"),
        // Gibbs update of a block in place
        GIBBS("gibbs", "resampleInPlace (whole; mutationHistoryGibbsOperator)"),
        // branch events
        BRANCHES("branches", "updateBranchesAroundBlock (whole; in local.after, path.after, gibbs)"),
        STORE_UPDATE_BRANCH("store.updateBranch", "  full scan of a reparented branch (in branches, parsimony)"),
        STORE_UPDATE_BRANCH_SITES("store.updateBranchSites", "  update of a branch at the changed sites (in branches)"),
        // the store
        STORE_JOURNAL_CLONE("store.journalClone", "copy-on-write clone of a node's states (in any phase that writes)"),
        STORE_FITCH_CLONE("store.fitchClone", "copy-on-write clone of a node's Fitch sets (in parsimony)"),
        STORE_SET_CATEGORY("store.setCategory", "change of a site's rate category (in categories.resample)"),
        STORE_COMMIT("store.commit", "commit of the journal (in the chain's store phase)"),
        STORE_REVERT("store.revert", "revert of the journal (in the chain's restore phase)"),
        // the likelihood
        LIK_UPDATE_MAT("lik.updateMAT", "updateMAT: topology check (in lik.calculate and the resampler)"),
        PARSIMONY("parsimony", "  parsimony update of the MAT (in lik.updateMAT)"),
        LIK_ENSURE_RATES("lik.ensureRates", "rates brought up to date (in lik.calculate and the resampler)"),
        LIK_CALCULATE("lik.calculate", "calculateLikelihood (whole; the chain's evaluate phase)"),
        LIK_BRANCH_TERMS("lik.branchTerms", "  recomputing the dirty branch terms (in lik.calculate)"),
        LIK_SUM("lik.sum", "  sum over all branch terms (in lik.calculate)"),
        LIK_SITE_LOG_LIKELIHOODS("lik.siteLogLikelihoods", "per-site log likelihoods (loggers)"),
        CATEGORIES_RESAMPLE("categories.resample", "Gibbs draw of the site rate categories (whole)");

        final String label;
        final String description;

        Phase(String label, String description) {
            this.label = label;
            this.description = description;
        }
    }

    public enum Counter {
        LOCAL_MOVES("local moves"),
        PATH_MOVES("path moves"),
        GIBBS_UPDATES("Gibbs block updates"),
        CANDIDATE_SITES("candidate sites (events around the block, local: substitutions + hot sites)"),
        CHANGED_SITES("sites whose states changed (per move or Gibbs update)"),
        HOT_SITES("hot sites of the local block (old + new chains)"),
        BRANCH_SCAN_PATTERNS("patterns scanned by full branch updates"),
        BRANCH_SITE_UPDATES("sites handled by partial branch updates"),
        JOURNAL_CLONE_BYTES("bytes cloned by copy-on-write of states"),
        RUNS_RECOMPUTED("nodes whose ambiguous runs were recomputed"),
        TERMS_RECOMPUTED("branch terms recomputed"),
        EVALUATIONS("likelihood evaluations"),
        TRANSITION_CACHE_MISSES("transition probability cache misses"),
        TRANSITION_CACHE_CLEARS("transition probability cache clears"),
        PARSIMONY_FULL("full parsimony rebuilds"),
        PARSIMONY_INCREMENTAL("incremental parsimony updates"),
        CATEGORY_CHANGES("site rate category changes");

        final String description;

        Counter(String description) {
            this.description = description;
        }
    }

    private static final long[] nanos = new long[Phase.values().length];
    private static final long[] calls = new long[Phase.values().length];
    private static final long[] counts = new long[Counter.values().length];

    // the same, attributed to the operator whose step was in progress (see ChainProfile.getCurrentOperatorName)
    private static final String NO_OPERATOR = "(no operator)";
    private static final java.util.Map<String, long[][]> byOperator = new java.util.LinkedHashMap<String, long[][]>();
    private static String lastOperator = null;
    private static long[][] lastOperatorTimes = null;

    static {
        if (ENABLED) {
            ChainProfile.addReporter(new Supplier<String>() {
                @Override
                public String get() {
                    return report();
                }
            });
            ChainProfile.addResetHook(new Runnable() {
                @Override
                public void run() {
                    reset();
                }
            });
        }
    }

    /**
     * @return the {nanos, calls, counts} arrays of the operator whose step is in progress
     */
    private static long[][] operatorTimes() {
        String operator = ChainProfile.getCurrentOperatorName();
        if (operator == null) {
            operator = NO_OPERATOR;
        }
        if (!operator.equals(lastOperator)) {
            long[][] times = byOperator.get(operator);
            if (times == null) {
                times = new long[][]{new long[Phase.values().length], new long[Phase.values().length],
                        new long[Counter.values().length]};
                byOperator.put(operator, times);
            }
            lastOperator = operator;
            lastOperatorTimes = times;
        }
        return lastOperatorTimes;
    }

    private MATProfile() {
    }

    /**
     * @return the start time of a phase (0 if profiling is off)
     */
    public static long start() {
        return ENABLED ? System.nanoTime() : 0L;
    }

    /**
     * Adds the time since start to the phase.
     */
    public static void stop(Phase phase, long start) {
        if (ENABLED) {
            final long elapsed = System.nanoTime() - start;
            nanos[phase.ordinal()] += elapsed;
            calls[phase.ordinal()]++;
            final long[][] times = operatorTimes();
            times[0][phase.ordinal()] += elapsed;
            times[1][phase.ordinal()]++;
        }
    }

    public static void count(Counter counter, long n) {
        if (ENABLED) {
            counts[counter.ordinal()] += n;
            operatorTimes()[2][counter.ordinal()] += n;
        }
    }

    public static void reset() {
        java.util.Arrays.fill(nanos, 0L);
        java.util.Arrays.fill(calls, 0L);
        java.util.Arrays.fill(counts, 0L);
        byOperator.clear();
        lastOperator = null;
        lastOperatorTimes = null;
    }

    public static String report() {
        final StringBuilder sb = new StringBuilder();
        sb.append("\nMAT profile (wall clock; the steps of the step profile, i.e. after the adaptation delay)\n");
        sb.append(String.format("%-26s %12s %12s %12s  %s%n", "Phase", "Calls", "Total(s)", "us/call", "Description"));
        for (Phase phase : Phase.values()) {
            final int i = phase.ordinal();
            if (calls[i] == 0) {
                continue;
            }
            sb.append(String.format("%-26s %12d %12.3f %12.2f  %s%n", phase.label, calls[i], nanos[i] / 1.0e9,
                    nanos[i] / 1.0e3 / calls[i], phase.description));
        }
        sb.append(String.format("%n%-50s %14s%n", "Counter", "Total"));
        for (Counter counter : Counter.values()) {
            sb.append(String.format("%-50s %14d  %s%n", counter.name().toLowerCase(), counts[counter.ordinal()],
                    counter.description));
        }
        for (java.util.Map.Entry<String, long[][]> entry : byOperator.entrySet()) {
            final String operator = entry.getKey();
            final long steps = ChainProfile.getSteps(operator);
            final long[] opNanos = entry.getValue()[0];
            final long[] opCalls = entry.getValue()[1];
            final long[] opCounts = entry.getValue()[2];
            sb.append(String.format("%nMAT profile of %s (%d steps; microseconds per step of the operator)%n", operator, steps));
            sb.append(String.format("%-26s %12s %12s %12s %12s%n", "Phase", "Calls", "Total(s)", "us/call", "us/step"));
            for (Phase phase : Phase.values()) {
                final int i = phase.ordinal();
                if (opCalls[i] == 0) {
                    continue;
                }
                sb.append(String.format("%-26s %12d %12.3f %12.2f %12.2f%n", phase.label, opCalls[i], opNanos[i] / 1.0e9,
                        opNanos[i] / 1.0e3 / opCalls[i], steps == 0 ? 0.0 : opNanos[i] / 1.0e3 / steps));
            }
            sb.append(String.format("%-50s %14s %12s%n", "Counter", "Total", "per step"));
            for (Counter counter : Counter.values()) {
                final int i = counter.ordinal();
                if (opCounts[i] == 0) {
                    continue;
                }
                sb.append(String.format("%-50s %14d %12.2f%n", counter.name().toLowerCase(), opCounts[i],
                        steps == 0 ? 0.0 : (double) opCounts[i] / steps));
            }
        }
        return sb.toString();
    }
}
