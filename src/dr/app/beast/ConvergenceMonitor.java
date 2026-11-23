/*
 * ConvergenceMonitor.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
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

package dr.app.beast;

import dr.inference.loggers.ChainDiagnostics;
import dr.inference.mcmcmc.ChainSet;
import dr.inference.trace.TraceCorrelation;
import dr.inference.trace.TraceType;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Compares independent chains while they run (or the two halves of a single chain, in which
 * case there is no ASDSF). The chains share one burn-in, estimated from
 * their agreement: the smallest burn-in (of at least 10%) after which the split R-hat of the posterior across
 * the chains stays below 1.01 (for that burn-in and every larger one tried, as the initial
 * transient can make R-hat look small when nothing is discarded). After that burn-in it gives
 * the effective sample size of each logged value (summed over the chains), the split potential
 * scale reduction factor (R-hat) and the average standard deviation of clade frequencies
 * (ASDSF). The diagnostics are asked for by the screen log each time every chain has reached
 * the next logged state, and are also written to a log file. A final report goes to the screen.
 * This only reports; it never stops the chains.
 *
 * @author Andrew Rambaut
 */
public class ConvergenceMonitor {

    public static final double MIN_CLADE_FREQUENCY = 0.1;

    // burn-ins tried, as proportions of the samples (the step matches the tree blocks)
    public static final double BURNIN_STEP = 1.0 / SampleCollector.TREE_BLOCK_COUNT;
    // at least this much is always left out, as the posterior can settle before everything else
    public static final double MIN_BURNIN = 0.1;
    public static final double MAX_BURNIN = 0.5;
    public static final double RHAT_THRESHOLD = 1.01;

    public static final String[] ASDSF_CITATIONS = {
            "Lakner C, van der Mark P, Huelsenbeck JP, Larget B, Ronquist F (2008) Efficiency of Markov chain Monte Carlo tree proposals in Bayesian phylogenetics. Systematic Biology 57, 86-103."
    };

    public static final String[] RHAT_CITATIONS = {
            "Gelman A, Rubin DB (1992) Inference from iterative simulation using multiple sequences. Statistical Science 7, 457-472."
    };

    /**
     * @param sets the sets the collectors are for, for their status in the final report (can
     *             be null)
     */
    public ConvergenceMonitor(List<ChainSet> sets, List<SampleCollector> collectors, File logFile) {
        this.sets = sets;
        this.collectors = collectors;
        this.logFile = logFile;
    }

    public void start() {
        if (logFile != null) {
            try {
                log = new PrintWriter(new FileWriter(logFile));
                log.println("state\tburnin\tposteriorRhat\tminESS\tminESS.parameter\tmaxRhat\tmaxRhat.parameter\tASDSF\tclades");
                log.flush();
            } catch (IOException e) {
                throw new RuntimeException("Unable to write the convergence log, " + logFile + ": " + e.getMessage());
            }
        }
    }

    /**
     * Writes the final report.
     */
    public void stop() {
        report("Final convergence diagnostics");
        if (log != null) {
            log.close();
        }
    }

    /**
     * The diagnostics of the samples up to and including a state, which every chain should
     * have reached. Adds a line to the convergence log.
     */
    public synchronized Diagnostics getDiagnostics(long state) {
        SampleCollector.Snapshot[] snapshots = new SampleCollector.Snapshot[collectors.size()];
        for (int k = 0; k < snapshots.length; k++) {
            snapshots[k] = collectors.get(k).snapshot(state);
        }
        Diagnostics diagnostics = getDiagnostics(snapshots);
        writeLog(state, diagnostics);
        return diagnostics;
    }

    private Diagnostics getDiagnostics(SampleCollector.Snapshot[] snapshots) {
        Burnin burnin = estimateBurnin(snapshots);
        // without a common burn-in, give the diagnostics after the largest one tried
        double proportion = burnin.found ? burnin.proportion : MAX_BURNIN;
        Diagnostics diagnostics = computeDiagnostics(collectors.get(0).getLabels(), snapshots, proportion);
        diagnostics.burnin = burnin;
        diagnostics.burninProportion = proportion;
        return diagnostics;
    }

    private void writeLog(long state, Diagnostics diagnostics) {
        if (log != null) {
            Burnin burnin = diagnostics.burnin;
            log.println(state + "\t" + (burnin.found ? Long.toString(burnin.states) : "NA") +
                    "\t" + format(burnin.posteriorRhat) +
                    "\t" + format(diagnostics.minESS) + "\t" + format(diagnostics.minESSLabel) +
                    "\t" + format(diagnostics.maxRhat) + "\t" + format(diagnostics.maxRhatLabel) +
                    "\t" + format(diagnostics.asdsf) + "\t" + diagnostics.cladeCount);
            log.flush();
        }
    }

    private synchronized void report(String title) {
        int chainCount = collectors.size();
        if (chainCount == 0) {
            // no chain was run
            return;
        }
        SampleCollector.Snapshot[] snapshots = new SampleCollector.Snapshot[chainCount];
        long[] states = new long[chainCount];
        String[] names = new String[chainCount];
        String[] statuses = new String[chainCount];
        long minState = Long.MAX_VALUE;
        for (int k = 0; k < chainCount; k++) {
            snapshots[k] = collectors.get(k).snapshot();
            // the collector does not see the very last state
            states[k] = (sets != null ? sets.get(k).getCurrentState() : snapshots[k].state);
            names[k] = collectors.get(k).getName();
            statuses[k] = (sets != null ? sets.get(k).getStatus() : "");
            minState = Math.min(minState, states[k]);
        }

        Diagnostics diagnostics = getDiagnostics(snapshots);

        // each chain on its own, after the same burn-in
        Diagnostics[] chainDiagnostics = new Diagnostics[chainCount];
        for (int k = 0; k < chainCount; k++) {
            chainDiagnostics[k] = computeDiagnostics(collectors.get(0).getLabels(),
                    new SampleCollector.Snapshot[]{snapshots[k]}, diagnostics.burninProportion);
        }

        // print in one go so the lines are not split by other output
        System.out.print(formatReport(title, names, states, statuses, snapshots, diagnostics, chainDiagnostics));
        System.out.flush();

        writeLog(minState, diagnostics);
    }

    // what the final report asks of the diagnostics
    public static final double ESS_TARGET = 200.0;
    public static final double ESS_BORDERLINE = 100.0;
    public static final double RHAT_BORDERLINE = 1.05;
    public static final double ASDSF_GOOD = 0.01;
    public static final double ASDSF_ACCEPTABLE = 0.05;

    private static final String GOOD = "good";
    private static final String BORDERLINE = "borderline";
    private static final String POOR = "poor";
    private static final String UNKNOWN = "-";

    private static String assessRhat(double rhat) {
        if (Double.isNaN(rhat)) {
            return UNKNOWN;
        }
        return rhat < RHAT_THRESHOLD ? GOOD : (rhat < RHAT_BORDERLINE ? BORDERLINE : POOR);
    }

    private static String assessESS(double ess) {
        if (Double.isNaN(ess)) {
            return UNKNOWN;
        }
        return ess >= ESS_TARGET ? GOOD : (ess >= ESS_BORDERLINE ? BORDERLINE : POOR);
    }

    private static String assessASDSF(double asdsf) {
        if (Double.isNaN(asdsf)) {
            return UNKNOWN;
        }
        return asdsf < ASDSF_GOOD ? GOOD : (asdsf < ASDSF_ACCEPTABLE ? "acceptable" : POOR);
    }

    private static String number(double value, String format) {
        return Double.isNaN(value) ? "-" : String.format(format, value);
    }

    private static String name(String label) {
        return label == null ? "" : label.trim();
    }

    /**
     * @return the number of states rounded up to two significant figures
     */
    static long roundUp(double states) {
        if (states <= 0) {
            return 0;
        }
        double unit = Math.pow(10, Math.floor(Math.log10(states)) - 1);
        return (long) (Math.ceil(states / unit - 1E-9) * unit);
    }

    static String formatReport(String title, String[] names, long[] states, String[] statuses,
                               SampleCollector.Snapshot[] snapshots,
                               Diagnostics diagnostics, Diagnostics[] chainDiagnostics) {
        final int chainCount = names.length;
        final Burnin burnin = diagnostics.burnin;
        final String nl = String.format("%n");

        long chainLength = Long.MAX_VALUE;
        int sampleCount = Integer.MAX_VALUE;
        for (int k = 0; k < chainCount; k++) {
            chainLength = Math.min(chainLength, states[k]);
            sampleCount = Math.min(sampleCount, snapshots[k].rows.size());
        }
        final int burninSamples = getBurninSamples(sampleCount, diagnostics.burninProportion);
        final long burninStates = (burnin.found ? burnin.states : burninSamples * snapshots[0].logEvery);
        final String all = (chainCount > 1 ? "every chain" : "the chain");

        StringBuilder sb = new StringBuilder();
        sb.append(nl).append(String.format("%s at state %d (%d chain%s)", title, chainLength, chainCount, chainCount > 1 ? "s" : "")).append(nl);

        // the burn-in
        sb.append(nl).append("Burn-in").append(nl);
        if (burnin.found) {
            sb.append(String.format("  The first %d states (%.0f%%) of %s are left out of everything below.",
                    burninStates, burnin.proportion * 100, all)).append(nl);
            sb.append(String.format("  After that the posterior%s (R-hat %.4f, which should be below %.2f).",
                    chainCount > 1 ? "s of the chains agree" : " is the same in the two halves of the chain",
                    burnin.posteriorRhat, RHAT_THRESHOLD)).append(nl);
        } else if (burnin.tried) {
            sb.append(String.format("  No burn-in of up to %.0f%% makes the posterior%s (R-hat %.4f after %.0f%%, which should be below %.2f).",
                    MAX_BURNIN * 100, chainCount > 1 ? "s of the chains agree" : " the same in the two halves of the chain",
                    burnin.posteriorRhat, MAX_BURNIN * 100, RHAT_THRESHOLD)).append(nl);
            sb.append(String.format("  The first %d states (%.0f%%) of %s are left out of everything below.",
                    burninStates, MAX_BURNIN * 100, all)).append(nl);
        } else {
            sb.append(String.format("  There are too few samples to estimate a burn-in. The first %.0f%% of %s is left out of everything below.",
                    MAX_BURNIN * 100, all)).append(nl);
        }

        // each chain
        int nameWidth = 9;
        int essWidth = 9;
        for (int k = 0; k < chainCount; k++) {
            nameWidth = Math.max(nameWidth, names[k].length());
            essWidth = Math.max(essWidth, name(chainDiagnostics[k].minESSLabel).length());
        }
        String rowFormat = "  %-" + nameWidth + "s  %-12s  %-14s  %-8s  %-8s  %-" + essWidth + "s  %-9s  %s";
        sb.append(nl).append(chainCount > 1 ? "Each chain on its own, after the burn-in" : "The chain, after the burn-in").append(nl);
        sb.append(String.format(rowFormat, "chain", "state", "posterior", "samples", "min ESS", "of", "max R-hat", "of").replaceAll("\\s+$", "")).append(nl);
        double lowestChainESS = Double.NaN;
        for (int k = 0; k < chainCount; k++) {
            Diagnostics chain = chainDiagnostics[k];
            String row = String.format(rowFormat, names[k], Long.toString(states[k]),
                    String.format("%.4f", snapshots[k].posterior),
                    Integer.toString(Math.max(0, snapshots[k].rows.size() - getBurninSamples(snapshots[k].rows.size(), diagnostics.burninProportion))),
                    number(chain.minESS, "%.0f"), name(chain.minESSLabel),
                    number(chain.maxRhat, "%.3f"), name(chain.maxRhatLabel));
            if (!statuses[k].isEmpty()) {
                row += "   " + statuses[k];
            }
            sb.append(row.replaceAll("\\s+$", "")).append(nl);
            if (!Double.isNaN(chain.minESS) && !(chain.minESS >= lowestChainESS)) {
                lowestChainESS = chain.minESS;
            }
        }
        sb.append("  Within a chain, R-hat compares the first half of the chain with the second.").append(nl);

        // the chains together
        final String posteriorAssessment = assessRhat(burnin.posteriorRhat);
        final String rhatAssessment = assessRhat(diagnostics.maxRhat);
        final String essAssessment = assessESS(diagnostics.minESS);
        final String asdsfAssessment = assessASDSF(diagnostics.asdsf);
        final boolean hasClades = chainCount > 1 && diagnostics.cladeCount > 0;

        int ofWidth = Math.max(10, Math.max(name(diagnostics.minESSLabel).length(), name(diagnostics.maxRhatLabel).length()));
        String tableFormat = "  %-15s  %-8s  %-" + ofWidth + "s  %-28s  %s";
        if (chainCount > 1) {
            sb.append(nl).append("All chains together, after the burn-in").append(nl);
            sb.append(String.format(tableFormat, "statistic", "value", "of", "wanted", "assessment")).append(nl);
            sb.append(String.format(tableFormat, "posterior R-hat", number(burnin.posteriorRhat, "%.4f"), "",
                    String.format("below %.2f", RHAT_THRESHOLD), posteriorAssessment)).append(nl);
            sb.append(String.format(tableFormat, "max R-hat", number(diagnostics.maxRhat, "%.4f"), name(diagnostics.maxRhatLabel),
                    String.format("below %.2f", RHAT_THRESHOLD), rhatAssessment)).append(nl);
            sb.append(String.format(tableFormat, "min ESS", number(diagnostics.minESS, "%.0f"), name(diagnostics.minESSLabel),
                    String.format("at least %.0f", ESS_TARGET), essAssessment)).append(nl);
            if (hasClades) {
                sb.append(String.format(tableFormat, "ASDSF", number(diagnostics.asdsf, "%.4f"), diagnostics.cladeCount + " clades",
                        String.format("below %.2f (%.2f acceptable)", ASDSF_GOOD, ASDSF_ACCEPTABLE), asdsfAssessment)).append(nl);
            }
            sb.append(nl);
            sb.append("  R-hat compares the chains with each other, and the two halves of each chain. Close to 1 means").append(nl);
            sb.append("  they are sampling the same distribution; larger means they are not, or not yet. The largest over").append(nl);
            sb.append("  everything in the log is given, with what it is of.").append(nl);
            sb.append("  ESS is the effective number of independent samples, added up over the chains. The smallest").append(nl);
            sb.append("  over everything in the log is given. It only means something if the chains agree.").append(nl);
            if (hasClades) {
                sb.append(String.format("  ASDSF is how much the chains differ in how often they sample each clade, averaged over the clades%n" +
                        "  in at least %.0f%% of some chain's trees. Close to 0 means the chains sample the same trees.", MIN_CLADE_FREQUENCY * 100)).append(nl);
            }
        }

        // what to do
        sb.append(nl).append("Recommendations").append(nl);

        final boolean agree = GOOD.equals(posteriorAssessment) && GOOD.equals(rhatAssessment) &&
                (!hasClades || !POOR.equals(asdsfAssessment));
        final boolean nearly = !agree && !POOR.equals(posteriorAssessment) && !UNKNOWN.equals(posteriorAssessment) &&
                !POOR.equals(rhatAssessment) && !UNKNOWN.equals(rhatAssessment) &&
                (!hasClades || !POOR.equals(asdsfAssessment));

        List<String> problems = new ArrayList<String>();
        if (!GOOD.equals(posteriorAssessment) && !UNKNOWN.equals(posteriorAssessment)) {
            problems.add(String.format("the R-hat of the posterior is %.3f", burnin.posteriorRhat));
        }
        if (!GOOD.equals(rhatAssessment) && !UNKNOWN.equals(rhatAssessment)) {
            problems.add(String.format("the R-hat of %s is %.3f", name(diagnostics.maxRhatLabel), diagnostics.maxRhat));
        }
        if (hasClades && POOR.equals(asdsfAssessment)) {
            problems.add(String.format("the ASDSF is %.3f", diagnostics.asdsf));
        }
        String because = "";
        for (int i = 0; i < problems.size(); i++) {
            because += (i == 0 ? "" : (i == problems.size() - 1 ? " and " : ", ")) + problems.get(i);
        }

        if (!burnin.tried || UNKNOWN.equals(rhatAssessment)) {
            sb.append("  There are too few samples to say whether " + (chainCount > 1 ? "the chains agree" : "the chain has settled") + ". Run for longer.").append(nl);
        } else if (chainCount > 1) {
            if (agree) {
                sb.append(String.format("  Combining: the chains agree, so their samples can be combined after removing the first %d%n" +
                        "  states of each (in LogCombiner, a burn-in of %d states).", burninStates, burninStates)).append(nl);
            } else if (nearly) {
                sb.append(String.format("  Combining: the chains nearly agree (%s). Their samples can be combined after%n" +
                        "  removing the first %d states of each, but running for longer would be safer.", because, burninStates)).append(nl);
            } else {
                sb.append(String.format("  Combining: the chains should not be combined: %s.%n" +
                        "  They are not sampling the same distribution, or not yet. Run them for longer; if they still%n" +
                        "  disagree, the posterior may have separate modes that the chains do not move between.", because)).append(nl);
            }
        } else {
            if (agree) {
                sb.append("  The two halves of the chain agree, so there is no sign that it is still settling.").append(nl);
                sb.append("  A single chain cannot show whether another start would end up in the same place: for that run several (-chains).").append(nl);
            } else {
                sb.append(String.format("  The two halves of the chain %s (%s), so it %s. Run for longer.",
                        nearly ? "nearly agree" : "differ", because, nearly ? "may still be settling" : "has not settled")).append(nl);
            }
        }

        if (burnin.tried && !UNKNOWN.equals(rhatAssessment) && !Double.isNaN(diagnostics.minESS)) {
            final double sampled = chainLength - burninStates;
            if (!agree && !nearly) {
                sb.append("  Chain length: the effective sample sizes cannot be trusted until " +
                        (chainCount > 1 ? "the chains agree" : "the chain has settled") + ", so they give no" + nl +
                        "  chain length. As a start, run for at least twice as long.").append(nl);
            } else if (diagnostics.minESS >= ESS_TARGET) {
                String line = String.format("  Chain length: enough. The smallest ESS is %.0f (%s)%s, above the%n  %.0f usually wanted.",
                        diagnostics.minESS, name(diagnostics.minESSLabel), chainCount > 1 ? " with the chains combined" : "", ESS_TARGET);
                sb.append(line).append(nl);
                if (chainCount > 1 && !Double.isNaN(lowestChainESS)) {
                    if (lowestChainESS >= ESS_TARGET) {
                        sb.append(String.format("  Each chain would also be enough on its own (the smallest ESS in any one chain is %.0f).", lowestChainESS)).append(nl);
                    } else {
                        sb.append(String.format("  No single chain would be enough on its own (the smallest ESS in any one chain is %.0f).", lowestChainESS)).append(nl);
                    }
                }
            } else if (sampled > 0 && diagnostics.minESS > 0) {
                final double factor = ESS_TARGET / diagnostics.minESS;
                final long recommended = roundUp(burninStates + sampled * factor);
                sb.append(String.format("  Chain length: too short. The smallest ESS is %.0f (%s)%s, below the%n  %.0f usually wanted.",
                        diagnostics.minESS, name(diagnostics.minESSLabel), chainCount > 1 ? " with the chains combined" : "", ESS_TARGET)).append(nl);
                sb.append(String.format("  ESS grows in proportion to the states sampled after the burn-in, so to reach %.0f run %s for about%n" +
                        "  %d states (%.1f times the sampling after the burn-in)%s.",
                        ESS_TARGET, chainCount > 1 ? "each chain" : "the chain", recommended, factor,
                        chainCount > 1 ? String.format(", or run about %d chains of this length", (int) Math.ceil(chainCount * factor - 1E-9)) : "")).append(nl);
            }
        }

        sb.append(nl).append("References").append(nl);
        sb.append("  R-hat, the split potential scale reduction factor:").append(nl);
        for (String citation : RHAT_CITATIONS) {
            sb.append("    ").append(citation).append(nl);
        }
        if (hasClades) {
            sb.append("  ASDSF, the average standard deviation of clade frequencies between the chains:").append(nl);
            for (String citation : ASDSF_CITATIONS) {
                sb.append("    ").append(citation).append(nl);
            }
        }
        return sb.toString();
    }

    private static String format(double value) {
        return Double.isNaN(value) ? "NA" : Double.toString(value);
    }

    private static String format(String label) {
        return label == null ? "NA" : label;
    }

    static class Burnin {
        boolean tried = false;
        boolean found = false;
        double proportion = Double.NaN;
        long states = 0;
        // at the burn-in found, or after MAX_BURNIN if none was
        double posteriorRhat = Double.NaN;
    }

    /**
     * Finds the smallest burn-in (from MIN_BURNIN to MAX_BURNIN in steps of BURNIN_STEP) for which the
     * split R-hat of the chains' posteriors is below RHAT_THRESHOLD, for that burn-in and every
     * larger one tried.
     */
    static Burnin estimateBurnin(SampleCollector.Snapshot[] snapshots) {
        Burnin burnin = new Burnin();
        int n = Integer.MAX_VALUE;
        for (SampleCollector.Snapshot snapshot : snapshots) {
            n = Math.min(n, snapshot.scores.size());
        }
        int steps = (int) Math.round(MAX_BURNIN / BURNIN_STEP);
        // need at least a few samples in each half chain after the largest burn-in
        if (n * (1.0 - MAX_BURNIN) < 8) {
            return burnin;
        }
        burnin.tried = true;

        double[] rhats = new double[steps + 1];
        for (int i = 0; i <= steps; i++) {
            int start = getBurninSamples(n, i * BURNIN_STEP);
            double[][] chains = new double[snapshots.length][n - start];
            for (int k = 0; k < snapshots.length; k++) {
                for (int s = start; s < n; s++) {
                    chains[k][s - start] = snapshots[k].scores.get(s);
                }
            }
            rhats[i] = splitRhat(chains);
        }

        int found = -1;
        int minStep = (int) Math.round(MIN_BURNIN / BURNIN_STEP);
        for (int i = steps; i >= minStep; i--) {
            if (!(rhats[i] < RHAT_THRESHOLD)) {
                break;
            }
            found = i;
        }
        if (found >= 0) {
            burnin.found = true;
            burnin.proportion = found * BURNIN_STEP;
            burnin.states = getBurninSamples(n, burnin.proportion) * snapshots[0].logEvery;
            burnin.posteriorRhat = rhats[found];
        } else {
            burnin.posteriorRhat = rhats[steps];
        }
        return burnin;
    }

    public static int getBurninSamples(int sampleCount, double proportion) {
        // round first so steps like 3 * 0.02 do not fall just above a whole number
        return (int) Math.ceil(Math.round(sampleCount * proportion * 1e6) / 1e6);
    }

    static class Diagnostics implements ChainDiagnostics {
        double minESS = Double.NaN;
        String minESSLabel = null;
        double maxRhat = Double.NaN;
        String maxRhatLabel = null;
        double asdsf = Double.NaN;
        int cladeCount = 0;
        Burnin burnin = new Burnin();
        // the burn-in the diagnostics are after: the common one, or MAX_BURNIN if there is none
        double burninProportion = MAX_BURNIN;

        @Override
        public double getMinESS() {
            return minESS;
        }

        @Override
        public String getMinESSLabel() {
            return minESSLabel;
        }

        @Override
        public double getMaxRhat() {
            return maxRhat;
        }

        @Override
        public String getMaxRhatLabel() {
            return maxRhatLabel;
        }

        @Override
        public double getASDSF() {
            return asdsf;
        }

        @Override
        public double getBurnin() {
            return burnin.found ? burnin.states : Double.NaN;
        }

        @Override
        public double getBurninProportion() {
            return burninProportion;
        }
    }

    static Diagnostics computeDiagnostics(String[] labels, SampleCollector.Snapshot[] snapshots, double burnin) {
        Diagnostics diagnostics = new Diagnostics();
        int chainCount = snapshots.length;

        // the same number of samples from each chain, after the burn-in
        int sampleCount = Integer.MAX_VALUE;
        for (SampleCollector.Snapshot snapshot : snapshots) {
            sampleCount = Math.min(sampleCount, snapshot.rows.size());
        }
        int start = getBurninSamples(sampleCount, burnin);
        int n = sampleCount - start;

        if (n >= 4) {
            for (int j = 0; j < labels.length; j++) {
                double[][] values = new double[chainCount][n];
                boolean constant = true;
                for (int k = 0; k < chainCount; k++) {
                    for (int s = 0; s < n; s++) {
                        values[k][s] = snapshots[k].rows.get(start + s)[j];
                        if (values[k][s] != values[0][0]) {
                            constant = false;
                        }
                    }
                }
                if (constant) {
                    continue;
                }

                double ess = 0.0;
                for (int k = 0; k < chainCount; k++) {
                    ess += getESS(values[k], snapshots[k].logEvery);
                }
                if (!Double.isNaN(ess) && (diagnostics.minESSLabel == null || ess < diagnostics.minESS)) {
                    diagnostics.minESS = ess;
                    diagnostics.minESSLabel = labels[j];
                }

                // also for a single chain, whose two halves are compared
                {
                    double rhat = splitRhat(values);
                    if (!Double.isNaN(rhat) && (diagnostics.maxRhatLabel == null || rhat > diagnostics.maxRhat)) {
                        diagnostics.maxRhat = rhat;
                        diagnostics.maxRhatLabel = labels[j];
                    }
                }
            }
        }

        if (chainCount > 1) {
            computeASDSF(snapshots, burnin, diagnostics);
        }
        return diagnostics;
    }

    public static double getESS(double[] values, long stepSize) {
        List<Double> list = new ArrayList<Double>(values.length);
        for (double value : values) {
            list.add(value);
        }
        return new TraceCorrelation(list, TraceType.REAL, Math.max(1, stepSize)).getESS();
    }

    /**
     * The split potential scale reduction factor (Gelman et al., Bayesian Data Analysis, 3rd
     * edition): each chain is split in half, so that a trend within a chain is detected as well
     * as disagreement between chains.
     */
    public static double splitRhat(double[][] chains) {
        int half = chains[0].length / 2;
        if (half < 2) {
            return Double.NaN;
        }
        int m = chains.length * 2;
        double[] means = new double[m];
        double[] variances = new double[m];
        int c = 0;
        for (double[] chain : chains) {
            for (int part = 0; part < 2; part++) {
                int offset = part * half;
                double mean = 0.0;
                for (int s = 0; s < half; s++) {
                    mean += chain[offset + s];
                }
                mean /= half;
                double variance = 0.0;
                for (int s = 0; s < half; s++) {
                    double d = chain[offset + s] - mean;
                    variance += d * d;
                }
                variance /= (half - 1);
                means[c] = mean;
                variances[c] = variance;
                c++;
            }
        }

        double grandMean = 0.0;
        double w = 0.0;
        for (int i = 0; i < m; i++) {
            grandMean += means[i];
            w += variances[i];
        }
        grandMean /= m;
        w /= m;

        double b = 0.0;
        for (int i = 0; i < m; i++) {
            double d = means[i] - grandMean;
            b += d * d;
        }
        b *= (double) half / (m - 1);

        if (w <= 0.0) {
            return Double.NaN;
        }
        double varPlus = ((half - 1.0) / half) * w + b / half;
        return Math.sqrt(varPlus / w);
    }

    /**
     * The average standard deviation of clade frequencies across the chains, over the clades
     * with a frequency of at least MIN_CLADE_FREQUENCY in some chain.
     */
    static void computeASDSF(SampleCollector.Snapshot[] snapshots, double burnin, Diagnostics diagnostics) {
        int chainCount = snapshots.length;
        List<Map<Long, Double>> frequencies = new ArrayList<Map<Long, Double>>();
        for (SampleCollector.Snapshot snapshot : snapshots) {
            int total = 0;
            for (int size : snapshot.treeBlockSizes) {
                total += size;
            }
            if (total == 0) {
                return;
            }
            // drop whole blocks up to the burn-in
            int skip = 0;
            int skipped = 0;
            while (skip < snapshot.treeBlockSizes.size() && skipped + snapshot.treeBlockSizes.get(skip) <= total * burnin + 1e-6) {
                skipped += snapshot.treeBlockSizes.get(skip);
                skip++;
            }
            int count = total - skipped;
            if (count == 0) {
                return;
            }
            Map<Long, Double> freq = new HashMap<Long, Double>();
            for (int b = skip; b < snapshot.treeBlocks.size(); b++) {
                for (Map.Entry<Long, Integer> entry : snapshot.treeBlocks.get(b).entrySet()) {
                    Double f = freq.get(entry.getKey());
                    freq.put(entry.getKey(), (f == null ? 0.0 : f) + entry.getValue());
                }
            }
            for (Map.Entry<Long, Double> entry : freq.entrySet()) {
                entry.setValue(entry.getValue() / count);
            }
            frequencies.add(freq);
        }

        Set<Long> clades = new HashSet<Long>();
        for (Map<Long, Double> freq : frequencies) {
            for (Map.Entry<Long, Double> entry : freq.entrySet()) {
                if (entry.getValue() >= MIN_CLADE_FREQUENCY) {
                    clades.add(entry.getKey());
                }
            }
        }
        if (clades.isEmpty()) {
            return;
        }

        double sum = 0.0;
        for (Long clade : clades) {
            double mean = 0.0;
            double[] f = new double[chainCount];
            for (int k = 0; k < chainCount; k++) {
                Double value = frequencies.get(k).get(clade);
                f[k] = (value == null ? 0.0 : value);
                mean += f[k];
            }
            mean /= chainCount;
            double variance = 0.0;
            for (int k = 0; k < chainCount; k++) {
                variance += (f[k] - mean) * (f[k] - mean);
            }
            variance /= (chainCount - 1);
            sum += Math.sqrt(variance);
        }
        diagnostics.asdsf = sum / clades.size();
        diagnostics.cladeCount = clades.size();
    }

    private final List<ChainSet> sets;
    private final List<SampleCollector> collectors;
    private final File logFile;

    private PrintWriter log = null;
}
