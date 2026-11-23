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
        long minState = Long.MAX_VALUE;
        for (int k = 0; k < chainCount; k++) {
            snapshots[k] = collectors.get(k).snapshot();
            minState = Math.min(minState, snapshots[k].state);
        }

        Diagnostics diagnostics = getDiagnostics(snapshots);
        Burnin burnin = diagnostics.burnin;

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%n%s at state %d (%d chain%s):%n", title, minState, chainCount, chainCount > 1 ? "s" : ""));
        for (int k = 0; k < chainCount; k++) {
            String status = (sets != null ? sets.get(k).getStatus() : "");
            sb.append(String.format("  %-10s state %-12d posterior %-16.4f %s%n",
                    collectors.get(k).getName(), snapshots[k].state, snapshots[k].posterior, status));
        }
        if (burnin.found) {
            sb.append(String.format("  common burn-in %d states (%.0f%%): posterior R-hat %.4f%n",
                    burnin.states, burnin.proportion * 100, burnin.posteriorRhat));
        } else if (burnin.tried) {
            sb.append(String.format("  no common burn-in yet: posterior R-hat %.4f after %.0f%% (needs < %.2f); below is after %.0f%%%n",
                    burnin.posteriorRhat, MAX_BURNIN * 100, RHAT_THRESHOLD, MAX_BURNIN * 100));
        } else {
            sb.append(String.format("  too few samples to estimate the burn-in (below is after %.0f%%)%n", MAX_BURNIN * 100));
        }
        if (diagnostics.minESSLabel != null) {
            sb.append(String.format("  min ESS %.1f (%s)", diagnostics.minESS, diagnostics.minESSLabel));
        } else {
            sb.append("  min ESS -");
        }
        if (diagnostics.maxRhatLabel != null) {
            sb.append(String.format("   max R-hat %.4f (%s)", diagnostics.maxRhat, diagnostics.maxRhatLabel));
        } else {
            sb.append("   max R-hat -");
        }
        if (diagnostics.cladeCount > 0) {
            sb.append(String.format("   ASDSF %.4f (%d clades)", diagnostics.asdsf, diagnostics.cladeCount));
        }
        sb.append(String.format("%n"));
        sb.append(String.format("  R-hat is the split potential scale reduction factor:%n"));
        for (String citation : RHAT_CITATIONS) {
            sb.append(String.format("    %s%n", citation));
        }

        // print in one go so the lines are not split by other output
        System.out.print(sb);
        System.out.flush();

        writeLog(minState, diagnostics);
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
