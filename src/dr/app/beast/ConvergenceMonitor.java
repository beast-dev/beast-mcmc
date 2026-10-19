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

import dr.inference.mcmcmc.ChainSet;
import dr.inference.trace.TraceCorrelation;
import dr.inference.trace.TraceType;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.TimeUnit;

/**
 * Periodically compares independent chains while they run: the effective sample size of each
 * logged value (summed over the chains), the split potential scale reduction factor (R-hat)
 * and the average standard deviation of clade frequencies (ASDSF). Reports go to the screen
 * and to a log file. This only reports; it never stops the chains.
 *
 * @author Andrew Rambaut
 */
public class ConvergenceMonitor {

    public static final double DEFAULT_BURNIN = 0.1;
    public static final double MIN_CLADE_FREQUENCY = 0.1;

    public ConvergenceMonitor(List<ChainSet> sets, List<SampleCollector> collectors,
                              long intervalSeconds, File logFile) {
        this.sets = sets;
        this.collectors = collectors;
        this.intervalSeconds = intervalSeconds;
        this.logFile = logFile;
        this.burnin = DEFAULT_BURNIN;
    }

    public void start() {
        if (logFile != null) {
            try {
                log = new PrintWriter(new FileWriter(logFile));
                log.println("state\tminESS\tminESS.parameter\tmaxRhat\tmaxRhat.parameter\tASDSF\tclades");
                log.flush();
            } catch (IOException e) {
                throw new RuntimeException("Unable to write the convergence log, " + logFile + ": " + e.getMessage());
            }
        }

        if (intervalSeconds > 0) {
            scheduler = Executors.newSingleThreadScheduledExecutor(r -> {
                Thread thread = new Thread(r, "convergence-monitor");
                thread.setDaemon(true);
                return thread;
            });
            scheduler.scheduleAtFixedRate(() -> {
                try {
                    report("Convergence");
                } catch (RuntimeException e) {
                    // a failed report must not stop the run
                    System.err.println("Convergence report failed: " + e.getMessage());
                }
            }, intervalSeconds, intervalSeconds, TimeUnit.SECONDS);
        }
    }

    /**
     * Stops the periodic reports and writes the final one.
     */
    public void stop() {
        if (scheduler != null) {
            scheduler.shutdownNow();
            try {
                scheduler.awaitTermination(10, TimeUnit.SECONDS);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }
        }
        report("Final convergence diagnostics");
        if (log != null) {
            log.close();
        }
    }

    private synchronized void report(String title) {
        int chainCount = collectors.size();
        SampleCollector.Snapshot[] snapshots = new SampleCollector.Snapshot[chainCount];
        long minState = Long.MAX_VALUE;
        for (int k = 0; k < chainCount; k++) {
            snapshots[k] = collectors.get(k).snapshot();
            minState = Math.min(minState, snapshots[k].state);
        }

        Diagnostics diagnostics = computeDiagnostics(collectors.get(0).getLabels(), snapshots, burnin);

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%n%s at state %d (%d chains, %.0f%% burn-in):%n",
                title, minState, chainCount, burnin * 100));
        for (int k = 0; k < chainCount; k++) {
            String status = sets.get(k).getStatus();
            sb.append(String.format("  %-10s state %-12d posterior %-16.4f %s%n",
                    collectors.get(k).getName(), snapshots[k].state, snapshots[k].posterior, status));
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

        // print in one go so the lines are not split by other output
        System.out.print(sb);
        System.out.flush();

        if (log != null) {
            log.println(minState + "\t" + format(diagnostics.minESS) + "\t" + format(diagnostics.minESSLabel) +
                    "\t" + format(diagnostics.maxRhat) + "\t" + format(diagnostics.maxRhatLabel) +
                    "\t" + format(diagnostics.asdsf) + "\t" + diagnostics.cladeCount);
            log.flush();
        }
    }

    private static String format(double value) {
        return Double.isNaN(value) ? "NA" : Double.toString(value);
    }

    private static String format(String label) {
        return label == null ? "NA" : label;
    }

    static class Diagnostics {
        double minESS = Double.NaN;
        String minESSLabel = null;
        double maxRhat = Double.NaN;
        String maxRhatLabel = null;
        double asdsf = Double.NaN;
        int cladeCount = 0;
    }

    static Diagnostics computeDiagnostics(String[] labels, SampleCollector.Snapshot[] snapshots, double burnin) {
        Diagnostics diagnostics = new Diagnostics();
        int chainCount = snapshots.length;

        // the same number of samples from each chain, after the burn-in
        int sampleCount = Integer.MAX_VALUE;
        for (SampleCollector.Snapshot snapshot : snapshots) {
            sampleCount = Math.min(sampleCount, snapshot.rows.size());
        }
        int start = (int) Math.ceil(sampleCount * burnin);
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

                if (chainCount > 1) {
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

    private static double getESS(double[] values, long stepSize) {
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
    static double splitRhat(double[][] chains) {
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
            while (skip < snapshot.treeBlockSizes.size() && skipped + snapshot.treeBlockSizes.get(skip) <= total * burnin) {
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
    private final long intervalSeconds;
    private final File logFile;
    private final double burnin;

    private ScheduledExecutorService scheduler = null;
    private PrintWriter log = null;
}
