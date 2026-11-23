/*
 * MultiChainScreenLog.java
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

import dr.inference.loggers.ChainDiagnosticColumn;
import dr.inference.loggers.ChainDiagnostics;
import dr.inference.loggers.ChainScreenFormatter;
import dr.inference.loggers.LogColumn;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.LongFunction;

/**
 * Shows the screen logs of several independent chains as one. Each chain's screen log is sent
 * here. When every chain has logged a state, that state is printed as a block: a row for each
 * chain under the usual column headings and, if the run has convergence diagnostics, a row with
 * the split R-hat of each column across the chains and a row with its effective sample size
 * (summed over the chains). Any {@link ChainDiagnosticColumn}s in the screen log are for the
 * run as a whole and are shown once, on the first row of the block.
 *
 * So the screen is updated by state, at the screen log's logEvery, and shows the slowest chain's
 * progress.
 *
 * @author Andrew Rambaut
 */
public class MultiChainScreenLog {

    public static final int STATE_WIDTH = 12;
    public static final int CHAIN_WIDTH = 5;
    public static final int DIAGNOSTIC_WIDTH = 8;

    // Once the diagnostics take a noticeable time to compute, they are not recomputed until
    // this many times as long as they last took has passed, so they cannot slow the chains much
    // however long the logs get. Until then the previous ones are shown.
    public static final int MIN_INTERVAL_FACTOR = 20;
    public static final long NOTICEABLE_NANOS = 10000000L;

    /**
     * @param hasDiagnostics true if convergence diagnostics will be available for the run
     */
    public MultiChainScreenLog(int chainCount, boolean hasDiagnostics, PrintStream out) {
        this.chainCount = chainCount;
        this.hasDiagnostics = hasDiagnostics;
        this.out = out;
        this.history = new ArrayList<List<double[]>>();
        for (int k = 0; k < chainCount; k++) {
            history.add(new ArrayList<double[]>());
        }
    }

    /**
     * @param source gives the diagnostics of the run up to a state that every chain has reached
     */
    public synchronized void setDiagnosticsSource(LongFunction<? extends ChainDiagnostics> source) {
        this.source = source;
    }

    /**
     * @return the destination for the screen log of a chain
     */
    public ChainScreenFormatter getFormatter(final int chain) {
        return new ChainScreenFormatter() {
            @Override
            public boolean hasDiagnostics() {
                return hasDiagnostics;
            }

            @Override
            public void setColumns(LogColumn[] columns) {
                MultiChainScreenLog.this.setColumns(columns);
            }

            @Override
            public void startLogging(String title) {
            }

            @Override
            public void logHeading(String heading) {
                // the same for every chain
                if (chain == 0 && heading != null) {
                    StringBuilder sb = new StringBuilder();
                    for (String line : heading.split("[\r\n]")) {
                        sb.append("# ").append(line).append(String.format("%n"));
                    }
                    print(sb.toString());
                }
            }

            @Override
            public void logLine(String line) {
                if (chain == 0) {
                    print(line + String.format("%n"));
                }
            }

            @Override
            public void logLabels(String[] labels) {
                setLabels(labels);
            }

            @Override
            public void logValues(String[] values) {
                addValues(chain, values);
            }

            @Override
            public void stopLogging() {
                out.flush();
            }
        };
    }

    private synchronized void setColumns(LogColumn[] columns) {
        // every chain has the same columns: those of the first are used for the diagnostics
        if (this.columns == null) {
            this.columns = columns;
        }
    }

    private synchronized void setLabels(String[] labels) {
        if (this.labels == null) {
            this.labels = labels;
        }
    }

    private synchronized void print(String text) {
        out.print(text);
        out.flush();
    }

    private synchronized void addValues(int chain, String[] values) {
        long state;
        try {
            state = Long.parseLong(values[0].trim());
        } catch (NumberFormatException e) {
            return;
        }
        String[][] rows = pending.get(state);
        if (rows == null) {
            rows = new String[chainCount][];
            pending.put(state, rows);
        }
        rows[chain] = values;

        // every chain logs the same states in order, so the earliest is the first to be complete
        while (!pending.isEmpty() && isComplete(pending.firstEntry().getValue())) {
            Map.Entry<Long, String[][]> entry = pending.pollFirstEntry();
            print(getBlock(entry.getKey(), entry.getValue()));
        }
    }

    private static boolean isComplete(String[][] rows) {
        for (String[] row : rows) {
            if (row == null) {
                return false;
            }
        }
        return true;
    }

    private boolean isDiagnostic(int column) {
        return columns != null && column < columns.length && columns[column] instanceof ChainDiagnosticColumn;
    }

    private String getBlock(long state, String[][] rows) {
        // the values after the state: the columns, then anything else (the time per million states)
        int valueCount = rows[0].length - 1;
        int columnCount = (labels != null ? Math.min(labels.length - 1, valueCount) : valueCount);

        StringBuilder sb = new StringBuilder();

        if (!headerShown && labels != null) {
            headerShown = true;
            List<String> cells = new ArrayList<String>();
            cells.add(pad("state", STATE_WIDTH));
            cells.add(pad("chain", CHAIN_WIDTH));
            for (int i = 0; i < columnCount; i++) {
                if (!isDiagnostic(i)) {
                    cells.add(labels[i + 1]);
                }
            }
            for (int i = 0; i < columnCount; i++) {
                if (isDiagnostic(i)) {
                    cells.add(pad(labels[i + 1], DIAGNOSTIC_WIDTH));
                }
            }
            appendLine(sb, cells);
        }

        for (int k = 0; k < chainCount; k++) {
            double[] numbers = new double[columnCount];
            for (int i = 0; i < columnCount; i++) {
                numbers[i] = isDiagnostic(i) ? Double.NaN : parse(rows[k][i + 1]);
            }
            history.get(k).add(numbers);
        }

        if (hasDiagnostics && source != null) {
            long now = System.nanoTime();
            if (summary == null || lastDuration < NOTICEABLE_NANOS ||
                    now - lastComputed >= MIN_INTERVAL_FACTOR * lastDuration) {
                summary = new Summary(state, columnCount);
                lastComputed = System.nanoTime();
                lastDuration = lastComputed - now;
            }
        }

        String state1 = pad(Long.toString(state), STATE_WIDTH);
        String blank = pad("", STATE_WIDTH);
        for (int k = 0; k < chainCount; k++) {
            List<String> cells = new ArrayList<String>();
            cells.add(k == 0 ? state1 : blank);
            cells.add(pad(Integer.toString(k + 1), CHAIN_WIDTH));
            for (int i = 0; i < columnCount; i++) {
                if (!isDiagnostic(i)) {
                    cells.add(rows[k][i + 1]);
                }
            }
            for (int i = 0; i < columnCount; i++) {
                if (isDiagnostic(i)) {
                    int width = Math.max(DIAGNOSTIC_WIDTH, labels[i + 1].length());
                    cells.add(pad(k == 0 && summary != null ? summary.diagnostics[i] : "", width));
                }
            }
            for (int i = columnCount; i < valueCount; i++) {
                cells.add(rows[k][i + 1]);
            }
            appendLine(sb, cells);
        }

        if (summary != null) {
            appendSummaryLine(sb, blank, "R-hat", summary.rhats, columnCount);
            appendSummaryLine(sb, blank, "ESS", summary.ess, columnCount);
        }

        return sb.toString();
    }

    private void appendSummaryLine(StringBuilder sb, String blank, String name, String[] values, int columnCount) {
        List<String> cells = new ArrayList<String>();
        cells.add(blank);
        cells.add(pad(name, CHAIN_WIDTH));
        for (int i = 0; i < columnCount; i++) {
            if (!isDiagnostic(i)) {
                cells.add(pad(values[i], labels != null ? labels[i + 1].length() : 0));
            }
        }
        appendLine(sb, cells);
    }

    private static void appendLine(StringBuilder sb, List<String> cells) {
        StringBuilder line = new StringBuilder();
        for (int i = 0; i < cells.size(); i++) {
            if (i > 0) {
                line.append('\t');
            }
            line.append(cells.get(i));
        }
        // no trailing white space
        int end = line.length();
        while (end > 0 && Character.isWhitespace(line.charAt(end - 1))) {
            end--;
        }
        sb.append(line, 0, end).append(String.format("%n"));
    }

    private static String pad(String text, int width) {
        StringBuilder sb = new StringBuilder(text);
        while (sb.length() < width) {
            sb.append(' ');
        }
        return sb.toString();
    }

    private static double parse(String value) {
        try {
            return Double.parseDouble(value.trim());
        } catch (NumberFormatException e) {
            return Double.NaN;
        }
    }

    /**
     * What is shown for the run as a whole: the diagnostic columns, and the split R-hat and
     * effective sample size of each of the other columns, from the values shown on screen after
     * the burn-in the chains share.
     */
    private class Summary {
        final String[] diagnostics;
        final String[] rhats;
        final String[] ess;

        Summary(long state, int columnCount) {
            ChainDiagnostics run = source.apply(state);

            diagnostics = new String[columnCount];
            rhats = new String[columnCount];
            ess = new String[columnCount];

            int sampleCount = history.get(0).size();
            int start = ConvergenceMonitor.getBurninSamples(sampleCount, run.getBurninProportion());
            int n = sampleCount - start;

            for (int i = 0; i < columnCount; i++) {
                rhats[i] = "-";
                ess[i] = "-";
                if (isDiagnostic(i)) {
                    diagnostics[i] = ((ChainDiagnosticColumn) columns[i]).format(run);
                    continue;
                }
                if (n < 4) {
                    continue;
                }
                double[][] values = new double[chainCount][n];
                boolean constant = true;
                boolean numeric = true;
                for (int k = 0; k < chainCount && numeric; k++) {
                    for (int s = 0; s < n; s++) {
                        double value = history.get(k).get(start + s)[i];
                        if (Double.isNaN(value) || Double.isInfinite(value)) {
                            numeric = false;
                            break;
                        }
                        values[k][s] = value;
                        if (value != values[0][0]) {
                            constant = false;
                        }
                    }
                }
                if (!numeric || constant) {
                    continue;
                }
                double total = 0.0;
                for (int k = 0; k < chainCount; k++) {
                    total += ConvergenceMonitor.getESS(values[k], 1);
                }
                if (!Double.isNaN(total)) {
                    ess[i] = String.format("%.0f", total);
                }
                double rhat = ConvergenceMonitor.splitRhat(values);
                if (!Double.isNaN(rhat)) {
                    rhats[i] = String.format("%.3f", rhat);
                }
            }
        }
    }

    private final int chainCount;
    private final boolean hasDiagnostics;
    private final PrintStream out;

    // all guarded by this
    private LongFunction<? extends ChainDiagnostics> source = null;
    private LogColumn[] columns = null;
    private String[] labels = null;
    private boolean headerShown = false;
    private final TreeMap<Long, String[][]> pending = new TreeMap<Long, String[][]>();
    // the values shown for each chain, by state then column (NaN where not a number)
    private final List<List<double[]>> history;
    private Summary summary = null;
    private long lastComputed = 0;
    private long lastDuration = 0;
}
