/*
 * MultiChainScreenLogTest.java
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

package test.dr.app.beast;

import dr.app.beast.ConvergenceMonitor;
import dr.app.beast.MultiChainScreenLog;
import dr.inference.loggers.ChainDiagnosticColumn;
import dr.inference.loggers.ChainDiagnostics;
import dr.inference.loggers.ChainScreenFormatter;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import junit.framework.TestCase;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

/**
 * @author Andrew Rambaut
 */
public class MultiChainScreenLogTest extends TestCase {

    private static final String[] LABELS = {"state", "posterior   ", "rate        "};

    private ByteArrayOutputStream bytes;

    private MultiChainScreenLog newLog(int chainCount, boolean hasDiagnostics, LogColumn[] columns, String[] labels) {
        bytes = new ByteArrayOutputStream();
        MultiChainScreenLog log = new MultiChainScreenLog(chainCount, hasDiagnostics, new PrintStream(bytes));
        for (int k = 0; k < chainCount; k++) {
            ChainScreenFormatter formatter = log.getFormatter(k);
            formatter.setColumns(columns);
            formatter.logLabels(labels);
        }
        return log;
    }

    private String[] lines() {
        String text = bytes.toString();
        return text.isEmpty() ? new String[0] : text.split("\r?\n");
    }

    private static LogColumn[] plainColumns() {
        return new LogColumn[]{column("posterior"), column("rate")};
    }

    private static LogColumn column(String label) {
        NumberColumn column = new NumberColumn(label) {
            @Override
            public double getDoubleValue() {
                return 0.0;
            }
        };
        column.setMinimumWidth(12);
        return column;
    }

    private static String[] row(long state, double posterior, double rate, String performance) {
        return new String[]{Long.toString(state), pad(Double.toString(posterior)), pad(Double.toString(rate)), performance};
    }

    private static String pad(String text) {
        StringBuilder sb = new StringBuilder(text);
        while (sb.length() < 12) {
            sb.append(' ');
        }
        return sb.toString();
    }

    private static ChainDiagnostics diagnostics(final double minESS, final double burninProportion) {
        return new ChainDiagnostics() {
            public double getMinESS() {
                return minESS;
            }

            public String getMinESSLabel() {
                return "kappa";
            }

            public double getMaxRhat() {
                return 1.2345;
            }

            public String getMaxRhatLabel() {
                return "rate";
            }

            public double getASDSF() {
                return Double.NaN;
            }

            public double getBurnin() {
                return 1000;
            }

            public double getBurninProportion() {
                return burninProportion;
            }
        };
    }

    /**
     * A state is shown once every chain has logged it, a row for each chain, in state order
     * however the chains interleave.
     */
    public void testBlocksWaitForEveryChain() {
        MultiChainScreenLog log = newLog(3, false, plainColumns(), LABELS);

        log.getFormatter(0).logValues(row(0, -10.5, 1.0, "-"));
        log.getFormatter(0).logValues(row(100, -9.5, 1.5, "1 hours"));
        log.getFormatter(2).logValues(row(0, -12.5, 3.0, "-"));
        assertEquals("nothing until every chain has reached a state", 0, lines().length);

        log.getFormatter(1).logValues(row(0, -11.5, 2.0, "-"));
        String[] lines = lines();
        assertEquals("the heading and a row for each chain", 4, lines.length);
        assertEquals("state       \tchain\tposterior   \trate", lines[0]);
        assertEquals("0           \t1    \t-10.5       \t1.0         \t-", lines[1]);
        assertEquals("            \t2    \t-11.5       \t2.0         \t-", lines[2]);
        assertEquals("            \t3    \t-12.5       \t3.0         \t-", lines[3]);

        // chain 3 gets ahead, then chain 2 catches up twice
        log.getFormatter(2).logValues(row(100, -8.5, 3.5, "3 hours"));
        log.getFormatter(2).logValues(row(200, -7.5, 3.5, "3 hours"));
        log.getFormatter(0).logValues(row(200, -7.0, 1.5, "1 hours"));
        assertEquals(4, lines().length);
        log.getFormatter(1).logValues(row(100, -9.0, 2.5, "2 hours"));
        assertEquals(7, lines().length);
        log.getFormatter(1).logValues(row(200, -7.2, 2.5, "2 hours"));
        lines = lines();
        assertEquals(10, lines.length);
        assertEquals("100         \t1    \t-9.5        \t1.5         \t1 hours", lines[4]);
        assertEquals("            \t2    \t-9.0        \t2.5         \t2 hours", lines[5]);
        assertEquals("200         \t1    \t-7.0        \t1.5         \t1 hours", lines[7]);
        assertEquals("            \t3    \t-7.5        \t3.5         \t3 hours", lines[9]);
    }

    /**
     * The heading and any other lines are shown once, not once for each chain.
     */
    public void testHeadingsShownOnce() {
        MultiChainScreenLog log = newLog(2, false, plainColumns(), LABELS);
        log.getFormatter(0).logHeading("BEAST\nseed");
        log.getFormatter(1).logHeading("BEAST\nseed");
        log.getFormatter(1).logLine("a line");
        log.getFormatter(0).logLine("a line");
        String[] lines = lines();
        assertEquals(3, lines.length);
        assertEquals("# BEAST", lines[0]);
        assertEquals("# seed", lines[1]);
        assertEquals("a line", lines[2]);
    }

    /**
     * Diagnostic columns are shown once for the run, after the chains' own columns wherever
     * they are in the log, and each of the other columns gets its split R-hat and its
     * effective sample size after the burn-in.
     */
    public void testDiagnostics() {
        ChainDiagnosticColumn minESS = new ChainDiagnosticColumn(ChainDiagnosticColumn.Type.MIN_ESS, true);
        minESS.setMinimumWidth(16);
        LogColumn[] columns = {column("posterior"), minESS, column("rate"),
                new ChainDiagnosticColumn(ChainDiagnosticColumn.Type.ASDSF, false)};
        String[] labels = {"state", "posterior   ", minESS.getLabel(), "rate        ", "ASDSF"};

        MultiChainScreenLog log = newLog(2, true, columns, labels);
        final long[] asked = {-1};
        log.setDiagnosticsSource(state -> {
            asked[0] = state;
            return diagnostics(123.4 + state, 0.5);
        });

        int sampleCount = 40;
        double[][][] values = new double[2][2][sampleCount];
        java.util.Random random = new java.util.Random(1);
        for (int s = 0; s < sampleCount; s++) {
            for (int k = 0; k < 2; k++) {
                values[k][0][s] = Math.round(1000 * (random.nextGaussian() + (s < 20 ? 50 : 0))) / 1000.0;
                // constant: no R-hat or effective sample size
                values[k][1][s] = 2.0;
                log.getFormatter(k).logValues(new String[]{Long.toString(s * 10L),
                        pad(Double.toString(values[k][0][s])), "NaN", pad("2.0"), "NaN", "fast"});
            }
            assertEquals("the diagnostics are for the state shown", s * 10L, asked[0]);
        }

        String[] lines = lines();
        assertEquals(1 + sampleCount * 4, lines.length);
        assertEquals("state       \tchain\tposterior   \trate        \tminESS          \tASDSF", lines[0]);

        // the first state: too few samples for the columns' own diagnostics
        assertEquals("0           \t1    \t" + pad(Double.toString(values[0][0][0])) + "\t2.0         \t123 (kappa)     \t-       \tfast", lines[1]);
        assertEquals("            \t2    \t" + pad(Double.toString(values[1][0][0])) + "\t2.0         \t                \t        \tfast", lines[2]);
        assertEquals("            \tR-hat\t-           \t-", lines[3]);
        assertEquals("            \tESS  \t-           \t-", lines[4]);

        // the last state: the first half (which is 50 higher) is burn-in
        int last = 1 + (sampleCount - 1) * 4;
        assertTrue(lines[last].startsWith("390         \t1    \t"));
        assertTrue(lines[last].contains("\t513 (kappa)     \t-       \tfast"));

        double[][] kept = new double[2][20];
        double ess = 0.0;
        for (int k = 0; k < 2; k++) {
            System.arraycopy(values[k][0], 20, kept[k], 0, 20);
            ess += ConvergenceMonitor.getESS(kept[k], 1);
        }
        String rhat = String.format("%.3f", ConvergenceMonitor.splitRhat(kept));
        assertTrue("the burn-in should be left out: R-hat " + rhat, Double.parseDouble(rhat) < 1.5);
        assertEquals("            \tR-hat\t" + pad(rhat) + "\t-", lines[last + 2]);
        assertEquals("            \tESS  \t" + pad(String.format("%.0f", ess)) + "\t-", lines[last + 3]);
    }

    /**
     * Without a source of diagnostics the chains' rows are still shown.
     */
    public void testNoDiagnosticsYet() {
        MultiChainScreenLog log = newLog(2, true, plainColumns(), LABELS);
        log.getFormatter(0).logValues(row(0, -10.5, 1.0, "-"));
        log.getFormatter(1).logValues(row(0, -11.5, 2.0, "-"));
        assertEquals(3, lines().length);
    }
}
