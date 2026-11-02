/*
 * LogCombinerBurninTest.java
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
 */

package test.dr.app.tools.logcombiner;

import dr.app.tools.logcombiner.LogCombiner;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

/**
 * Combining log and tree files with a burn-in and a lower sampling frequency.
 *
 * Each sample carries its own state number in a data column (and in a tree comment), so that every row of
 * the combined file can be traced back to the sample it was taken from. The samples taken must be evenly
 * spaced at the resampling frequency whatever burn-in is applied, and the state numbers written must run
 * consecutively through the combined output.
 *
 * See <a href="https://github.com/beast-dev/beast-mcmc/issues/857">issue #857</a>, where a burn-in shifted
 * the samples taken by one sampling interval and the samples taken from the second file depended on how
 * many had been taken from the first.
 */
public class LogCombinerBurninTest {

    private static final long STEP = 10000;
    private static final int SAMPLES = 200;

    /**
     * A log file whose samples run from state 0 to 2000000 every 10000. The third column of each row holds
     * the state the row was written for.
     */
    private File log(String name) throws Exception {
        File file = File.createTempFile(name, ".log");
        file.deleteOnExit();
        PrintWriter writer = new PrintWriter(file);
        writer.println("state\tposterior\ttrueState");
        for (int i = 0; i <= SAMPLES; i++) {
            writer.println((i * STEP) + "\t" + (-1000.0 - i) + "\t" + (i * STEP));
        }
        writer.close();
        return file;
    }

    private File trees(String name) throws Exception {
        File file = File.createTempFile(name, ".trees");
        file.deleteOnExit();
        PrintWriter writer = new PrintWriter(file);
        writer.println("#NEXUS");
        writer.println("Begin trees;");
        for (int i = 0; i <= SAMPLES; i++) {
            writer.println("tree STATE_" + (i * STEP) + " [&trueState=" + (i * STEP) + "] = [&R] (A:1.0,B:1.0);");
        }
        writer.println("End;");
        writer.close();
        return file;
    }

    private File output(String name, String suffix) throws Exception {
        File file = File.createTempFile(name, suffix);
        file.deleteOnExit();
        file.delete();
        return file;
    }

    /** the state number written against each row */
    private List<Long> writtenStates(File file) throws Exception {
        return column(file, 0);
    }

    /** the state each row was actually taken from */
    private List<Long> sourceStates(File file) throws Exception {
        return column(file, 2);
    }

    private List<Long> column(File file, int index) throws Exception {
        List<Long> values = new ArrayList<Long>();
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line = reader.readLine(); // the column headings
        line = reader.readLine();
        while (line != null) {
            values.add(Long.parseLong(line.split("\t")[index]));
            line = reader.readLine();
        }
        reader.close();
        return values;
    }

    private List<Long> treeStates(File file, String prefix) throws Exception {
        List<Long> values = new ArrayList<Long>();
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line = reader.readLine();
        while (line != null) {
            int start = line.indexOf(prefix);
            if (start >= 0) {
                start += prefix.length();
                int end = start;
                while (end < line.length() && Character.isDigit(line.charAt(end))) {
                    end++;
                }
                values.add(Long.parseLong(line.substring(start, end)));
            }
            line = reader.readLine();
        }
        reader.close();
        return values;
    }

    private List<Long> expected(long from, long by, int count) {
        List<Long> values = new ArrayList<Long>();
        for (int i = 0; i < count; i++) {
            values.add(from + (i * by));
        }
        return values;
    }

    private void combine(File[] inputs, long[] burnins, long resample, boolean renumber, File out, boolean treeFiles)
            throws Exception {
        String[] names = new String[inputs.length];
        for (int i = 0; i < inputs.length; i++) {
            names[i] = inputs[i].getPath();
        }
        new LogCombiner(burnins, LogCombiner.BurninType.STATES, resample, names, out.getPath(), treeFiles,
                false, false, renumber, false, 1.0);
    }

    @Test
    public void resamplingWithoutABurninTakesEveryNthState() throws Exception {
        File in = log("input");
        File out = output("combined", ".log");
        combine(new File[]{in}, new long[]{0}, 100000, false, out, false);

        assertEquals(expected(0, 100000, 21), sourceStates(out));
        assertEquals(expected(0, 100000, 21), writtenStates(out));
    }

    @Test
    public void aBurninDoesNotShiftTheStatesTaken() throws Exception {
        File in = log("input");
        File out = output("combined", ".log");
        combine(new File[]{in}, new long[]{1000000}, 100000, false, out, false);

        // the samples taken start at the first one kept and are a resampling frequency apart; before the fix
        // for issue #857 they were one sampling interval short, starting at 1090000
        assertEquals(expected(1000000, 100000, 11), sourceStates(out));
        assertEquals(expected(0, 100000, 11), writtenStates(out));
    }

    @Test
    public void aBurninOnOneFileDoesNotShiftTheStatesTakenFromAnother() throws Exception {
        File first = log("first");
        File second = log("second");
        File out = output("combined", ".log");
        combine(new File[]{first, second}, new long[]{0, 1000000}, 100000, false, out, false);

        List<Long> sources = sourceStates(out);
        assertEquals("every sample of the first file", expected(0, 100000, 21), sources.subList(0, 21));
        // before the fix the second file started at 1090000, because the first file had shifted the count
        assertEquals("every sample of the second file", expected(1000000, 100000, 11), sources.subList(21, 32));

        // the states written run consecutively across the two files so the result reads as one chain
        assertEquals(expected(0, 100000, 32), writtenStates(out));
    }

    @Test
    public void statesWrittenRunConsecutivelyWhenNotResampling() throws Exception {
        File first = log("first");
        File second = log("second");
        File out = output("combined", ".log");
        combine(new File[]{first, second}, new long[]{0, 0}, -1, false, out, false);

        assertEquals(expected(0, STEP, 2 * (SAMPLES + 1)), writtenStates(out));
    }

    @Test
    public void renumberingWithoutResamplingNumbersEverySample() throws Exception {
        File in = log("input");
        File out = output("combined", ".log");
        // this used to fail with an ArithmeticException
        combine(new File[]{in}, new long[]{0}, -1, true, out, false);

        assertEquals(expected(0, 1, SAMPLES + 1), writtenStates(out));
        assertEquals(expected(0, STEP, SAMPLES + 1), sourceStates(out));
    }

    @Test
    public void renumberingWithResamplingNumbersTheSamplesTaken() throws Exception {
        File in = log("input");
        File out = output("combined", ".log");
        combine(new File[]{in}, new long[]{1000000}, 100000, true, out, false);

        assertEquals(expected(0, 1, 11), writtenStates(out));
        assertEquals(expected(1000000, 100000, 11), sourceStates(out));
    }

    @Test
    public void treeFilesAreResampledLikeLogFiles() throws Exception {
        File in = trees("input");
        File out = output("combined", ".trees");
        combine(new File[]{in}, new long[]{1000000}, 100000, false, out, true);

        assertEquals(expected(1000000, 100000, 11), treeStates(out, "[&trueState="));
        assertEquals(expected(0, 100000, 11), treeStates(out, "tree STATE_"));
    }
}
