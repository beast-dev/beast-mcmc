/*
 * ResumedLogFileTest.java
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

package test.dr.inference.loggers;

import dr.inference.loggers.ResumedLogFile;
import junit.framework.TestCase;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

/**
 * Tests the naming of the log files of a resumed chain and the checks made before appending
 * to an existing log or trees file.
 */
public class ResumedLogFileTest extends TestCase {

    private static final String LOG_HEADER = "# BEAST X\n# Generated\nstate\tposterior\tprior\n";
    private static final String TREES_HEADER = "#NEXUS\n\nBegin taxa;\n\tDimensions ntax=2;\n\tTaxlabels\n" +
            "\t\tA\n\t\tB\n\t\t;\nEnd;\n\nBegin trees;\n\tTranslate\n\t\t1 A,\n\t\t2 B\n\t\t;\n";

    private File directory;

    @Override
    protected void setUp() throws IOException {
        directory = Files.createTempDirectory("resumed-log-test").toFile();
    }

    @Override
    protected void tearDown() {
        final File[] files = directory.listFiles();
        if (files != null) {
            for (File file : files) {
                file.delete();
            }
        }
        directory.delete();
    }

    public void testRenamed() {
        assertEquals("run.from3000.log", ResumedLogFile.renamed(new File("run.log"), 3000).getName());
        assertEquals("run.chain1.from3000.trees",
                ResumedLogFile.renamed(new File("run.chain1.trees"), 3000).getName());
        assertEquals("run.from3000.log.zst", ResumedLogFile.renamed(new File("run.log.zst"), 3000).getName());
        assertEquals("run.from0", ResumedLogFile.renamed(new File("run"), 0).getName());
        assertEquals(new File("out", "run.from5.log"), ResumedLogFile.renamed(new File("out", "run.log"), 5));
    }

    public void testExpectedLastState() {
        assertEquals(3000, ResumedLogFile.expectedLastState(3000, 500));
        assertEquals(3000, ResumedLogFile.expectedLastState(3499, 500));
        assertEquals(0, ResumedLogFile.expectedLastState(499, 500));
    }

    public void testLogThatMatches() throws IOException {
        final String content = LOG_HEADER + logLines(0, 3000, 500);
        final File file = write("match.log", content);

        assertEquals(0, ResumedLogFile.prepareForAppend(file, 3000, 500, false));
        assertEquals(content, read(file));

        // a saved state between two logged states carries on from the earlier one
        assertEquals(0, ResumedLogFile.prepareForAppend(file, 3200, 500, false));
        assertEquals(content, read(file));
    }

    public void testLogThatGoesPast() throws IOException {
        final String content = LOG_HEADER + logLines(0, 5000, 500);
        final File file = write("ahead.log", content);

        try {
            ResumedLogFile.prepareForAppend(file, 3000, 500, false);
            fail("A log that goes past the saved state should not be appended to");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage().contains("5000"));
            assertTrue(expected.getMessage().contains("--force_append"));
        }
        assertEquals("A refused file is left as it is", content, read(file));

        assertEquals(4, ResumedLogFile.prepareForAppend(file, 3000, 500, true));
        assertEquals(LOG_HEADER + logLines(0, 3000, 500), read(file));
    }

    public void testLogThatStopsShort() throws IOException {
        final String content = LOG_HEADER + logLines(0, 2000, 500);
        final File file = write("behind.log", content);

        for (boolean force : new boolean[]{false, true}) {
            try {
                ResumedLogFile.prepareForAppend(file, 3000, 500, force);
                fail("A log with states missing should not be appended to");
            } catch (IllegalStateException expected) {
                assertTrue(expected.getMessage().contains("missing"));
            }
            assertEquals(content, read(file));
        }
    }

    public void testLogWithoutTheStateResumedFrom() throws IOException {
        // logged every 1000 before, every 500 now: state 2500 was never logged
        final String content = LOG_HEADER + logLines(0, 4000, 1000);
        final File file = write("other.log", content);

        try {
            ResumedLogFile.prepareForAppend(file, 2500, 500, true);
            fail("A log without the state to carry on from should not be appended to");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage().contains("2500"));
        }
        assertEquals(content, read(file));
    }

    public void testLogWithNoStates() throws IOException {
        final File file = write("empty.log", LOG_HEADER);
        try {
            ResumedLogFile.prepareForAppend(file, 0, 500, true);
            fail("A log with no states should not be appended to");
        } catch (IllegalStateException expected) {
            // as it should be
        }

        final File nothing = write("nothing.log", "");
        try {
            ResumedLogFile.prepareForAppend(nothing, 0, 500, true);
            fail("An empty file should not be appended to");
        } catch (IllegalStateException expected) {
            // as it should be
        }
    }

    public void testIncompleteLastLine() throws IOException {
        final String content = LOG_HEADER + logLines(0, 3000, 500) + "3500\t-12.";
        final File file = write("partial.log", content);

        try {
            ResumedLogFile.prepareForAppend(file, 3000, 500, false);
            fail("A log that ends with an incomplete line should not be appended to");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage().contains("incomplete"));
        }
        assertEquals(content, read(file));

        assertEquals(1, ResumedLogFile.prepareForAppend(file, 3000, 500, true));
        assertEquals(LOG_HEADER + logLines(0, 3000, 500), read(file));
    }

    public void testTreesFile() throws IOException {
        final String trees = TREES_HEADER + treeLines(0, 3000, 1000);

        // a finished file has its closing End; removed, to be written again at the end
        File file = write("finished.trees", trees + "End;\n");
        assertEquals(0, ResumedLogFile.prepareForAppend(file, 3000, 1000, false));
        assertEquals(trees, read(file));

        // a file of a run that was stopped has none
        file = write("stopped.trees", trees);
        assertEquals(0, ResumedLogFile.prepareForAppend(file, 3000, 1000, false));
        assertEquals(trees, read(file));

        file = write("ahead.trees", TREES_HEADER + treeLines(0, 6000, 1000) + "End;\n");
        try {
            ResumedLogFile.prepareForAppend(file, 3000, 1000, false);
            fail("A trees file that goes past the saved state should not be appended to");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage().contains("6000"));
        }
        assertEquals(3, ResumedLogFile.prepareForAppend(file, 3000, 1000, true));
        assertEquals(trees, read(file));
    }

    public void testLongLines() throws IOException {
        // lines longer than the blocks the file is read in
        final StringBuilder padding = new StringBuilder();
        for (int i = 0; i < 200000; i++) {
            padding.append('x');
        }
        final StringBuilder content = new StringBuilder(TREES_HEADER);
        for (int state = 0; state <= 2000; state += 1000) {
            content.append("tree STATE_").append(state).append(" [&lnP=-1.0] = [&R] (1:").append(padding).append(",2:1.0);\n");
        }
        final String kept = content.toString();
        content.append("tree STATE_3000 = [&R] (1:").append(padding).append(",2:1.0);\nEnd;\n");

        final File file = write("long.trees", content.toString());
        assertEquals(1, ResumedLogFile.prepareForAppend(file, 2000, 1000, true));
        assertEquals(kept, read(file));
    }

    private static String logLines(long from, long to, long every) {
        final StringBuilder sb = new StringBuilder();
        for (long state = from; state <= to; state += every) {
            sb.append(state).append("\t-100.5\t-3.25\n");
        }
        return sb.toString();
    }

    private static String treeLines(long from, long to, long every) {
        final StringBuilder sb = new StringBuilder();
        for (long state = from; state <= to; state += every) {
            sb.append("tree STATE_").append(state).append(" [&lnP=-100.5] = [&R] (1:1.0,2:1.0);\n");
        }
        return sb.toString();
    }

    private File write(String name, String content) throws IOException {
        final File file = new File(directory, name);
        Files.write(file.toPath(), content.getBytes(StandardCharsets.ISO_8859_1));
        return file;
    }

    private static String read(File file) throws IOException {
        return new String(Files.readAllBytes(file.toPath()), StandardCharsets.ISO_8859_1);
    }
}
