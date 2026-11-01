/*
 * OutputFileHelperTest.java
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

package test.dr.util;

import dr.util.OutputFileHelper;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import static org.junit.Assert.*;

public class OutputFileHelperTest {

    private File tempDir;

    @Before
    public void setUp() throws IOException {
        OutputFileHelper.clearOutputDirectory();
        System.clearProperty(OutputFileHelper.OUTPUT_DIR_PROPERTY);
        tempDir = Files.createTempDirectory("beast_output_test").toFile();
    }

    @After
    public void tearDown() {
        OutputFileHelper.clearOutputDirectory();
        System.clearProperty(OutputFileHelper.OUTPUT_DIR_PROPERTY);
        deleteRecursively(tempDir);
    }

    @Test
    public void nothingSetGivesNull() {
        assertNull(OutputFileHelper.getExplicitOutputDirectory());
        assertNull(OutputFileHelper.getOutputDirectory());
    }

    @Test
    public void propertyIsUsedWhenNothingSetProgrammatically() {
        System.setProperty(OutputFileHelper.OUTPUT_DIR_PROPERTY, tempDir.getPath());
        assertEquals(tempDir, OutputFileHelper.getExplicitOutputDirectory());
        assertEquals(tempDir, OutputFileHelper.getOutputDirectory());
    }

    @Test
    public void explicitDirectoryBeatsInferredDirectory() {
        File inferred = new File(tempDir, "logs");
        OutputFileHelper.setDefaultOutputDirectory(inferred);
        assertNull(OutputFileHelper.getExplicitOutputDirectory());
        assertEquals(inferred, OutputFileHelper.getOutputDirectory());

        System.setProperty(OutputFileHelper.OUTPUT_DIR_PROPERTY, tempDir.getPath());
        assertEquals(tempDir, OutputFileHelper.getOutputDirectory());
    }

    @Test
    public void firstInferredDirectoryWins() {
        File first = new File(tempDir, "first");
        File second = new File(tempDir, "second");
        OutputFileHelper.setDefaultOutputDirectory(first);
        OutputFileHelper.setDefaultOutputDirectory(second);
        assertEquals(first, OutputFileHelper.getOutputDirectory());
    }

    @Test
    public void resolveRelativeNameIntoDirectory() {
        File resolved = OutputFileHelper.resolveOutputFile("run.log", null, tempDir);
        assertEquals(new File(tempDir, "run.log"), resolved);
    }

    @Test
    public void resolveStripsLeadingDotSlashAndAppliesPrefix() {
        File resolved = OutputFileHelper.resolveOutputFile("./run.log", "pre_", tempDir);
        assertEquals(new File(tempDir, "pre_run.log"), resolved);
    }

    @Test
    public void resolveKeepsRelativeSubdirectory() {
        File resolved = OutputFileHelper.resolveOutputFile("sub/run.log", "pre_", tempDir);
        assertEquals(new File(new File(tempDir, "sub"), "pre_run.log"), resolved);
    }

    @Test
    public void resolveLeavesAbsoluteNameAlone() {
        File absolute = new File(tempDir, "elsewhere.log").getAbsoluteFile();
        File resolved = OutputFileHelper.resolveOutputFile(absolute.getPath(), null, new File(tempDir, "other"));
        assertEquals(absolute, resolved);
    }

    @Test
    public void resolveWithoutDirectoryFallsBackToWorkingDirectory() {
        File resolved = OutputFileHelper.resolveOutputFile("run.log", null, null);
        assertEquals(new File(System.getProperty("user.dir"), "run.log"), resolved);
    }

    @Test
    public void validateAcceptsWritableDirectory() throws IOException {
        OutputFileHelper.validateOutputDirectory(tempDir);
        assertEquals("probe file must be removed", 0, tempDir.list().length);
    }

    @Test
    public void validateRejectsMissingDirectory() {
        try {
            OutputFileHelper.validateOutputDirectory(new File(tempDir, "missing"));
            fail("expected IOException");
        } catch (IOException e) {
            assertTrue(e.getMessage(), e.getMessage().startsWith("Output directory does not exist"));
        }
    }

    @Test
    public void validateRejectsFile() throws IOException {
        File file = new File(tempDir, "a.txt");
        assertTrue(file.createNewFile());
        try {
            OutputFileHelper.validateOutputDirectory(file);
            fail("expected IOException");
        } catch (IOException e) {
            assertTrue(e.getMessage(), e.getMessage().startsWith("Output path is not a directory"));
        }
    }

    @Test
    public void validateRejectsUnwritableDirectory() throws IOException {
        File readOnly = new File(tempDir, "ro");
        assertTrue(readOnly.mkdir());
        if (!readOnly.setWritable(false, false) || readOnly.canWrite()) {
            // running as root or on a file system without POSIX permissions
            return;
        }
        try {
            OutputFileHelper.validateOutputDirectory(readOnly);
            fail("expected IOException");
        } catch (IOException e) {
            assertTrue(e.getMessage(), e.getMessage().startsWith("Output directory is not writable"));
        } finally {
            readOnly.setWritable(true, false);
        }
    }

    private static void deleteRecursively(File file) {
        if (file == null || !file.exists()) {
            return;
        }
        file.setWritable(true, false);
        File[] children = file.listFiles();
        if (children != null) {
            for (File child : children) {
                deleteRecursively(child);
            }
        }
        file.delete();
    }
}
