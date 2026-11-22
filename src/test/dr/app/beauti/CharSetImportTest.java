/*
 * CharSetImportTest.java
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

package test.dr.app.beauti;

import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.util.BeautiImporter;
import dr.evolution.io.Importer.ImportException;
import org.junit.Test;

import java.io.File;
import java.io.PrintWriter;

import static org.junit.Assert.*;

/**
 * Importing a NEXUS file whose character sets do not fit the alignment they describe.
 *
 * Such a file used to fail while the sites were being taken, with an index error that said nothing about
 * which character set was at fault, and in the application the import simply did not happen.
 *
 * See <a href="https://github.com/beast-dev/beast-mcmc/issues/1066">issue #1066</a>.
 */
public class CharSetImportTest {

    private static final int SITES = 60;
    private static final String[] TAXA = {"Alpha", "Bravo", "Charlie", "Delta"};

    /**
     * A small protein alignment of 60 residues, followed by the given assumptions block.
     */
    private File nexus(String assumptions) throws Exception {
        StringBuilder nexus = new StringBuilder();
        nexus.append("#NEXUS\n\nbegin data;\n");
        nexus.append("\tdimensions ntax=").append(TAXA.length).append(" nchar=").append(SITES).append(";\n");
        nexus.append("\tformat datatype=protein missing=? gap=-;\n\tmatrix\n");
        for (int i = 0; i < TAXA.length; i++) {
            StringBuilder sequence = new StringBuilder();
            for (int j = 0; j < SITES; j++) {
                // a different residue per taxon at each site, so the sites are not all identical
                sequence.append("ACDEFGHIKLMNPQRSTVWY".charAt((i + j) % 20));
            }
            nexus.append("\t").append(TAXA[i]).append("\t").append(sequence).append("\n");
        }
        nexus.append("\t;\nend;\n\n");
        if (assumptions != null) {
            nexus.append(assumptions);
        }

        File file = File.createTempFile("charsets", ".nex");
        file.deleteOnExit();
        PrintWriter writer = new PrintWriter(file, "UTF-8");
        writer.print(nexus);
        writer.close();
        return file;
    }

    private BeautiOptions importFile(File file) throws Exception {
        BeautiOptions options = BeautiOptions.newInstance(BeautiHeadlessHarness.COMPONENTS);
        new BeautiImporter(null, options).importFromFile(file);
        return options;
    }

    /** @return the message of the ImportException, failing the test if the import succeeded */
    private String importFailure(String assumptions) throws Exception {
        try {
            importFile(nexus(assumptions));
        } catch (ImportException expected) {
            return expected.getMessage();
        }
        fail("the import should have been refused");
        return null;
    }

    @Test
    public void characterSetsThatFitAreImported() throws Exception {
        BeautiOptions options = importFile(nexus(
                "begin assumptions;\n\tcharset first = 1-30;\n\tcharset second = 31-60;\nend;\n"));

        assertEquals(2, options.dataPartitions.size());
        for (AbstractPartitionData partition : options.dataPartitions) {
            assertEquals(30, partition.getSiteCount());
            assertEquals("amino acid", partition.getDataType().getDescription());
        }
    }

    @Test
    public void aProteinAlignmentWithoutCharacterSetsIsStillOnePartition() throws Exception {
        BeautiOptions options = importFile(nexus(null));

        assertEquals(1, options.dataPartitions.size());
        assertEquals(SITES, options.dataPartitions.get(0).getSiteCount());
    }

    @Test
    public void aCharacterSetRunningPastTheEndIsReported() throws Exception {
        String message = importFailure(
                "begin assumptions;\n\tcharset first = 1-30;\n\tcharset second = 31-90;\nend;\n");

        assertTrue(message, message.contains("second"));
        assertTrue("it should give the range asked for", message.contains("90"));
        assertTrue("and the length of the alignment", message.contains("60"));
    }

    @Test
    public void aCharacterSetSteppingPastTheEndIsReported() throws Exception {
        String message = importFailure("begin assumptions;\n\tcharset every = 1-70\\3;\nend;\n");

        assertTrue(message, message.contains("every"));
        assertTrue(message, message.contains("70"));
    }

    @Test
    public void aCharacterSetStartingBeforeTheFirstSiteIsReported() throws Exception {
        String message = importFailure("begin assumptions;\n\tcharset zero = 0-30;\nend;\n");

        assertTrue(message, message.contains("zero"));
        assertTrue("character sets are numbered from 1", message.contains("1"));
    }

    @Test
    public void aCharacterSetEndingBeforeItStartsIsReported() throws Exception {
        String message = importFailure("begin assumptions;\n\tcharset backwards = 40-20;\nend;\n");

        assertTrue(message, message.contains("backwards"));
        assertTrue(message, message.contains("ends before it starts"));
    }

    @Test
    public void aCharacterSetRunningToTheEndIsAllowed() throws Exception {
        // '.' means the last site, and is given to the importer as a 'to' of zero
        BeautiOptions options = importFile(nexus(
                "begin assumptions;\n\tcharset first = 1-30;\n\tcharset rest = 31-.;\nend;\n"));

        assertEquals(2, options.dataPartitions.size());
        assertEquals(30, options.dataPartitions.get(1).getSiteCount());
    }
}
