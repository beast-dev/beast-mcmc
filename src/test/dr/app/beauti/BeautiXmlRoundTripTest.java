/*
 * BeautiXmlRoundTripTest.java
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

package test.dr.app.beauti;

import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.util.BeautiImporter;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;

/**
 * Checks how much of an analysis survives being written out as a BEAST XML and read back in: each case is
 * generated, imported into empty options, and generated again, and the two files are compared. Any
 * difference is a setting the import did not recover, and the failure names the line it first shows up on.
 * <p>
 * Cases are listed here as they begin to round-trip; the ones not yet listed are waiting on later stages of
 * the importer.
 *
 * @author Andrew Rambaut
 */
public class BeautiXmlRoundTripTest {

    /**
     * Generates the case, imports the file, generates again and compares.
     */
    private static void roundTrip(String caseName) throws Exception {
        if (BeautiHeadlessHarness.isSkipped(caseName)) {
            return;
        }

        BeautiOptions original = BeautiHeadlessHarness.getCase(caseName).build();
        File file = BeautiHeadlessHarness.generateToFile(caseName, original);
        try {
            String before = BeautiHeadlessHarness.normalised(BeautiHeadlessHarness.read(file));

            BeautiOptions reopened = BeautiOptions.newInstance(BeautiHeadlessHarness.COMPONENTS);
            new BeautiImporter(null, reopened).importFromFile(file);

            String after = BeautiHeadlessHarness.regenerate(caseName, reopened);
            BeautiHeadlessHarness.assertSameXML(caseName, before, after);
        } finally {
            file.delete();
        }
    }

    // ---- cases that survive a round trip in full --------------------------------------------------

    @Test
    public void apes_jc_strict_constant() throws Exception {
        roundTrip("apes_jc_strict_constant");
    }

    @Test
    public void apes_hky_strict_constant() throws Exception {
        roundTrip("apes_hky_strict_constant");
    }

    @Test
    public void apes_hky_gamma_strict_constant() throws Exception {
        roundTrip("apes_hky_gamma_strict_constant");
    }

    @Test
    public void apes_hky_gamma_invar_strict_constant() throws Exception {
        roundTrip("apes_hky_gamma_invar_strict_constant");
    }

    @Test
    public void apes_gtr_gamma_strict_constant() throws Exception {
        roundTrip("apes_gtr_gamma_strict_constant");
    }

    @Test
    public void apes_hky_codon112_strict_constant() throws Exception {
        roundTrip("apes_hky_codon112_strict_constant");
    }

    @Test
    public void primates2_linked_hky_strict_constant() throws Exception {
        roundTrip("primates2_linked_hky_strict_constant");
    }

    @Test
    public void primates2_unlinked_gtr_hky_strict_constant() throws Exception {
        roundTrip("primates2_unlinked_gtr_hky_strict_constant");
    }

    @Test
    public void h1n1_hky_gamma_strict_constant() throws Exception {
        roundTrip("h1n1_hky_gamma_strict_constant");
    }

    @Test
    public void dengue4_hky_gamma_strict_constant() throws Exception {
        roundTrip("dengue4_hky_gamma_strict_constant");
    }

    @Test
    public void dengue4aa_jtt_strict_constant() throws Exception {
        roundTrip("dengue4aa_jtt_strict_constant");
    }

    // ---- the MCMC settings, and the rule about a document that is not empty -----------------------

    /**
     * The chain length, the log intervals and the file names should come back whatever the models are.
     */
    @Test
    public void mcmcSettingsSurviveARoundTrip() throws Exception {
        BeautiOptions original = BeautiHeadlessHarness.getCase("apes_hky_strict_constant").build();
        original.chainLength = 12345678;
        original.logEvery = 4321;
        original.echoEvery = 8642;
        original.checkpointEvery = 999000;
        original.allowOverwriteLog = true;
        original.operatorAnalysis = true;

        File file = BeautiHeadlessHarness.generateToFile("round_trip_mcmc", original);
        try {
            BeautiOptions reopened = BeautiOptions.newInstance(BeautiHeadlessHarness.COMPONENTS);
            new BeautiImporter(null, reopened).importFromFile(file);

            Assert.assertEquals("chain length", 12345678, reopened.chainLength);
            Assert.assertEquals("log every", 4321, reopened.logEvery);
            Assert.assertEquals("echo every", 8642, reopened.echoEvery);
            Assert.assertEquals("checkpoint every", 999000, reopened.checkpointEvery);
            Assert.assertTrue("overwrite the log", reopened.allowOverwriteLog);
            Assert.assertTrue("operator analysis", reopened.operatorAnalysis);
            Assert.assertEquals("file name stem", "round_trip_mcmc", reopened.fileNameStem);
            Assert.assertEquals("log file name", "round_trip_mcmc.log", reopened.logFileName);
        } finally {
            file.delete();
        }
    }

    /**
     * Importing into a document that already has data must leave its settings alone.
     */
    @Test
    public void settingsAreNotTakenIntoANonEmptyDocument() throws Exception {
        BeautiOptions original = BeautiHeadlessHarness.getCase("apes_hky_strict_constant").build();
        original.chainLength = 777000;
        File file = BeautiHeadlessHarness.generateToFile("round_trip_not_empty", original);
        try {
            BeautiOptions loaded = BeautiHeadlessHarness.load("examples/Data/apes.nex");
            int before = loaded.chainLength;

            BeautiImporter importer = new BeautiImporter(null, loaded);
            importer.importFromFile(file);

            Assert.assertEquals("the chain length of the loaded document was replaced", before, loaded.chainLength);
            Assert.assertFalse("the import should say it left the settings alone",
                    importer.getReport().isComplete());
        } finally {
            file.delete();
        }
    }
}
