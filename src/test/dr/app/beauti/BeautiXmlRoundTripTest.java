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
import dr.app.beauti.options.Operator;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.types.OperatorSetType;
import dr.app.beauti.types.PriorType;
import dr.app.beauti.util.BeautiImporter;
import dr.app.beauti.util.ImportReport;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.FileWriter;

/**
 * Checks how much of an analysis survives being written out as a BEAST XML and read back in: each case is
 * generated, imported into empty options, and generated again, and the two files are compared. Any
 * difference is a setting the import did not recover, and the failure names the line it first shows up on.
 * <p>
 * Every one of the harness's cases is listed here: all of them come back unchanged.
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

    // ---- every case the harness builds --------------------------------------------------

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
    public void apes_gtr_gamma_strict_constant_classic_priors() throws Exception {
        roundTrip("apes_gtr_gamma_strict_constant_classic_priors");
    }

    @Test
    public void apes_hky_codon112_strict_constant() throws Exception {
        roundTrip("apes_hky_codon112_strict_constant");
    }

    @Test
    public void apes_hky_codon123_gamma_strict_yule() throws Exception {
        roundTrip("apes_hky_codon123_gamma_strict_yule");
    }

    @Test
    public void apes_hky_strict_yule() throws Exception {
        roundTrip("apes_hky_strict_yule");
    }

    @Test
    public void apes_hky_strict_birthdeath() throws Exception {
        roundTrip("apes_hky_strict_birthdeath");
    }

    @Test
    public void apes_hky_strict_exponential() throws Exception {
        roundTrip("apes_hky_strict_exponential");
    }

    @Test
    public void apes_hky_ucld_constant() throws Exception {
        roundTrip("apes_hky_ucld_constant");
    }

    @Test
    public void apes_hky_uced_constant() throws Exception {
        roundTrip("apes_hky_uced_constant");
    }

    @Test
    public void apes_hky_strict_constant_classic_operators() throws Exception {
        roundTrip("apes_hky_strict_constant_classic_operators");
    }

    @Test
    public void apes_hky_strict_constant_fixed_topology() throws Exception {
        roundTrip("apes_hky_strict_constant_fixed_topology");
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
    public void primates2_unlinked_codon112_ucld_yule() throws Exception {
        roundTrip("primates2_unlinked_codon112_ucld_yule");
    }

    @Test
    public void h1n1_hky_gamma_strict_constant() throws Exception {
        roundTrip("h1n1_hky_gamma_strict_constant");
    }

    @Test
    public void h1n1_hky_gamma_strict_constant_tipdates() throws Exception {
        roundTrip("h1n1_hky_gamma_strict_constant_tipdates");
    }

    @Test
    public void h1n1_hky_gamma_ucld_skygrid_tipdates() throws Exception {
        roundTrip("h1n1_hky_gamma_ucld_skygrid_tipdates");
    }

    @Test
    public void dengue4_hky_gamma_strict_constant() throws Exception {
        roundTrip("dengue4_hky_gamma_strict_constant");
    }

    @Test
    public void dengue4_hky_gamma_strict_constant_tipdates() throws Exception {
        roundTrip("dengue4_hky_gamma_strict_constant_tipdates");
    }

    @Test
    public void dengue4_gtr_gamma_ucld_exponential_tipdates() throws Exception {
        roundTrip("dengue4_gtr_gamma_ucld_exponential_tipdates");
    }

    @Test
    public void dengue4_hky_gamma_uced_skyride_tipdates() throws Exception {
        roundTrip("dengue4_hky_gamma_uced_skyride_tipdates");
    }

    @Test
    public void dengue4_hky_gamma_strict_skygrid_tipdates() throws Exception {
        roundTrip("dengue4_hky_gamma_strict_skygrid_tipdates");
    }

    @Test
    public void dengue4_hky_gamma_strict_constant_tipdates_classic() throws Exception {
        roundTrip("dengue4_hky_gamma_strict_constant_tipdates_classic");
    }

    @Test
    public void yfv_hky_gamma_strict_constant_tipdates() throws Exception {
        roundTrip("yfv_hky_gamma_strict_constant_tipdates");
    }

    @Test
    public void yfv_hky_gamma_ucld_skygrid_tipdates() throws Exception {
        roundTrip("yfv_hky_gamma_ucld_skygrid_tipdates");
    }

    @Test
    public void yfv_hky_gamma_invar_strict_yule() throws Exception {
        roundTrip("yfv_hky_gamma_invar_strict_yule");
    }

    @Test
    public void carnivores_hky_gamma_strict_yule() throws Exception {
        roundTrip("carnivores_hky_gamma_strict_yule");
    }

    @Test
    public void carnivores_gtr_codon123_ucld_birthdeath() throws Exception {
        roundTrip("carnivores_gtr_codon123_ucld_birthdeath");
    }

    @Test
    public void dengue4aa_jtt_strict_constant() throws Exception {
        roundTrip("dengue4aa_jtt_strict_constant");
    }

    @Test
    public void dengue4aa_blosum62_gamma_ucld_yule() throws Exception {
        roundTrip("dengue4aa_blosum62_gamma_ucld_yule");
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
     * Settings the user changed by hand, which no default would reproduce: a prior of another kind with its
     * own numbers and a truncation, a starting value, and the weight and tuning of two operators.
     */
    @Test
    public void editedPriorsValuesAndOperatorsSurviveARoundTrip() throws Exception {
        BeautiOptions original = BeautiHeadlessHarness.getCase("apes_hky_strict_constant").build();

        Parameter kappa = parameter(original, "kappa");
        kappa.priorType = PriorType.NORMAL_PRIOR;
        kappa.mean = 5.0;
        kappa.stdev = 2.5;
        kappa.isTruncated = true;
        kappa.truncationLower = 0.1;
        kappa.truncationUpper = 100.0;
        kappa.setPriorEdited(true);

        Parameter popSize = parameter(original, "constant.popSize");
        popSize.setInitial(17.5);
        popSize.setInitialSet(true);

        operator(original, "kappa").setWeight(7.0);
        operator(original, "kappa").setTuning(0.33);
        operator(original, "frequencies").setWeight(11.0);
        operator(original, "frequencies").setTuning(0.07);

        File file = BeautiHeadlessHarness.generateToFile("round_trip_edited", original);
        try {
            String before = BeautiHeadlessHarness.normalised(BeautiHeadlessHarness.read(file));

            BeautiOptions reopened = BeautiOptions.newInstance(BeautiHeadlessHarness.COMPONENTS);
            new BeautiImporter(null, reopened).importFromFile(file);

            Parameter read = parameter(reopened, "kappa");
            Assert.assertEquals("the kind of prior", PriorType.NORMAL_PRIOR, read.priorType);
            Assert.assertEquals("the mean of the prior", 5.0, read.mean, 1e-12);
            Assert.assertEquals("the standard deviation of the prior", 2.5, read.stdev, 1e-12);
            Assert.assertTrue("the prior should still be truncated", read.isTruncated);
            Assert.assertEquals("the lower truncation", 0.1, read.truncationLower, 1e-12);
            Assert.assertEquals("the upper truncation", 100.0, read.truncationUpper, 1e-12);

            Assert.assertEquals("the starting value", 17.5,
                    parameter(reopened, "constant.popSize").getInitial(), 1e-12);

            Assert.assertEquals("the weight of the kappa operator", 7.0,
                    operator(reopened, "kappa").getWeight(), 1e-12);
            Assert.assertEquals("the tuning of the kappa operator", 0.33,
                    operator(reopened, "kappa").getTuning(), 1e-12);
            Assert.assertEquals("the weight of the frequencies operator", 11.0,
                    operator(reopened, "frequencies").getWeight(), 1e-12);
            Assert.assertEquals("the tuning of the frequencies operator", 0.07,
                    operator(reopened, "frequencies").getTuning(), 1e-12);

            String after = BeautiHeadlessHarness.regenerate("round_trip_edited", reopened);
            BeautiHeadlessHarness.assertSameXML("round_trip_edited", before, after);
        } finally {
            file.delete();
        }
    }

    private static Parameter parameter(BeautiOptions options, String name) {
        for (Parameter parameter : options.selectParameters()) {
            if (parameter.getName().equals(name)) {
                return parameter;
            }
        }
        throw new IllegalArgumentException("no parameter called " + name);
    }

    private static Operator operator(BeautiOptions options, String name) {
        for (Operator operator : options.selectOperators()) {
            if (operator.getName().equals(name) && operator.isUsed()) {
                return operator;
            }
        }
        throw new IllegalArgumentException("no operator called " + name);
    }

    /**
     * A file BEAUti did not write, carrying a model it cannot represent, must still import its data and say
     * plainly which settings it could not take.
     */
    @Test
    public void anUnrecognisedModelIsReportedRatherThanFailing() throws Exception {
        File file = File.createTempFile("round_trip_foreign", ".xml");
        FileWriter writer = new FileWriter(file);
        writer.write("<?xml version='1.0'?>\n" +
                "<beast>\n" +
                "  <taxa id='taxa'>\n" +
                "    <taxon id='A'/><taxon id='B'/><taxon id='C'/><taxon id='D'/>\n" +
                "  </taxa>\n" +
                "  <alignment id='alignment' dataType='nucleotide'>\n" +
                "    <sequence><taxon idref='A'/>ACGTACGTAC</sequence>\n" +
                "    <sequence><taxon idref='B'/>ACGTACGTAG</sequence>\n" +
                "    <sequence><taxon idref='C'/>ACGTACGTAT</sequence>\n" +
                "    <sequence><taxon idref='D'/>ACGTACGTAA</sequence>\n" +
                "  </alignment>\n" +
                "  <patterns id='patterns' from='1'>\n" +
                "    <alignment idref='alignment'/>\n" +
                "  </patterns>\n" +
                "  <covarionHKYModel id='covarion'/>\n" +
                "  <siteModel id='siteModel'>\n" +
                "    <substitutionModel><covarionHKYModel idref='covarion'/></substitutionModel>\n" +
                "  </siteModel>\n" +
                "</beast>\n");
        writer.close();

        try {
            BeautiOptions options = BeautiOptions.newInstance(BeautiHeadlessHarness.COMPONENTS);
            BeautiImporter importer = new BeautiImporter(null, options);
            importer.importFromFile(file);

            Assert.assertEquals("the sequences should still be read", 1, options.dataPartitions.size());
            Assert.assertEquals("the taxa should still be read", 4, options.taxonList.getTaxonCount());

            ImportReport report = importer.getReport();
            Assert.assertFalse("the import should not claim to be complete", report.isComplete());
            Assert.assertTrue("the unrecognised model should be named",
                    containing(report.getNotRecovered(), "covarionHKYModel"));
            Assert.assertTrue("the missing MCMC settings should be named",
                    containing(report.getNotRecovered(), "MCMC"));
        } finally {
            file.delete();
        }
    }

    private static boolean containing(java.util.List<String> lines, String text) {
        for (String line : lines) {
            if (line.contains(text)) {
                return true;
            }
        }
        return false;
    }

    /**
     * One of the benchmark files, written by hand years before the current BEAUti, to check that a file it
     * did not write still gives up most of its settings and says what it could not.
     */
    @Test
    public void aHandWrittenBenchmarkFileGivesUpWhatItCan() throws Exception {
        File file = new File("examples/Benchmarks/benchmark1.xml");
        Assert.assertTrue("the benchmark file should be in the repository", file.exists());

        BeautiOptions options = BeautiOptions.newInstance(BeautiHeadlessHarness.COMPONENTS);
        BeautiImporter importer = new BeautiImporter(null, options);
        importer.importFromFile(file);

        Assert.assertEquals("the sequences should be read", 1, options.dataPartitions.size());
        Assert.assertEquals("the taxa should be read", 1441, options.taxonList.getTaxonCount());
        Assert.assertEquals("the operator set", OperatorSetType.CLASSIC, options.operatorSetType);
        Assert.assertTrue("the older parameterisation", options.useClassicOperatorsAndPriors);
        Assert.assertEquals("the chain length", 100000, options.chainLength);

        ImportReport report = importer.getReport();
        Assert.assertTrue("the operators it could not tell apart should be reported",
                containing(report.getNotRecovered(), "could not be told apart"));
        Assert.assertTrue("the operators it could match should be taken",
                containing(report.getRecovered(), "weights and tunings"));
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
