/*
 * BeautiHeadlessHarness.java
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

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.ancestralstates.AncestralStatesComponentFactory;
import dr.app.beauti.components.continuous.ContinuousComponentFactory;
import dr.app.beauti.components.discrete.DiscreteTraitsComponentFactory;
import dr.app.beauti.components.hpm.HierarchicalModelComponentFactory;
import dr.app.beauti.components.linkedparameters.LinkedParameterComponentFactory;
import dr.app.beauti.components.marginalLikelihoodEstimation.MarginalLikelihoodEstimationComponentFactory;
import dr.app.beauti.components.sequenceerror.SequenceErrorModelComponentFactory;
import dr.app.beauti.components.tipdatesampling.TipDateSamplingComponentFactory;
import dr.app.beauti.generator.BeastGenerator;
import dr.app.beauti.options.*;
import dr.app.beauti.types.ClockDistributionType;
import dr.app.beauti.types.ClockType;
import dr.app.beauti.types.OperatorSetType;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.util.BEAUTiImporter;
import dr.evomodel.substmodel.aminoacid.AminoAcidModelType;
import dr.evomodel.substmodel.nucleotide.NucModelType;
import org.junit.Assert;
import org.junit.Assume;

import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Drives BEAUti without its user interface: imports a data file from examples/, applies a set of model
 * choices, generates the BEAST XML and compares it with a stored snapshot. Shared by
 * {@link BeautiXmlSnapshotTest} (the XML a given set of choices should produce) and
 * {@link BeautiSaveOpenTest} (that the same XML survives a save and an open).
 * <p>
 * Regenerate the snapshots after an intended change with
 * <code>ant -Dbeauti.snapshots.update=true junit</code>, and run one case with
 * <code>ant -Dbeauti.snapshots.only=&lt;case name&gt; junit</code>.
 *
 * @author Andrew Rambaut
 */
public final class BeautiHeadlessHarness {

    /**
     * The components BeautiFrame installs (see its constructor).
     */
    public static final ComponentFactory[] COMPONENTS = new ComponentFactory[]{
            AncestralStatesComponentFactory.INSTANCE,
            ContinuousComponentFactory.INSTANCE,
            DiscreteTraitsComponentFactory.INSTANCE,
            LinkedParameterComponentFactory.INSTANCE,
            HierarchicalModelComponentFactory.INSTANCE,
            MarginalLikelihoodEstimationComponentFactory.INSTANCE,
            SequenceErrorModelComponentFactory.INSTANCE,
            TipDateSamplingComponentFactory.INSTANCE
    };

    private static final File ROOT = findRoot();
    private static final File SNAPSHOT_DIR = new File(ROOT, "src/test/dr/app/beauti/snapshots");
    private static final File ACTUAL_DIR = new File(ROOT, "build/beauti-snapshots");

    private static final boolean UPDATE = Boolean.getBoolean("beauti.snapshots.update");
    private static final Set<String> ONLY = only();

    static {
        System.setProperty("java.awt.headless", "true");
    }

    private BeautiHeadlessHarness() {
    }

    /**
     * One set of model choices applied to one data file.
     */
    public interface Settings {
        void apply(BeautiOptions options) throws Exception;
    }

    public static final class Case {
        public final String name;
        public final String dataFile;
        private final Settings settings;

        Case(String name, String dataFile, Settings settings) {
            this.name = name;
            this.dataFile = dataFile;
            this.settings = settings;
        }

        /**
         * @return fresh options with the data imported and the choices applied
         */
        public BeautiOptions build() throws Exception {
            BeautiOptions options = load(dataFile);
            settings.apply(options);
            return options;
        }
    }

    private static final Map<String, Case> CASES = new LinkedHashMap<String, Case>();

    private static void add(String name, String dataFile, Settings settings) {
        CASES.put(name, new Case(name, dataFile, settings));
    }

    public static Case getCase(String name) {
        Case testCase = CASES.get(name);
        if (testCase == null) {
            throw new IllegalArgumentException("Unknown case '" + name + "'; known cases: " + CASES.keySet());
        }
        return testCase;
    }

    /**
     * Imports a data file from the repository's examples directory into fresh options.
     */
    public static BeautiOptions load(String path) throws Exception {
        BeautiOptions options = BeautiOptions.newInstance(COMPONENTS);
        new BEAUTiImporter(null, options).importFromFile(new File(ROOT, path));
        return options;
    }

    /**
     * Takes each taxon's date from its name, as the Tip Dates panel's "Guess Dates" does, and calibrates the
     * tips with it.
     *
     * @param order    which number in the name is the date
     * @param fromLast whether to count the numbers from the end of the name
     */
    public static void enableTipDates(BeautiOptions options, int order, boolean fromLast) throws Exception {
        DateGuesser guesser = options.dateGuesser;
        guesser.guessDates = true;
        guesser.guessType = DateGuesser.GuessType.ORDER;
        guesser.order = order;
        guesser.fromLast = fromLast;
        guesser.guessDates(options.taxonList);

        options.useTipDates = true;
        options.calculateTipHeights();
        for (PartitionTreeModel treeModel : options.getPartitionTreeModels()) {
            treeModel.setTipCalibrations(true);
        }
    }

    /**
     * Gives each partition its own substitution model, as the Partitions panel's "Unlink Subst Models" does.
     */
    public static void unlinkSubstitutionModels(BeautiOptions options) {
        for (AbstractPartitionData partition : options.dataPartitions) {
            PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();
            if (!model.getName().equals(partition.getName())) {
                PartitionSubstitutionModel newModel =
                        new PartitionSubstitutionModel(options, partition.getName(), model);
                partition.setPartitionSubstitutionModel(newModel);
            }
        }
        options.clearDataPartitionCaches();
    }

    public static PartitionSubstitutionModel substitutionModel(BeautiOptions options) {
        return options.getPartitionSubstitutionModels().get(0);
    }

    public static PartitionClockModel clockModel(BeautiOptions options) {
        return options.getPartitionClockModels().get(0);
    }

    public static PartitionTreePrior treePrior(BeautiOptions options) {
        return options.getPartitionTreePriors().get(0);
    }

    /**
     * Generates the BEAST XML, with the output file names BEAUti would derive from the case name.
     *
     * @return the XML with its version-bearing lines normalised
     */
    public static String generate(String caseName, BeautiOptions options) throws Exception {
        // (the importer replaces the stem with the data file's name, so set it afterwards)
        options.fileNameStem = caseName;
        options.updateFileNames(false);

        BeastGenerator generator = new BeastGenerator(options, COMPONENTS);
        generator.checkOptions();

        File file = File.createTempFile(caseName, ".xml");
        file.deleteOnExit();
        generator.generateXML(file);
        String xml = new String(Files.readAllBytes(file.toPath()), StandardCharsets.UTF_8);
        file.delete();
        return normalise(xml);
    }

    /**
     * Prepares the XML for comparison: rewrites the two lines carrying the BEAST version, which differs
     * between a development build and a release build, and replaces each sequence with its length and
     * checksum, which keeps the snapshots small while still showing up any change to the data written.
     */
    private static String normalise(String xml) {
        StringBuilder builder = new StringBuilder();
        boolean inSequence = false;
        for (String line : xml.replace("\r\n", "\n").split("\n", -1)) {
            String trimmed = line.trim();
            if (line.startsWith("<!-- Generated by BEAUTi ")) {
                line = "<!-- Generated by BEAUTi (version normalised)                              -->";
            } else if (line.matches("<beast version=\"[^\"]*\">")) {
                line = "<beast version=\"(version normalised)\">";
            } else if (trimmed.equals("<sequence>")) {
                inSequence = true;
            } else if (trimmed.equals("</sequence>")) {
                inSequence = false;
            } else if (inSequence && !trimmed.startsWith("<")) {
                line = line.substring(0, line.length() - trimmed.length()) +
                        "(sequence of " + trimmed.length() + " states, checksum " +
                        Integer.toHexString(trimmed.hashCode()) + ")";
            }
            builder.append(line).append('\n');
        }
        return builder.toString();
    }

    /**
     * Compares the XML with the case's stored snapshot (or, in update mode, replaces it).
     */
    public static void assertSnapshot(String caseName, String actual) throws IOException {
        File expectedFile = new File(SNAPSHOT_DIR, caseName + ".xml");

        if (UPDATE) {
            if (!SNAPSHOT_DIR.isDirectory() && !SNAPSHOT_DIR.mkdirs()) {
                throw new IOException("Could not create " + SNAPSHOT_DIR);
            }
            Files.write(expectedFile.toPath(), actual.getBytes(StandardCharsets.UTF_8));
            return;
        }

        Assert.assertTrue("No snapshot for '" + caseName + "' at " + expectedFile +
                        "\n  create it with: ant -Dbeauti.snapshots.update=true junit",
                expectedFile.isFile());
        String expected = new String(Files.readAllBytes(expectedFile.toPath()), StandardCharsets.UTF_8)
                .replace("\r\n", "\n");
        if (!expected.equals(actual)) {
            if (ACTUAL_DIR.isDirectory() || ACTUAL_DIR.mkdirs()) {
                Files.write(new File(ACTUAL_DIR, caseName + ".xml").toPath(),
                        actual.getBytes(StandardCharsets.UTF_8));
            }
            Assert.fail("The XML of '" + caseName + "' differs from its snapshot at line " +
                    firstDifference(expected, actual) + "\n  diff " +
                    relative(expectedFile) + " " + relative(new File(ACTUAL_DIR, caseName + ".xml")));
        }
    }

    /**
     * Builds the case, generates its XML and checks it against the snapshot. Cases other than the one named
     * by beauti.snapshots.only are skipped.
     */
    public static void run(String caseName) throws Exception {
        Assume.assumeTrue(ONLY.isEmpty() || ONLY.contains(caseName));
        Case testCase = getCase(caseName);
        assertSnapshot(caseName, generate(caseName, testCase.build()));
    }

    /**
     * Writes the options to bytes and reads them back, as saving and opening a .beauti file does.
     */
    public static BeautiOptions saveAndOpen(BeautiOptions options) throws Exception {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        ObjectOutputStream out = new ObjectOutputStream(bytes);
        try {
            out.writeObject(options);
        } finally {
            out.close();
        }
        ObjectInputStream in = new ObjectInputStream(new ByteArrayInputStream(bytes.toByteArray()));
        try {
            BeautiOptions copy = (BeautiOptions) in.readObject();
            copy.registerComponents(COMPONENTS);
            BeautiOptions.setInstance(copy);
            return copy;
        } finally {
            in.close();
        }
    }

    public static boolean isSkipped(String caseName) {
        return !ONLY.isEmpty() && !ONLY.contains(caseName);
    }

    private static int firstDifference(String expected, String actual) {
        String[] expectedLines = expected.split("\n", -1);
        String[] actualLines = actual.split("\n", -1);
        for (int i = 0; i < Math.min(expectedLines.length, actualLines.length); i++) {
            if (!expectedLines[i].equals(actualLines[i])) {
                return i + 1;
            }
        }
        return Math.min(expectedLines.length, actualLines.length) + 1;
    }

    private static String relative(File file) {
        String root = ROOT.getAbsolutePath() + File.separator;
        String path = file.getAbsolutePath();
        return path.startsWith(root) ? path.substring(root.length()) : path;
    }

    /**
     * @return the repository root, found by walking up from the working directory
     */
    private static File findRoot() {
        File directory = new File(System.getProperty("user.dir")).getAbsoluteFile();
        while (directory != null) {
            if (new File(directory, "examples/Data/apes.nex").isFile()) {
                return directory;
            }
            directory = directory.getParentFile();
        }
        throw new IllegalStateException("Could not find the repository root (with examples/Data) from " +
                System.getProperty("user.dir"));
    }

    private static Set<String> only() {
        String property = System.getProperty("beauti.snapshots.only", "").trim();
        Set<String> names = new HashSet<String>();
        if (!property.isEmpty()) {
            names.addAll(Arrays.asList(property.split("\\s*,\\s*")));
        }
        return names;
    }

    /**
     * @return the names of every case, in the order they were defined
     */
    public static List<String> caseNames() {
        return new ArrayList<String>(CASES.keySet());
    }

    static {
        // ---- apes.nex (14 taxa, no dates): the substitution and site models ----
        add("apes_jc_strict_constant", "examples/Data/apes.nex", options ->
                substitutionModel(options).setNucSubstitutionModel(NucModelType.JC));
        add("apes_hky_strict_constant", "examples/Data/apes.nex", options -> {
        });
        add("apes_hky_gamma_strict_constant", "examples/Data/apes.nex", options ->
                substitutionModel(options).setGammaHetero(true));
        add("apes_hky_gamma_invar_strict_constant", "examples/Data/apes.nex", options -> {
            substitutionModel(options).setGammaHetero(true);
            substitutionModel(options).setInvarHetero(true);
        });
        add("apes_gtr_gamma_strict_constant", "examples/Data/apes.nex", options -> {
            substitutionModel(options).setNucSubstitutionModel(NucModelType.GTR);
            substitutionModel(options).setGammaHetero(true);
        });
        add("apes_gtr_gamma_strict_constant_classic_priors", "examples/Data/apes.nex", options -> {
            substitutionModel(options).setNucSubstitutionModel(NucModelType.GTR);
            substitutionModel(options).setGammaHetero(true);
            options.useClassicOperatorsAndPriors = true;
        });
        add("apes_hky_codon112_strict_constant", "examples/Data/apes.nex", options ->
                substitutionModel(options).setCodonHeteroPattern("112"));
        add("apes_hky_codon123_gamma_strict_yule", "examples/Data/apes.nex", options -> {
            substitutionModel(options).setCodonHeteroPattern("123");
            substitutionModel(options).setGammaHetero(true);
            treePrior(options).setNodeHeightPrior(TreePriorType.YULE);
        });

        // ---- apes.nex: the tree priors, clocks and operator mixes ----
        add("apes_hky_strict_yule", "examples/Data/apes.nex", options ->
                treePrior(options).setNodeHeightPrior(TreePriorType.YULE));
        add("apes_hky_strict_birthdeath", "examples/Data/apes.nex", options ->
                treePrior(options).setNodeHeightPrior(TreePriorType.BIRTH_DEATH));
        add("apes_hky_strict_exponential", "examples/Data/apes.nex", options ->
                treePrior(options).setNodeHeightPrior(TreePriorType.EXPONENTIAL));
        add("apes_hky_ucld_constant", "examples/Data/apes.nex", options ->
                clockModel(options).setClockType(ClockType.UNCORRELATED));
        add("apes_hky_uced_constant", "examples/Data/apes.nex", options -> {
            clockModel(options).setClockType(ClockType.UNCORRELATED);
            clockModel(options).setClockDistributionType(ClockDistributionType.EXPONENTIAL);
        });
        add("apes_hky_strict_constant_classic_operators", "examples/Data/apes.nex", options ->
                options.operatorSetType = OperatorSetType.CLASSIC);
        add("apes_hky_strict_constant_fixed_topology", "examples/Data/apes.nex", options ->
                options.operatorSetType = OperatorSetType.FIXED_TREE_TOPOLOGY);

        // ---- Primates_two_partition.nex (12 taxa, two charsets): linked and unlinked partitions ----
        add("primates2_linked_hky_strict_constant", "examples/Data/Primates_two_partition.nex", options -> {
        });
        add("primates2_unlinked_gtr_hky_strict_constant", "examples/Data/Primates_two_partition.nex", options -> {
            unlinkSubstitutionModels(options);
            substitutionModel(options).setNucSubstitutionModel(NucModelType.GTR);
        });
        add("primates2_unlinked_codon112_ucld_yule", "examples/Data/Primates_two_partition.nex", options -> {
            unlinkSubstitutionModels(options);
            for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
                model.setCodonHeteroPattern("112");
            }
            clockModel(options).setClockType(ClockType.UNCORRELATED);
            treePrior(options).setNodeHeightPrior(TreePriorType.YULE);
        });

        // ---- H1N1_HA.nex (two charsets, dates in the taxon names) ----
        add("h1n1_hky_gamma_strict_constant", "examples/Data/H1N1_HA.nex", options ->
                substitutionModel(options).setGammaHetero(true));
        add("h1n1_hky_gamma_strict_constant_tipdates", "examples/Data/H1N1_HA.nex", options -> {
            substitutionModel(options).setGammaHetero(true);
            enableTipDates(options, 0, true);
        });
        add("h1n1_hky_gamma_ucld_skygrid_tipdates", "examples/Data/H1N1_HA.nex", options -> {
            substitutionModel(options).setGammaHetero(true);
            enableTipDates(options, 0, true);
            clockModel(options).setClockType(ClockType.UNCORRELATED);
            treePrior(options).setNodeHeightPrior(TreePriorType.SKYGRID);
            treePrior(options).setSkyGridCount(50);
            treePrior(options).setSkyGridInterval(100.0);
        });

        // ---- Dengue4.env.nex (17 taxa, two-digit years in the taxon names) ----
        add("dengue4_hky_gamma_strict_constant", "examples/Data/Dengue4.env.nex", options ->
                substitutionModel(options).setGammaHetero(true));
        add("dengue4_hky_gamma_strict_constant_tipdates", "examples/Data/Dengue4.env.nex", options -> {
            substitutionModel(options).setGammaHetero(true);
            enableTipDates(options, 0, true);
        });
        add("dengue4_gtr_gamma_ucld_exponential_tipdates", "examples/Data/Dengue4.env.nex", options -> {
            substitutionModel(options).setNucSubstitutionModel(NucModelType.GTR);
            substitutionModel(options).setGammaHetero(true);
            enableTipDates(options, 0, true);
            clockModel(options).setClockType(ClockType.UNCORRELATED);
            treePrior(options).setNodeHeightPrior(TreePriorType.EXPONENTIAL);
        });
        add("dengue4_hky_gamma_uced_skyride_tipdates", "examples/Data/Dengue4.env.nex", options -> {
            substitutionModel(options).setGammaHetero(true);
            enableTipDates(options, 0, true);
            clockModel(options).setClockType(ClockType.UNCORRELATED);
            clockModel(options).setClockDistributionType(ClockDistributionType.EXPONENTIAL);
            treePrior(options).setNodeHeightPrior(TreePriorType.GMRF_SKYRIDE);
        });
        add("dengue4_hky_gamma_strict_skygrid_tipdates", "examples/Data/Dengue4.env.nex", options -> {
            substitutionModel(options).setGammaHetero(true);
            enableTipDates(options, 0, true);
            treePrior(options).setNodeHeightPrior(TreePriorType.SKYGRID);
            treePrior(options).setSkyGridCount(50);
            treePrior(options).setSkyGridInterval(50.0);
        });
        add("dengue4_hky_gamma_strict_constant_tipdates_classic", "examples/Data/Dengue4.env.nex", options -> {
            substitutionModel(options).setGammaHetero(true);
            enableTipDates(options, 0, true);
            options.operatorSetType = OperatorSetType.CLASSIC;
        });

        // ---- YFV.nex (71 taxa, four-digit years at the start of the taxon names) ----
        add("yfv_hky_gamma_strict_constant_tipdates", "examples/Data/YFV.nex", options -> {
            substitutionModel(options).setGammaHetero(true);
            enableTipDates(options, 0, false);
        });
        add("yfv_hky_gamma_ucld_skygrid_tipdates", "examples/Data/YFV.nex", options -> {
            substitutionModel(options).setGammaHetero(true);
            enableTipDates(options, 0, false);
            clockModel(options).setClockType(ClockType.UNCORRELATED);
            treePrior(options).setNodeHeightPrior(TreePriorType.SKYGRID);
            treePrior(options).setSkyGridCount(50);
            treePrior(options).setSkyGridInterval(100.0);
        });
        add("yfv_hky_gamma_invar_strict_yule", "examples/Data/YFV.nex", options -> {
            substitutionModel(options).setGammaHetero(true);
            substitutionModel(options).setInvarHetero(true);
            treePrior(options).setNodeHeightPrior(TreePriorType.YULE);
        });

        // ---- FASTA files, nucleotide and amino acid ----
        add("carnivores_hky_gamma_strict_yule", "examples/Data/carnivores.fasta", options -> {
            substitutionModel(options).setGammaHetero(true);
            treePrior(options).setNodeHeightPrior(TreePriorType.YULE);
        });
        add("carnivores_gtr_codon123_ucld_birthdeath", "examples/Data/carnivores.fasta", options -> {
            substitutionModel(options).setNucSubstitutionModel(NucModelType.GTR);
            substitutionModel(options).setCodonHeteroPattern("123");
            clockModel(options).setClockType(ClockType.UNCORRELATED);
            treePrior(options).setNodeHeightPrior(TreePriorType.BIRTH_DEATH);
        });
        add("dengue4aa_jtt_strict_constant", "examples/Data/Dengue4.env.aa.fasta", options ->
                substitutionModel(options).setAaSubstitutionModel(AminoAcidModelType.JTT));
        add("dengue4aa_blosum62_gamma_ucld_yule", "examples/Data/Dengue4.env.aa.fasta", options -> {
            substitutionModel(options).setGammaHetero(true);
            clockModel(options).setClockType(ClockType.UNCORRELATED);
            treePrior(options).setNodeHeightPrior(TreePriorType.YULE);
        });
    }
}
