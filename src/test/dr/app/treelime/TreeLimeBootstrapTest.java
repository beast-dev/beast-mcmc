/*
 * TreeLimeBootstrapTest.java
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

package test.dr.app.treelime;

import dr.app.beast.BeastParser;
import dr.app.beast.BeastVersion;
import dr.app.treelime.generator.BootstrapGenerator;
import dr.app.treelime.options.PartitionData;
import dr.app.treelime.options.PartitionTreeModel;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.types.StartingTreeType;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SiteResampler;
import dr.evolution.tree.Tree;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.File;
import java.io.StringReader;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Bootstrap replicate files (plan 8 S3g): seeded and reproducible, columns resampled within each model partition,
 * the identity replicate gives the original likelihood, the replicate's pattern weights are its counts, and a padded
 * replicate's redrawn constant-site counts are the constant columns it drew.
 */
public class TreeLimeBootstrapTest extends TestCase {

    private static final String DATA = "examples/Data/Dengue4.env.nex";

    public TreeLimeBootstrapTest(String name) {
        super(name);
    }

    public void testSeededFilesAreIdentical() throws Exception {
        TreeLimeOptions options = options();
        if (options == null) {
            return;
        }
        File first = directory();
        File second = directory();
        BootstrapGenerator generator = new BootstrapGenerator(options, TreeLimeHarness.COMPONENTS);
        generator.generate(first, 3, 42);
        generator.generate(second, 3, 42);
        for (String name : new String[]{"rep1.xml", "rep2.xml", "rep3.xml", BootstrapGenerator.COUNTS_FILE,
                BootstrapGenerator.SCRIPT_FILE}) {
            assertTrue(name, Arrays.equals(bytes(first, name), bytes(second, name)));
        }
        assertTrue(new File(first, BootstrapGenerator.SCRIPT_FILE).canExecute());

        // one replicate on its own is the same file
        File single = directory();
        generator.writeReplicate(single, 2, 42);
        assertTrue(Arrays.equals(bytes(first, "rep2.xml"), bytes(single, "rep2.xml")));
        // another seed gives other data
        File other = directory();
        generator.writeReplicate(other, 2, 43);
        assertFalse(Arrays.equals(bytes(first, "rep2.xml"), bytes(other, "rep2.xml")));

        String counts = new String(bytes(first, BootstrapGenerator.COUNTS_FILE));
        assertTrue(counts.contains("master_seed=42"));
        assertEquals("one row per replicate and partition", 3, countRows(counts));
        String xml = new String(bytes(first, "rep1.xml"));
        assertTrue(xml.contains("<alignment id=\"bootstrap.alignment\""));
        assertTrue(xml.contains("<alignment idref=\"bootstrap.alignment\"/>"));
        assertTrue(xml.contains("rep1.max.txt"));
    }

    public void testCountsSumAndMoments() throws Exception {
        TreeLimeOptions options = options();
        if (options == null) {
            return;
        }
        BootstrapGenerator generator = new BootstrapGenerator(options, TreeLimeHarness.COMPONENTS);
        for (int b = 1; b <= 5; b++) {
            for (BootstrapGenerator.PartitionDraw draw : generator.draw(b, 7)) {
                int sum = 0;
                for (int count : draw.siteCounts) {
                    sum += count;
                }
                assertEquals(BootstrapGenerator.sites(draw.partition).length, sum);
            }
        }

        // every column is drawn once on average, with variance 1 - 1/n
        final int n = 20;
        final int draws = 20000;
        double[] sum = new double[n];
        double[] sumOfSquares = new double[n];
        for (int r = 0; r < draws; r++) {
            int[] counts = SiteResampler.resample(n, new Random(SiteResampler.replicateSeed(11, r)));
            for (int i = 0; i < n; i++) {
                sum[i] += counts[i];
                sumOfSquares[i] += counts[i] * counts[i];
            }
        }
        for (int i = 0; i < n; i++) {
            double mean = sum[i] / draws;
            double variance = sumOfSquares[i] / draws - mean * mean;
            assertEquals("mean of column " + i, 1.0, mean, 0.03);
            assertEquals("variance of column " + i, 1.0 - 1.0 / n, variance, 0.06);
        }
    }

    public void testIdentityReplicateIsTheOriginal() throws Exception {
        TreeLimeOptions options = options();
        if (options == null) {
            return;
        }
        double original = logL(TreeLimeHarness.generate(options));
        BootstrapGenerator generator = new BootstrapGenerator(options, TreeLimeHarness.COMPONENTS);
        File directory = directory();
        generator.writeReplicate(directory, 1, generator.identity());
        double identity = logL(new String(bytes(directory, "rep1.xml")));
        assertEquals(original, identity, 1E-10 * Math.abs(original));
    }

    public void testPatternWeightsAreTheCounts() throws Exception {
        TreeLimeOptions options = options();
        if (options == null) {
            return;
        }
        BootstrapGenerator generator = new BootstrapGenerator(options, TreeLimeHarness.COMPONENTS);
        List<BootstrapGenerator.PartitionDraw> draws = generator.draw(4, 99);
        File directory = directory();
        generator.writeReplicate(directory, 4, draws);
        PatternList patterns = (PatternList) parse(new String(bytes(directory, "rep4.xml"))).get("patterns")
                .getNativeObject();

        // the weight of each pattern implied by the counts of the original columns
        BootstrapGenerator.PartitionDraw draw = draws.get(0);
        Alignment alignment = draw.partition.getAlignment();
        int[] sites = BootstrapGenerator.sites(draw.partition);
        Map<String, Double> expected = new HashMap<String, Double>();
        for (int i = 0; i < sites.length; i++) {
            if (draw.siteCounts[i] > 0) {
                expected.merge(Arrays.toString(alignment.getSitePattern(sites[i])), (double) draw.siteCounts[i], Double::sum);
            }
        }
        Map<String, Double> actual = new HashMap<String, Double>();
        for (int p = 0; p < patterns.getPatternCount(); p++) {
            actual.merge(Arrays.toString(patterns.getPattern(p)), patterns.getPatternWeight(p), Double::sum);
        }
        assertEquals(expected, actual);
    }

    // a column drawn for codon position k always fills a position-k slot; under 112 positions 1 and 2 are pooled
    public void testCodonPositionsKeepTheirStructure() throws Exception {
        for (String pattern : new String[]{"123", "112"}) {
            TreeLimeOptions options = options();
            if (options == null) {
                return;
            }
            PartitionData partition = (PartitionData) options.dataPartitions.get(0);
            partition.getPartitionSubstitutionModel().setCodonHeteroPattern(pattern);
            BootstrapGenerator generator = new BootstrapGenerator(options, TreeLimeHarness.COMPONENTS);
            Alignment original = partition.getAlignment();
            int[] sites = BootstrapGenerator.sites(partition);
            for (int b = 1; b <= 3; b++) {
                BootstrapGenerator.PartitionDraw draw = generator.draw(b, 5).get(0);
                Alignment replicate = BootstrapGenerator.replicateAlignment(draw);
                assertEquals(sites.length, replicate.getSiteCount());
                for (int slot = 0; slot < replicate.getSiteCount(); slot++) {
                    int[] column = replicate.getSitePattern(slot);
                    boolean found = false;
                    for (int i = 0; i < sites.length && !found; i++) {
                        found = sameGroup(pattern, i % 3, slot % 3) && draw.siteCounts[i] > 0 &&
                                Arrays.equals(column, original.getSitePattern(sites[i]));
                    }
                    assertTrue(pattern + ": slot " + slot + " has a column from another codon position", found);
                }
            }
            // and the XML still reads the positions as before
            String xml = new String(bytes(writeOne(generator, 1, 5), "rep1.xml"));
            assertTrue(xml.contains("every=\"3\""));
            parse(xml);
        }
    }

    // a padded partition's constant sites are drawn with its columns: the replicate's XML, which pads the redrawn
    // constant-site counts, gives the likelihood of its drawn columns with those constant sites written out
    public void testPaddedReplicateIsItsDraw() throws Exception {
        TreeLimeOptions options = options();
        if (options == null) {
            return;
        }
        PartitionData partition = (PartitionData) options.dataPartitions.get(0);
        partition.setConstantSiteCounts(new double[]{30, 20, 25, 15});
        BootstrapGenerator generator = new BootstrapGenerator(options, TreeLimeHarness.COMPONENTS);
        File directory = directory();
        generator.generate(directory, 2, 21);
        String xml = new String(bytes(directory, "rep2.xml"));
        assertTrue(xml.contains("<constantPatterns>"));

        // replicate 2's draw with its constant sites written out as columns, as a FASTA alignment of its own
        BootstrapGenerator.PartitionDraw draw = generator.draw(2, 21).get(0);
        double drawn = 0;
        for (double count : draw.constantCounts) {
            drawn += count;
        }
        int columns = 0;
        for (int count : draw.siteCounts) {
            columns += count;
        }
        assertEquals("columns and constant sites drawn", BootstrapGenerator.sites(partition).length + 90,
                columns + (int) drawn);
        Alignment replicate = BootstrapGenerator.replicateAlignment(draw);
        StringBuilder fasta = new StringBuilder();
        for (int t = 0; t < replicate.getSequenceCount(); t++) {
            fasta.append('>').append(replicate.getTaxonId(t)).append('\n').append(replicate.getAlignedSequenceString(t));
            for (int state = 0; state < draw.constantCounts.length; state++) {
                for (int k = 0; k < (int) draw.constantCounts[state]; k++) {
                    fasta.append(replicate.getDataType().getChar(state));
                }
            }
            fasta.append('\n');
        }
        File fastaFile = new File(directory, "rep2.fasta");
        Files.write(fastaFile.toPath(), fasta.toString().getBytes());
        TreeLimeOptions written = TreeLimeHarness.importData(fastaFile.getPath());
        TreeLimeHarness.setTipDatesFromNames(written);
        Tree start = (Tree) parse(xml).get("startingTree").getNativeObject();
        for (PartitionTreeModel tree : written.getPartitionTreeModels()) {
            tree.setStartingTreeType(StartingTreeType.USER);
            tree.setUserStartingTree(start);
        }
        written.maximizationOptions.printToScreen = false;
        double fromColumns = logL(TreeLimeHarness.generate(written));
        double fromXml = logL(xml);
        assertEquals(fromXml, fromColumns, 1E-8 * Math.abs(fromXml));
    }

    public void testPriorOnlyIsRefused() throws Exception {
        TreeLimeOptions options = options();
        if (options == null) {
            return;
        }
        options.samplePriorOnly = true;
        try {
            new BootstrapGenerator(options, TreeLimeHarness.COMPONENTS).generate(directory(), 2, 1);
            fail("expected a refusal");
        } catch (dr.app.treelime.generator.Generator.GeneratorException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("prior only"));
        }
    }

    private static boolean sameGroup(String pattern, int a, int b) {
        if (pattern.equals("112")) {
            return (a == 2) == (b == 2);
        }
        return a == b;
    }

    private static int countRows(String counts) {
        int rows = 0;
        for (String line : counts.split("\n")) {
            if (!line.startsWith("#") && !line.startsWith("replicate")) {
                rows++;
            }
        }
        return rows;
    }

    private static File writeOne(BootstrapGenerator generator, int replicate, long seed) throws Exception {
        File directory = directory();
        generator.writeReplicate(directory, replicate, seed);
        return directory;
    }

    // the script runs the replicates JOBS at a time (the option's value unless JOBS is set), then gathers every
    // replicate's time tree in one tree log and every replicate's estimates in one table, leaving out a replicate that
    // wrote nothing; BEAST itself is replaced by the outputs it would write
    public void testScriptGathersTheReplicates() throws Exception {
        TreeLimeOptions options = options();
        if (options == null || System.getProperty("os.name").toLowerCase().startsWith("windows")) {
            return;
        }
        File directory = directory();
        options.maximizationOptions.setBootstrapJobs(3);
        options.maximizationOptions.bootstrapThreads = 2;
        BootstrapGenerator generator = new BootstrapGenerator(options, TreeLimeHarness.COMPONENTS);
        generator.generate(directory, 3, 42);
        String script = new String(bytes(directory, BootstrapGenerator.SCRIPT_FILE));
        assertTrue(script.contains("JOBS=\"${JOBS:-3}\""));
        assertTrue(script.contains("xargs -n 1 -P \"$JOBS\""));
        assertTrue(script.contains("THREADS=\"${THREADS:-2}\""));
        for (int b : new int[]{1, 3}) { // replicate 2 failed
            write(directory, "dengue.rep" + b + ".params.tsv", "name\tvalue\ttransformed\n" +
                    "kappa\t" + (10 + b) + "\t0\nclock.rate\t0.00" + b + "\t0\n");
            write(directory, "dengue.rep" + b + ".timetree.nex", "#NEXUS\n\nBegin taxa;\n\tDimensions ntax=3;\n" +
                    "\tTaxlabels\n\t\tA\n\t\tB\n\t\tC\n\t\t;\nEnd;\n\nBegin trees;\n\tTranslate\n\t\t1 A,\n" +
                    "\t\t2 B,\n\t\t3 C\n\t\t;\ntree STATE_0 [&lnP=-1.0] = [&R] ((1:1.0,2:1.0):" + b + ".0,3:" +
                    (1 + b) + ".0);\nEnd;\n");
        }
        ProcessBuilder process = new ProcessBuilder("sh", BootstrapGenerator.SCRIPT_FILE);
        process.directory(directory);
        // the command each replicate's BEAST runs: it records its arguments, and fails for replicate 2
        write(directory, "fake_beast.sh", "echo \"$@\" >> calls.txt; case \"$@\" in *rep2.xml*) exit 1;; esac\n");
        process.environment().put("BEAST", "sh fake_beast.sh");
        process.environment().put("JOBS", "2");
        process.redirectErrorStream(true);
        Process running = process.start();
        String output = new String(readAll(running.getInputStream()));
        assertEquals(0, running.waitFor());
        assertTrue(output, output.contains("replicate 2 failed: see rep2.log"));
        List<String> calls = Files.readAllLines(new File(directory, "calls.txt").toPath());
        java.util.Collections.sort(calls);
        assertEquals(Arrays.asList("-threads 2 -beagle_threads 2 -overwrite -seed 1 rep1.xml",
                "-threads 2 -beagle_threads 2 -overwrite -seed 2 rep2.xml",
                "-threads 2 -beagle_threads 2 -overwrite -seed 3 rep3.xml"), calls);

        List<String> table = Files.readAllLines(new File(directory, generator.getBootstrapParametersFile()).toPath());
        assertEquals(Arrays.asList("replicate\tkappa\tclock.rate", "1\t11\t0.001", "3\t13\t0.003"), table);

        File trees = new File(directory, generator.getBootstrapTreesFile(options.getPartitionTreeModels().get(0)));
        assertEquals("dengue.bootstrap.trees", trees.getName());
        java.io.FileReader reader = new java.io.FileReader(trees);
        dr.evolution.io.NexusImporter importer = new dr.evolution.io.NexusImporter(reader);
        Tree first = importer.importNextTree();
        Tree second = importer.importNextTree();
        assertFalse(importer.hasTree());
        reader.close();
        assertEquals("REPLICATE_1", first.getId());
        assertEquals("REPLICATE_3", second.getId());
        assertEquals(4.0, second.getNodeHeight(second.getRoot()), 1E-12);
        for (File file : directory.listFiles()) {
            file.deleteOnExit();
        }
    }

    private static byte[] readAll(java.io.InputStream stream) throws Exception {
        java.io.ByteArrayOutputStream bytes = new java.io.ByteArrayOutputStream();
        byte[] buffer = new byte[4096];
        for (int n; (n = stream.read(buffer)) > 0; ) {
            bytes.write(buffer, 0, n);
        }
        return bytes.toByteArray();
    }

    private static void write(File directory, String name, String text) throws Exception {
        Files.write(new File(directory, name).toPath(), text.getBytes("UTF-8"));
    }

    private static TreeLimeOptions options() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(DATA);
        if (options == null) {
            return null;
        }
        TreeLimeHarness.setTipDatesFromNames(options);
        for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
            tree.setStartingTreeType(StartingTreeType.UPGMA);
        }
        options.fileNameStem = "dengue";
        options.maximizationOptions.printToScreen = false;
        options.selectParameters();
        options.selectOperators();
        return options;
    }

    private static File directory() throws Exception {
        File directory = Files.createTempDirectory("treelimeBootstrap").toFile();
        directory.deleteOnExit();
        return directory;
    }

    private static byte[] bytes(File directory, String name) throws Exception {
        return Files.readAllBytes(new File(directory, name).toPath());
    }

    private static Map<String, XMLObject> parse(String xml) throws Exception {
        xml = xml.replace("<maximizeWrtParameter id=\"maximizer\"",
                "<maximizeWrtParameter id=\"maximizer\" maximizeOnParse=\"false\"");
        xml = xml.substring(0, xml.indexOf("<!-- The optimum")) + "</beast>\n";
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        return parser.parse(new StringReader(xml), false);
    }

    private static double logL(String xml) throws Exception {
        return ((dr.inference.model.Likelihood) parse(xml).get("likelihood").getNativeObject()).getLogLikelihood();
    }
}
