/*
 * TreeLimeConstantSitesTest.java
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
import dr.app.treelime.options.MaximizationWarnings;
import dr.app.treelime.options.PartitionData;
import dr.app.treelime.options.PartitionTreeModel;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.types.StartingTreeType;
import dr.evolution.alignment.Alignment;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.File;
import java.io.StringReader;
import java.nio.file.Files;
import java.util.Map;

/**
 * Constant-site padding in TreeLime (plan 8 M4a): a SNP alignment (the variable columns) with counts of the constant
 * columns gives the full alignment's likelihood and optimum on the same starting tree.
 */
public class TreeLimeConstantSitesTest extends TestCase {

    private static final String FULL = "examples/Data/Dengue4.env.nex";

    public TreeLimeConstantSitesTest(String name) {
        super(name);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
        dr.math.MathUtils.setSeed(666);
    }

    public void testSnpAlignmentWithCountsEqualsTheFullAlignment() throws Exception {
        TreeLimeOptions full = TreeLimeHarness.importData(FULL);
        if (full == null) {
            return;
        }
        Alignment alignment = ((PartitionData) full.dataPartitions.get(0)).getAlignment();

        // the variable columns as a FASTA file, and the counts of the constant ones
        double[] counts = new double[4];
        StringBuilder[] sequences = new StringBuilder[alignment.getSequenceCount()];
        for (int i = 0; i < sequences.length; i++) {
            sequences[i] = new StringBuilder();
        }
        for (int site = 0; site < alignment.getSiteCount(); site++) {
            int state = alignment.getState(0, site);
            boolean constant = !alignment.getDataType().isAmbiguousState(state);
            for (int i = 1; i < sequences.length && constant; i++) {
                constant = alignment.getState(i, site) == state;
            }
            if (constant) {
                counts[state]++;
            } else {
                for (int i = 0; i < sequences.length; i++) {
                    sequences[i].append(alignment.getDataType().getChar(alignment.getState(i, site)));
                }
            }
        }
        assertTrue(counts[0] + counts[1] + counts[2] + counts[3] > 100);
        File fasta = File.createTempFile("snps", ".fasta");
        fasta.deleteOnExit();
        StringBuilder text = new StringBuilder();
        for (int i = 0; i < sequences.length; i++) {
            text.append('>').append(alignment.getTaxonId(i)).append('\n').append(sequences[i]).append('\n');
        }
        Files.write(fasta.toPath(), text.toString().getBytes());

        TreeLimeOptions snps = TreeLimeHarness.importData(fasta.getPath());
        // the FASTA has no dates: give both the same ones
        TreeLimeHarness.setTipDatesFromNames(full);
        TreeLimeHarness.setTipDatesFromNames(snps);
        PartitionData snpPartition = (PartitionData) snps.dataPartitions.get(0);
        assertEquals(0, snpPartition.getConstantColumnCount());
        snpPartition.setConstantSiteCounts(counts);
        assertTrue(MaximizationWarnings.get(snps).toString(), MaximizationWarnings.get(snps).isEmpty());

        // the same starting tree for both: the full alignment's UPGMA tree
        String fullXml = TreeLimeHarness.generate(quiet(full));
        Tree start = (Tree) parse(fullXml).get("startingTree").getNativeObject();
        for (PartitionTreeModel tree : snps.getPartitionTreeModels()) {
            tree.setStartingTreeType(StartingTreeType.USER);
            tree.setUserStartingTree(start);
        }
        String snpXml = TreeLimeHarness.generate(quiet(snps));
        // BEAST's own form: the counts inside the partition's <patterns>
        String patterns = snpXml.substring(snpXml.indexOf("<patterns id=\"patterns\""), snpXml.indexOf("</patterns>"));
        assertTrue(patterns, patterns.contains("<constantPatterns>") && patterns.contains("<parameter value=\"" +
                (long) counts[0] + " " + (long) counts[1] + " " + (long) counts[2] + " " + (long) counts[3] + "\"/>"));
        assertFalse(snpXml.contains("mergePatterns"));

        Run fullRun = run(fullXml);
        Run snpRun = run(snpXml);
        assertEquals("the likelihood at the start", fullRun.likelihoodAtStart, snpRun.likelihoodAtStart,
                1E-8 * Math.abs(fullRun.likelihoodAtStart));
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, fullRun.maximizer.getReason());
        assertEquals(MaximizerWrtParameter.Reason.CONVERGED, snpRun.maximizer.getReason());
        assertEquals("the optimum", fullRun.end, snpRun.end, 1E-6 * Math.abs(fullRun.end));
    }

    public void testChecksAndWarnings() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(FULL);
        if (options == null) {
            return;
        }
        PartitionData partition = (PartitionData) options.dataPartitions.get(0);
        assertTrue(partition.getConstantColumnCount() > 0);
        // the full alignment already has constant columns: padding it counts them twice
        partition.setConstantSiteCounts(new double[]{10, 20, 30, 40});
        assertTrue(MaximizationWarnings.get(options).toString(),
                MaximizationWarnings.get(options).get(0).contains("counted twice"));

        partition.getPartitionSubstitutionModel().setCodonHeteroPattern("123");
        assertRefused(options, "codon position");
        partition.getPartitionSubstitutionModel().setCodonHeteroPattern(null);

        TreeLimeOptions aminoAcids = TreeLimeHarness.importData("examples/Data/YFV.aa.fasta");
        ((PartitionData) aminoAcids.dataPartitions.get(0)).setConstantSiteCounts(new double[]{1, 2, 3, 4});
        assertRefused(aminoAcids, "nucleotide");
    }

    public void testParseCounts() {
        assertNull(PartitionData.parseConstantSiteCounts("  "));
        double[] counts = PartitionData.parseConstantSiteCounts("1234,567,890,1011\n"); // snp-sites -C
        assertEquals(1011.0, counts[3], 0.0);
        assertEquals(567.0, PartitionData.parseConstantSiteCounts("1234 567\t890 1011")[1], 0.0);
        for (String bad : new String[]{"1 2 3", "1 2 3 x", "1 2 -3 4", "1 2 3 4 5", "1 2 3.5 4"}) {
            try {
                PartitionData.parseConstantSiteCounts(bad);
                fail(bad);
            } catch (IllegalArgumentException e) {
                // expected
            }
        }
    }

    private static void assertRefused(TreeLimeOptions options, String expected) {
        try {
            new dr.app.treelime.generator.BeastGenerator(options, TreeLimeHarness.COMPONENTS).checkOptions();
            fail("expected a refusal mentioning " + expected);
        } catch (dr.app.treelime.generator.Generator.GeneratorException e) {
            assertTrue(e.getMessage(), e.getMessage().contains(expected));
        }
    }

    private static TreeLimeOptions quiet(TreeLimeOptions options) {
        for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
            if (tree.getStartingTreeType() != StartingTreeType.USER) {
                tree.setStartingTreeType(StartingTreeType.UPGMA);
            }
        }
        options.maximizationOptions.nIterations = 500;
        options.maximizationOptions.printToScreen = false;
        return options;
    }

    private static Map<String, XMLObject> parse(String xml) throws Exception {
        xml = xml.replace("<maximizeWrtParameter id=\"maximizer\"",
                "<maximizeWrtParameter id=\"maximizer\" maximizeOnParse=\"false\"");
        xml = xml.substring(0, xml.indexOf("<!-- The optimum")) + "</beast>\n";
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        return parser.parse(new StringReader(xml), false);
    }

    private static class Run {
        MaximizerWrtParameter maximizer;
        double likelihoodAtStart;
        double end;
    }

    private static Run run(String xml) throws Exception {
        Map<String, XMLObject> objects = parse(xml);
        Run run = new Run();
        run.maximizer = (MaximizerWrtParameter) objects.get("maximizer").getNativeObject();
        run.likelihoodAtStart = ((dr.inference.model.Likelihood) objects.get("likelihood").getNativeObject()).getLogLikelihood();
        run.maximizer.maximize();
        run.end = run.maximizer.getLikelihood().getLogLikelihood();
        return run;
    }
}
