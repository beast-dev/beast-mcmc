/*
 * TreeLimeSearchGeneratorTest.java
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
import dr.app.treelime.generator.BeastGenerator;
import dr.app.treelime.generator.Generator.GeneratorException;
import dr.app.treelime.options.MaximizationOptions;
import dr.app.treelime.options.Parameter;
import dr.app.treelime.options.PartitionClockModel;
import dr.app.treelime.options.TopologySearchOptions;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.types.ClockType;
import dr.evolution.util.Taxa;
import dr.evomodel.tree.search.TimeTreeSearch;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.File;
import java.io.StringReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

/**
 * The topology search in TreeLime's generator (plan 7 M2, V12): with the search off the XML is exactly what TreeLime
 * wrote before the search existed (golden files, written once from the generator before the change); with it on, the
 * golden files per mode, a strict parse, and the blockers. Rewrite the golden files with
 * -Dtreelime.golden.update=true after a deliberate change to the generator.
 */
public class TreeLimeSearchGeneratorTest extends TestCase {

    static final String DENGUE = "examples/Data/Dengue4.env.nex";
    private static final File GOLDEN = new File("src/test/dr/app/treelime/golden");
    private static final String STEM = "golden";

    public TreeLimeSearchGeneratorTest(String name) {
        super(name);
    }

    // V12: the search off writes today's XML, byte for byte (but the version lines)
    public void testSearchOffIsUnchanged() throws Exception {
        assertGolden("off_strict_dated.xml", dengue(false));
        assertGolden("off_hmc_dated.xml", dengue(true));
        TreeLimeOptions ml = dengue(false);
        if (ml != null) {
            ml.maximizationOptions.objective = dr.app.treelime.options.MaximizationOptions.Objective.ML;
        }
        assertGolden("off_strict_dated_ml.xml", ml);
        assertGolden("off_two_partitions.xml", TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS));
    }

    // the search on: golden files per clock, every id resolves
    public void testSearchOnGolden() throws Exception {
        TreeLimeOptions strict = search(dengue(false));
        TreeLimeOptions hmc = search(dengue(true));
        if (strict == null) {
            return;
        }
        assertGolden("nni_strict_dated.xml", strict);
        assertGolden("nni_hmc_dated.xml", hmc);
        TreeLimeOptions spr = search(dengue(true));
        spr.maximizationOptions.getTopologySearch().mode = TopologySearchOptions.Mode.NNI_SPR;
        assertGolden("nni_spr_hmc_dated.xml", spr);
        assertTrue(generate(spr).contains("<timeTreeSearch id=\"treeSearch\" seed=\"1\" moves=\"nni,spr,root\""));
        // perturbation: nStop and the candidate trees' file, only when on
        assertFalse(generate(spr).contains("nStop="));
        spr.maximizationOptions.getTopologySearch().nStop = 20;
        String perturbed = generate(spr);
        assertTrue(perturbed.contains("nStop=\"20\" candidatesFileName=\"golden.candidates.trees\""));
        TreeLimeHarness.assertIdsResolve(perturbed);
        String xml = generate(hmc);
        TreeLimeHarness.assertIdsResolve(xml);
        assertTrue(xml.contains("<maximizeWrtParameter id=\"maximizer\" maximizeOnParse=\"false\""));
        // the round maximizer re-fits the node heights and the HMC clock's branch rates only
        String round = between(xml, "<maximizeWrtParameter id=\"roundMaximizer\"", "</maximizeWrtParameter>");
        assertTrue(round.contains("<compoundGradient idref=\"searchGradient\"/>"));
        assertTrue(round.contains("<nodeHeightTransform idref=\"nodeHeightTransform\"/>"));
        String gradient = between(xml, "<compoundGradient id=\"searchGradient\">", "</compoundGradient>");
        assertTrue(gradient.contains("idref=\"gradient.branchRatesHMC.rates\""));
        assertTrue(gradient.contains("idref=\"gradient.nodeHeights\""));
        assertEquals(2, gradient.split("idref=").length - 1);
        // both reports include the search
        String reports = xml.substring(xml.indexOf("<!-- The optimum"));
        assertEquals(2, reports.split("<timeTreeSearch idref=\"treeSearch\"/>").length - 1);
        // the search runs before the reports and the export
        assertTrue(xml.indexOf("<timeTreeSearch id=\"treeSearch\"") < xml.indexOf("<!-- The optimum"));
    }

    // the generated file runs end to end: the search, the reports and the export, with a search section in the report
    // and the final F at least the start's; no parser warning names an element of the search
    public void testSearchOnRuns() throws Exception {
        for (boolean relaxed : new boolean[]{false, true}) {
            TreeLimeOptions options = search(dengue(relaxed));
            if (options == null) {
                return;
            }
            options.maximizationOptions.printToScreen = false;
            File directory = Files.createTempDirectory("treeLimeSearch").toFile();
            options.fileNameStem = new File(directory, "run").getAbsolutePath();
            String xml = generateTo(options, new File(directory, "run.xml"));

            List<String> warnings = new ArrayList<String>();
            Handler handler = collect(warnings);
            Map<String, XMLObject> objects;
            try {
                BeastParser parser = new BeastParser(new String[]{}, null, false, true, false, new BeastVersion());
                objects = parser.parse(new StringReader(xml), false);
            } finally {
                Logger.getLogger("dr.xml").removeHandler(handler);
            }
            for (String warning : warnings) {
                for (String id : new String[]{"treeSearch", "roundMaximizer", "searchGradient", "searchRates",
                        "timeTreeSearch"}) {
                    assertFalse(warning, warning.contains(id));
                }
            }

            TimeTreeSearch search = (TimeTreeSearch) objects.get("treeSearch").getNativeObject();
            assertTrue(search.wasExecuted());
            assertTrue(search.getFinalValue() >= search.getStartValue());
            MaximizerWrtParameter maximizer = (MaximizerWrtParameter) objects.get("maximizer").getNativeObject();
            assertTrue(maximizer.wasExecuted());
            String report = new String(Files.readAllBytes(new File(directory, "run.max.txt").toPath()),
                    StandardCharsets.UTF_8);
            assertTrue(report.contains("Topology search (rooted NNI hill-climb)"));
            assertTrue(report.contains("StopReason:"));
            for (String file : new String[]{"run.search.trees", "run.search.log", "run.timetree.nex"}) {
                assertTrue(file, new File(directory, file).exists());
            }
        }
    }

    // run as BEAST runs it, the file has no runnable (the maximizer and the search run while parsing), yet the
    // citations are logged once at the end, as an MCMC run logs them before it starts, the ratio transform's included
    public void testCitationsWithoutARunnable() throws Exception {
        TreeLimeOptions options = search(dengue(false));
        if (options == null) {
            return;
        }
        options.maximizationOptions.printToScreen = false;
        File directory = Files.createTempDirectory("treeLimeCitations").toFile();
        options.fileNameStem = new File(directory, "run").getAbsolutePath();
        String xml = generateTo(options, new File(directory, "run.xml"));

        final List<String> lines = new ArrayList<String>();
        Handler handler = new Handler() {
            public void publish(LogRecord record) {
                lines.add(record.getMessage());
            }

            public void flush() {
            }

            public void close() {
            }
        };
        Logger citations = Logger.getLogger("dr.util.citations");
        citations.addHandler(handler);
        try {
            BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
            parser.parse(new StringReader(xml), true);
        } finally {
            citations.removeHandler(handler);
        }
        String text = String.join("\n", lines);
        assertEquals(text, 1, count(text, "Citations for this analysis"));
        assertTrue(text, text.contains("BEAST primary citation"));
        assertTrue(text, text.contains("HKY nucleotide substitution model"));
        // the node heights are optimized as ratios
        assertTrue(text, text.contains("Using the ratio transformation of internal node heights"));
        assertTrue(text, text.contains("Scalable Bayesian divergence time estimation with ratio transformations"));
        // and kappa's and the frequencies' gradients
        assertTrue(text, text.contains("Using substitution model parameter gradients"));
        assertTrue(text, text.contains("Detecting evolutionary change-points with branch-specific substitution models"));
    }

    private static int count(String text, String part) {
        int n = 0;
        for (int at = text.indexOf(part); at >= 0; at = text.indexOf(part, at + 1)) {
            n++;
        }
        return n;
    }

    // a tree prior on a subtree keeps the subtree a clade: it is a constraint of the search (the Trees tab offers
    // monophyletic taxon sets only, so the objective's monophyly statistic names it too)
    public void testExcludedSubtreeIsAConstraint() throws Exception {
        TreeLimeOptions options = search(dengue(false));
        if (options == null) {
            return;
        }
        Taxa taxa = new Taxa("subtree");
        for (int i = 0; i < 4; i++) {
            taxa.addTaxon(options.taxonList.getTaxon(i));
        }
        options.taxonSets.add(taxa);
        options.taxonSetsMono.put(taxa, true);
        options.taxonSetsIncludeStem.put(taxa, false);
        options.taxonSetsTreeModel.put(taxa, options.getPartitionTreeModels().get(0));
        options.getPartitionTreeModels().get(0).getPartitionTreePrior().setSubtreeTaxonSet(taxa);
        options.getPartitionTreeModels().get(0).getPartitionTreePrior()
                .setSubtreePrior(dr.app.treelime.types.TreePriorType.CONSTANT);
        String xml = generate(options);
        TreeLimeHarness.assertIdsResolve(xml);
        String search = between(xml, "<timeTreeSearch id=\"treeSearch\"", "</timeTreeSearch>");
        assertTrue(search.contains("<constraints>"));
        assertTrue(search.contains("<taxa idref=\"subtree\"/>"));

        // parsed without running: one constraint, from both sources
        xml = xml.replace("<timeTreeSearch id=\"treeSearch\"", "<timeTreeSearch id=\"treeSearch\" searchOnParse=\"false\"");
        xml = xml.substring(0, xml.indexOf("<!-- The optimum")) + "</beast>\n";
        TimeTreeSearch parsed = (TimeTreeSearch) parse(xml, false).get("treeSearch").getNativeObject();
        assertEquals(1, parsed.getConstraints().size());
        assertTrue(parsed.getConstraints().areSatisfied(parsed.getTree()));
    }

    // the blockers: the generator refuses the search, and the options name the reason
    public void testBlockers() throws Exception {
        TreeLimeOptions options = search(dengue(false));
        if (options == null) {
            return;
        }
        assertNull(TopologySearchOptions.getBlocker(options));

        options.maximizationOptions.objective = MaximizationOptions.Objective.PRIOR;
        assertBlocked(options, "prior-only");
        options.maximizationOptions.objective = MaximizationOptions.Objective.MAP;

        options.maximizationOptions.setGradientMode(MaximizationOptions.GradientMode.NUMERIC);
        assertBlocked(options, "analytic");
        options.maximizationOptions.setGradientMode(MaximizationOptions.GradientMode.ANALYTIC_WHERE_AVAILABLE);

        for (Parameter parameter : options.selectParameters()) {
            if (parameter.isNodeHeight && !parameter.isStatistic) {
                options.maximizationOptions.setExcluded(parameter, true);
            }
        }
        assertBlocked(options, "node heights optimized");

        TreeLimeOptions relaxed = search(dengue(true));
        relaxed.maximizationOptions.objective = MaximizationOptions.Objective.ML;
        assertBlocked(relaxed, "use MAP or a strict clock");

        TreeLimeOptions local = search(dengue(false));
        for (PartitionClockModel clock : local.getPartitionClockModels()) {
            clock.setClockType(ClockType.SHRINKAGE_LOCAL_CLOCK);
        }
        assertBlocked(local, "shrinkage local clock");

        TreeLimeOptions unlinked = search(TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS));
        TreeLimeHarness.unlinkTreeModels(unlinked);
        assertBlocked(unlinked, "unlinked trees");

        // ML with a strict clock is fine
        TreeLimeOptions ml = search(dengue(false));
        ml.maximizationOptions.objective = MaximizationOptions.Objective.ML;
        assertNull(TopologySearchOptions.getBlocker(ml));
        TreeLimeHarness.assertIdsResolve(generate(ml));
    }

    // -------------------------------------------------------------------------------------------------------------

    static TreeLimeOptions search(TreeLimeOptions options) {
        if (options != null) {
            options.maximizationOptions.getTopologySearch().mode = TopologySearchOptions.Mode.NNI;
        }
        return options;
    }

    private static void assertBlocked(TreeLimeOptions options, String reason) throws Exception {
        options.selectParameters();
        String blocker = TopologySearchOptions.getBlocker(options);
        assertNotNull(reason, blocker);
        assertTrue(blocker, blocker.contains(reason));
        try {
            generate(options);
            fail("generated with " + blocker);
        } catch (GeneratorException e) {
            assertTrue(e.getMessage(), e.getMessage().contains(reason));
        }
    }

    private static Handler collect(final List<String> warnings) {
        Handler handler = new Handler() {
            public void publish(LogRecord record) {
                if (record.getLevel().intValue() >= Level.WARNING.intValue()) {
                    warnings.add(record.getMessage());
                }
            }

            public void flush() {
            }

            public void close() {
            }
        };
        Logger.getLogger("dr.xml").addHandler(handler);
        return handler;
    }

    static String between(String text, String start, String end) {
        int from = text.indexOf(start);
        assertTrue("no " + start, from >= 0);
        return text.substring(from, text.indexOf(end, from) + end.length());
    }

    static String generateTo(TreeLimeOptions options, File xml) throws Exception {
        options.selectParameters();
        options.selectOperators();
        new BeastGenerator(options, TreeLimeHarness.COMPONENTS).generateXML(xml);
        return new String(Files.readAllBytes(xml.toPath()), StandardCharsets.UTF_8);
    }

    static TreeLimeOptions dengue(boolean relaxed) throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(DENGUE);
        if (options == null) {
            return null;
        }
        TreeLimeHarness.setTipDatesFromNames(options);
        if (relaxed) {
            for (dr.app.treelime.options.PartitionClockModel clock : options.getPartitionClockModels()) {
                clock.setClockType(ClockType.HMC_CLOCK);
            }
        }
        return options;
    }

    // the XML with a fixed file-name stem, written into a temporary directory
    static String generate(TreeLimeOptions options) throws Exception {
        File directory = Files.createTempDirectory("treeLimeGolden").toFile();
        directory.deleteOnExit();
        File xml = new File(directory, STEM + ".xml");
        xml.deleteOnExit();
        options.fileNameStem = STEM;
        options.selectParameters();
        options.selectOperators();
        new BeastGenerator(options, TreeLimeHarness.COMPONENTS).generateXML(xml);
        return new String(Files.readAllBytes(xml.toPath()), StandardCharsets.UTF_8);
    }

    // without the version: the line that names it dropped, the beast element's attribute removed
    static String stripVersion(String xml) {
        StringBuilder sb = new StringBuilder();
        for (String line : xml.split("\n", -1)) {
            if (line.startsWith("<beast version=")) {
                sb.append("<beast>\n");
            } else if (!line.contains("Generated by TreeLime")) {
                sb.append(line).append('\n');
            }
        }
        return sb.toString();
    }

    static void assertGolden(String name, TreeLimeOptions options) throws Exception {
        if (options == null) {
            return;
        }
        final String xml = stripVersion(generate(options));
        final File file = new File(GOLDEN, name);
        if (Boolean.getBoolean("treelime.golden.update") || !file.exists() && Boolean.getBoolean("treelime.golden.create")) {
            Files.write(file.toPath(), xml.getBytes(StandardCharsets.UTF_8));
            return;
        }
        assertTrue("no golden file " + file + " (create it with -Dtreelime.golden.create=true)", file.exists());
        final String golden = new String(Files.readAllBytes(file.toPath()), StandardCharsets.UTF_8);
        if (!golden.equals(xml)) {
            String[] a = golden.split("\n", -1);
            String[] b = xml.split("\n", -1);
            int line = 0;
            while (line < Math.min(a.length, b.length) && a[line].equals(b[line])) {
                line++;
            }
            fail(name + " differs at line " + (line + 1) + ":\n golden: " + (line < a.length ? a[line] : "<end>") +
                    "\n now:    " + (line < b.length ? b[line] : "<end>"));
        }
    }

    static Map<String, XMLObject> parse(String xml, boolean strict) throws Exception {
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, strict, new BeastVersion());
        return parser.parse(new StringReader(xml), false);
    }
}
