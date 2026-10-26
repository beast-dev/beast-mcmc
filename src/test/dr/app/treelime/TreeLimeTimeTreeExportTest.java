/*
 * TreeLimeTimeTreeExportTest.java
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
import dr.app.treelime.options.MaximizationOptions;
import dr.app.treelime.options.PartitionTreeModel;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.types.StartingTreeType;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evomodel.tree.TimeTreeExporter;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.model.Likelihood;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.File;
import java.io.FileReader;
import java.io.StringReader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * The files timeTreeExport writes after a TreeLime maximization (plan 8 O1): the trees BEAST's tree loggers write
 * read back with the model's heights and rates; the tables agree with them; the run record's log likelihoods add up and match the
 * maximizer's objective.
 */
public class TreeLimeTimeTreeExportTest extends TestCase {

    private static final double TOLERANCE = 1E-10;

    public TreeLimeTimeTreeExportTest(String name) {
        super(name);
    }

    public void testIsoDate() {
        assertEquals("2013-01-01", TimeTreeExporter.isoDate(2013.0));
        assertEquals("2013-07-02", TimeTreeExporter.isoDate(2013.5));
        assertEquals("2012-12-31", TimeTreeExporter.isoDate(2012.999));
        assertEquals(TimeTreeExporter.NA, TimeTreeExporter.isoDate(Double.NaN));
    }

    // dated tips, MAP, node heights optimized: the trees, the dates and rates tables, and logL + logPrior = -F
    public void testDatedTree() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData("examples/Data/Dengue4.env.nex");
        if (options == null) {
            return;
        }
        TreeLimeHarness.setTipDatesFromNames(options);
        Export export = run(quiet(options));

        assertEquals("MAP", export.exporter.objectiveName());
        assertEquals(-export.maximizer.getMinimumValue(), export.logL + export.logPrior, 1E-8);
        assertTrue(export.summary.contains("Objective: MAP"));
        assertTrue(export.summary.contains(TimeTreeExporter.STATEMENT));

        Tree model = export.tree("treeModel");
        Tree time = read(export.file(".timetree.nex"));
        Tree subst = read(export.file(".subst.nex"));
        assertEquals(model.getNodeCount(), time.getNodeCount());
        // the time tree's heights, relative to the youngest tip (D4ElSal94, at height 0), are the model's
        Map<java.util.Set<String>, Double> modelHeights = heightsByClade(model);
        Map<java.util.Set<String>, Double> heights = heightsByClade(time);
        double youngest = heights.get(java.util.Collections.singleton("D4ElSal94"));
        assertEquals(modelHeights.keySet(), heights.keySet());
        for (Map.Entry<java.util.Set<String>, Double> clade : modelHeights.entrySet()) {
            double height = clade.getValue();
            assertEquals(height, heights.get(clade.getKey()) - youngest, TOLERANCE * Math.max(1.0, height));
        }
        // the substitution tree's lengths are the time tree's, times the rates it is annotated with
        double length = 0.0;
        for (int i = 0; i < time.getNodeCount(); i++) {
            NodeRef node = time.getNode(i);
            if (!time.isRoot(node)) {
                length += (Double) time.getNodeAttribute(node, "rate") *
                        (time.getNodeHeight(time.getParent(node)) - time.getNodeHeight(node));
            }
        }
        assertEquals(length, treeLength(subst), 1E-12);

        // numeric date + height is the date of height 0 everywhere
        List<String[]> dates = table(export.file(".dates.tsv"));
        assertEquals("#node", dates.get(0)[0]);
        assertEquals(model.getNodeCount() + 1, dates.size());
        boolean foundRoot = false;
        for (String[] row : dates.subList(1, dates.size())) {
            assertEquals(1994.0, Double.parseDouble(row[3]) + Double.parseDouble(row[4]), 1E-9);
            assertEquals(TimeTreeExporter.isoDate(Double.parseDouble(row[3])), row[2]);
            if (row[1].split(",").length == model.getExternalNodeCount()) {
                assertEquals(model.getNodeHeight(model.getRoot()), Double.parseDouble(row[4]), 1E-12);
                foundRoot = true;
            }
        }
        assertTrue(foundRoot);
        assertEquals(model.getNodeCount(), table(export.file(".rates.tsv")).size()); // a header, no root row

        List<String[]> parameters = table(export.file(".params.tsv"));
        assertEquals(export.maximizer.getParameter().getDimension() + 1, parameters.size());
        int dependent = 0;
        for (String[] row : parameters.subList(1, parameters.size())) {
            if (row[4].equals("dependent")) {
                assertEquals("UnitSimplexTransform", row[3]);
                dependent++;
            } else {
                assertTrue(row[4], row[4].equals("optimized") || row[4].equals("held_at_bound"));
            }
        }
        assertEquals(1, dependent); // the frequencies' last coordinate
        assertEquals(export.maximizer.getParameter().getDimension() - 1, export.maximizer.getFreeCoordinateCount());
        assertTrue(export.maximizer.getReport().contains("Optimum:\n  kappa = "));
    }

    // no tip dates, two partitions, ML: one log likelihood per partition, logL = -F, dates NA, a taxon set's tMRCA
    public void testPartitionsUnderMaximumLikelihood() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        Taxa clade = new Taxa("clade");
        clade.addTaxon(options.taxonList.getTaxon(0));
        clade.addTaxon(options.taxonList.getTaxon(1));
        options.taxonSets.add(clade);
        options.taxonSetsTreeModel.put(clade, options.getPartitionTreeModels().get(0));
        options.taxonSetsMono.put(clade, false);
        options.taxonSetsIncludeStem.put(clade, false);
        options.maximizationOptions.objective = MaximizationOptions.Objective.ML;
        Export export = run(quiet(options));

        assertEquals("ML", export.exporter.objectiveName());
        assertEquals(-export.maximizer.getMinimumValue(), export.logL, 1E-8);

        List<TimeTreeExporter.PartitionTerm> terms = export.exporter.getPartitionTerms();
        assertEquals(2, terms.size());
        assertEquals("firsthalf.patterns", terms.get(0).name);
        assertEquals("secondhalf.patterns", terms.get(1).name);
        assertEquals(export.logL, terms.get(0).logLikelihood + terms.get(1).logLikelihood, 1E-8);
        assertTrue(export.summary.contains("n (weighted sites in the likelihood): " +
                (long) (terms.get(0).sites + terms.get(1).sites)));

        List<String[]> dates = table(export.file(".dates.tsv"));
        String[] last = dates.get(dates.size() - 1);
        assertEquals("tmrca(clade)", last[0]);
        assertEquals(TimeTreeExporter.NA, last[2]);
        assertEquals(TimeTreeExporter.NA, last[3]);
        // the clock rate is fixed without tip dates
        assertTrue(export.summary.contains("fixed in the Priors tab: clock.rate"));
    }

    // -------------------------------------------------------------------------------------------------------------

    private static TreeLimeOptions quiet(TreeLimeOptions options) {
        for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
            tree.setStartingTreeType(StartingTreeType.UPGMA);
        }
        options.maximizationOptions.nIterations = 300;
        options.maximizationOptions.printToScreen = false;
        return options;
    }

    private static class Export {
        File directory;
        Map<String, XMLObject> objects;
        MaximizerWrtParameter maximizer;
        TimeTreeExporter exporter;
        double logL;
        double logPrior;
        String summary;

        File file(String suffix) {
            return new File(directory, "run" + suffix);
        }

        Tree tree(String id) {
            return (Tree) objects.get(id).getNativeObject();
        }
    }

    // writes the XML with its files in a new directory and runs all of it, the export included
    private static Export run(TreeLimeOptions options) throws Exception {
        String xml = TreeLimeHarness.generate(options);
        Matcher stem = Pattern.compile("fileNameStem=\"([^\"]+)\"").matcher(xml);
        assertTrue(stem.find());
        Export export = new Export();
        export.directory = Files.createTempDirectory("timeTreeExport").toFile();
        export.directory.deleteOnExit();
        xml = xml.replace(stem.group(1), new File(export.directory, "run").getAbsolutePath());

        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        export.objects = parser.parse(new StringReader(xml), false);
        export.maximizer = (MaximizerWrtParameter) export.objects.get("maximizer").getNativeObject();
        export.exporter = (TimeTreeExporter) export.objects.get("timeTreeExport").getNativeObject();
        export.logL = ((Likelihood) export.objects.get("likelihood").getNativeObject()).getLogLikelihood();
        export.logPrior = ((Likelihood) export.objects.get("prior").getNativeObject()).getLogLikelihood();
        export.summary = new String(Files.readAllBytes(export.file(".summary.txt").toPath()));
        for (File file : export.directory.listFiles()) {
            file.deleteOnExit();
        }
        return export;
    }

    private static Tree read(File file) throws Exception {
        FileReader reader = new FileReader(file);
        Tree tree = new NexusImporter(reader).importNextTree();
        reader.close();
        return tree;
    }

    private static Map<java.util.Set<String>, Double> heightsByClade(Tree tree) {
        Map<java.util.Set<String>, Double> heights = new java.util.HashMap<java.util.Set<String>, Double>();
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            heights.put(TreeUtils.getDescendantLeaves(tree, node), tree.getNodeHeight(node));
        }
        return heights;
    }

    private static double treeLength(Tree tree) {
        double length = 0.0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                length += tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node);
            }
        }
        return length;
    }

    private static List<String[]> table(File file) throws Exception {
        List<String[]> rows = new ArrayList<String[]>();
        for (String line : Files.readAllLines(file.toPath())) {
            rows.add(line.split("\t", -1));
        }
        return rows;
    }
}
