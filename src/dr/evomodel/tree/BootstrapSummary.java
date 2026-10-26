/*
 * BootstrapSummary.java
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

package dr.evomodel.tree;

import dr.app.tools.NexusExporter;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Percentile intervals from a fixed-topology bootstrap (plan 8 S3g): the replicates' estimates (the gathered
 * bootstrap.params.tsv) and time trees (bootstrap.trees) give intervals that are written into the original analysis's
 * tables (params.tsv, dates.tsv, rates.tsv) and as height_CI, date_CI and rate_CI annotations on its time tree. Nodes
 * are matched by clade, which the fixed topology keeps the same in every replicate.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class BootstrapSummary {

    public static final String CI_METHOD = "bootstrap_percentile";
    private static final String NA = TimeTreeExporter.NA;

    /**
     * The files of one tree: the original analysis's time tree and tables, and the replicates' gathered trees.
     */
    public static final class TreeFiles {
        final File timeTree;
        final File bootstrapTrees;
        final File dates;
        final File rates; // null without a clock

        public TreeFiles(File timeTree, File bootstrapTrees, File dates, File rates) {
            this.timeTree = timeTree;
            this.bootstrapTrees = bootstrapTrees;
            this.dates = dates;
            this.rates = rates;
        }
    }

    private final File parameters;
    private final File bootstrapParameters;
    private final File summary; // null for none
    private final List<TreeFiles> trees;
    private final double level;
    private final List<String> notes = new ArrayList<String>();

    /**
     * @param level the intervals' coverage, e.g. 0.95
     */
    public BootstrapSummary(File parameters, File bootstrapParameters, File summary, List<TreeFiles> trees,
                            double level) {
        this.parameters = parameters;
        this.bootstrapParameters = bootstrapParameters;
        this.summary = summary;
        this.trees = trees;
        this.level = level;
    }

    /**
     * @return the number of replicates summarized
     */
    public int summarize() throws IOException, Importer.ImportException {
        int replicates = summarizeParameters();
        for (TreeFiles files : trees) {
            summarizeTree(files);
        }
        if (summary != null) {
            appendToSummary(replicates);
        }
        return replicates;
    }

    // ---------------------------------------------------------------------------------------------------------------

    /**
     * The percentile at p of the sorted values, interpolated between order statistics (R's type 7).
     */
    public static double percentile(double[] sorted, double p) {
        if (sorted.length == 0) {
            return Double.NaN;
        }
        double position = p * (sorted.length - 1);
        int below = (int) Math.floor(position);
        int above = Math.min(below + 1, sorted.length - 1);
        return sorted[below] + (position - below) * (sorted[above] - sorted[below]);
    }

    private double[] interval(List<Double> values) {
        double[] sorted = new double[values.size()];
        for (int i = 0; i < sorted.length; ++i) {
            sorted[i] = values.get(i);
        }
        Arrays.sort(sorted);
        return new double[]{percentile(sorted, (1.0 - level) / 2.0), percentile(sorted, (1.0 + level) / 2.0)};
    }

    private static double standardDeviation(List<Double> values) {
        if (values.size() < 2) {
            return Double.NaN;
        }
        double mean = 0.0;
        for (double value : values) {
            mean += value;
        }
        mean /= values.size();
        double sum = 0.0;
        for (double value : values) {
            sum += (value - mean) * (value - mean);
        }
        return Math.sqrt(sum / (values.size() - 1));
    }

    // ---------------------------------------------------------------------------------------------------------------
    // Parameters

    private int summarizeParameters() throws IOException {
        List<String[]> bootstrap = read(bootstrapParameters);
        if (bootstrap.size() < 2) {
            throw new IOException(bootstrapParameters + " has no replicates");
        }
        String[] names = bootstrap.get(0);
        Map<String, List<Double>> columns = new HashMap<String, List<Double>>();
        for (int j = 1; j < names.length; ++j) {
            List<Double> column = new ArrayList<Double>();
            for (String[] row : bootstrap.subList(1, bootstrap.size())) {
                if (j < row.length && isNumber(row[j])) {
                    column.add(Double.parseDouble(row[j]));
                }
            }
            columns.put(names[j], column);
        }

        List<String[]> table = read(parameters);
        Columns header = new Columns(table.get(0));
        for (String[] row : table.subList(1, table.size())) {
            List<Double> values = columns.get(row[header.index("name")]);
            if (row[header.index("status")].startsWith("fixed") || values == null || values.isEmpty()) {
                continue;
            }
            double[] interval = interval(values);
            row[header.index("se")] = format(standardDeviation(values));
            row[header.index("lower")] = Double.toString(interval[0]);
            row[header.index("upper")] = Double.toString(interval[1]);
            row[header.index("ci_method")] = CI_METHOD;
        }
        write(parameters, table);
        return bootstrap.size() - 1;
    }

    // ---------------------------------------------------------------------------------------------------------------
    // Trees

    // per clade: the heights, relative to the youngest tip, and the rates of the branch above
    private static final class CladeValues {
        final Map<Set<String>, Double> heights = new HashMap<Set<String>, Double>();
        final Map<Set<String>, Double> rates = new HashMap<Set<String>, Double>();
        final Map<Set<String>, NodeRef> nodes = new HashMap<Set<String>, NodeRef>();

        CladeValues(Tree tree) {
            double youngest = Double.POSITIVE_INFINITY;
            for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
                youngest = Math.min(youngest, tree.getNodeHeight(tree.getExternalNode(i)));
            }
            for (int i = 0; i < tree.getNodeCount(); ++i) {
                NodeRef node = tree.getNode(i);
                Set<String> clade = TreeUtils.getDescendantLeaves(tree, node);
                nodes.put(clade, node);
                heights.put(clade, tree.getNodeHeight(node) - youngest);
                Object rate = tree.getNodeAttribute(node, "rate");
                if (rate instanceof Number) {
                    rates.put(clade, ((Number) rate).doubleValue());
                }
            }
        }
    }

    private void summarizeTree(TreeFiles files) throws IOException, Importer.ImportException {
        Tree imported = readTrees(files.timeTree).get(0);
        FlexibleTree original = new FlexibleTree(imported, true);
        // the tree's own annotations (the log density) and name, which the copy does not keep
        java.util.Iterator<String> names = imported.getAttributeNames();
        while (names != null && names.hasNext()) {
            String name = names.next();
            original.setAttribute(name, imported.getAttribute(name));
        }
        original.setId(imported.getId());
        CladeValues estimate = new CladeValues(original);
        List<CladeValues> replicates = new ArrayList<CladeValues>();
        for (Tree tree : readTrees(files.bootstrapTrees)) {
            replicates.add(new CladeValues(tree));
        }

        // dates: an interval on the date (decimal year) or the age when the tips are dated, on the height otherwise
        List<String[]> dates = read(files.dates);
        Columns header = new Columns(dates.get(0));
        boolean ages = header.has("age");
        int dateColumn = header.index(ages ? "age" : "numeric_date");
        Map<String, Set<String>> clades = new HashMap<String, Set<String>>();
        Map<Set<String>, double[]> dateIntervals = new HashMap<Set<String>, double[]>();
        Map<Set<String>, double[]> heightIntervals = new HashMap<Set<String>, double[]>();
        int missing = 0;
        for (String[] row : dates.subList(1, dates.size())) {
            Set<String> clade = new HashSet<String>(Arrays.asList(row[header.index("clade")].split(",")));
            clades.put(row[header.index("#node")], clade);
            Double height = estimate.heights.get(clade);
            if (height == null) {
                continue;
            }
            List<Double> heights = new ArrayList<Double>();
            for (CladeValues replicate : replicates) {
                Double value = replicate.heights.get(clade);
                if (value == null) {
                    missing++;
                } else {
                    heights.add(value);
                }
            }
            if (heights.isEmpty()) {
                continue;
            }
            double[] heightInterval = interval(heights);
            heightIntervals.put(clade, heightInterval);
            double[] written = heightInterval;
            if (isNumber(row[dateColumn])) {
                // a date falls as the height rises, an age rises with it
                double date = Double.parseDouble(row[dateColumn]);
                double[] dateInterval = ages ?
                        new double[]{date + heightInterval[0] - height, date + heightInterval[1] - height} :
                        new double[]{date - (heightInterval[1] - height), date - (heightInterval[0] - height)};
                dateIntervals.put(clade, dateInterval);
                written = dateInterval;
            }
            row[header.index("lower")] = Double.toString(written[0]);
            row[header.index("upper")] = Double.toString(written[1]);
            row[header.index("ci_method")] = CI_METHOD;
        }
        write(files.dates, dates);
        if (missing > 0) {
            notes.add(missing + " node(s) of " + files.timeTree.getName() + " were missing from a replicate tree");
        }

        // rates, on the branch above each node
        Map<Set<String>, double[]> rateIntervals = new HashMap<Set<String>, double[]>();
        if (files.rates != null) {
            List<String[]> rates = read(files.rates);
            Columns rateHeader = new Columns(rates.get(0));
            for (String[] row : rates.subList(1, rates.size())) {
                Set<String> clade = clades.get(row[rateHeader.index("branch")]);
                if (clade == null) {
                    continue;
                }
                List<Double> values = new ArrayList<Double>();
                for (CladeValues replicate : replicates) {
                    Double value = replicate.rates.get(clade);
                    if (value != null) {
                        values.add(value);
                    }
                }
                if (values.isEmpty()) {
                    continue;
                }
                double[] rateInterval = interval(values);
                rateIntervals.put(clade, rateInterval);
                row[rateHeader.index("rate_lower")] = Double.toString(rateInterval[0]);
                row[rateHeader.index("rate_upper")] = Double.toString(rateInterval[1]);
                row[rateHeader.index("ci_method")] = CI_METHOD;
            }
            write(files.rates, rates);
        }

        // the original time tree, annotated
        for (Map.Entry<Set<String>, NodeRef> entry : estimate.nodes.entrySet()) {
            NodeRef node = entry.getValue();
            annotate(original, node, "height_CI", heightIntervals.get(entry.getKey()));
            annotate(original, node, ages ? "age_CI" : "date_CI", dateIntervals.get(entry.getKey()));
            annotate(original, node, "rate_CI", rateIntervals.get(entry.getKey()));
        }
        PrintStream out = new PrintStream(files.timeTree, "UTF-8");
        NexusExporter exporter = new NexusExporter(out);
        exporter.setSortedTranslationTable(true);
        exporter.exportTree(original);
        out.close();
    }

    private static void annotate(FlexibleTree tree, NodeRef node, String name, double[] interval) {
        if (interval != null) {
            tree.setNodeAttribute(node, name, new Object[]{interval[0], interval[1]});
        }
    }

    private static List<Tree> readTrees(File file) throws IOException, Importer.ImportException {
        FileReader reader = new FileReader(file);
        try {
            List<Tree> trees = new ArrayList<Tree>();
            NexusImporter importer = new NexusImporter(reader);
            while (importer.hasTree()) {
                trees.add(importer.importNextTree());
            }
            if (trees.isEmpty()) {
                throw new IOException(file + " has no trees");
            }
            return trees;
        } finally {
            reader.close();
        }
    }

    // ---------------------------------------------------------------------------------------------------------------
    // Run record and tables

    private void appendToSummary(int replicates) throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(summary, true));
        writer.println();
        writer.println("Bootstrap: " + replicates + " replicates (fixed topology); " +
                formatLevel() +
                "% percentile intervals (" + CI_METHOD + ") and bootstrap standard errors in the parameter table, " +
                "date (or height) intervals in the dates table, rate intervals in the rates table, and height_CI, " +
                "date_CI and rate_CI on the time tree.");
        for (String note : notes) {
            writer.println("  " + note);
        }
        writer.close();
    }

    // 95 for 0.95, 97.5 for 0.975
    private String formatLevel() {
        double percent = level * 100.0;
        return percent == Math.rint(percent) ? Long.toString((long) percent) : Double.toString(percent);
    }

    private static final class Columns {
        final List<String> names;

        Columns(String[] header) {
            names = Arrays.asList(header);
        }

        boolean has(String name) {
            return names.contains(name);
        }

        int index(String name) {
            int index = names.indexOf(name);
            if (index < 0) {
                throw new IllegalArgumentException("No column " + name + " in " + names);
            }
            return index;
        }
    }

    private static boolean isNumber(String value) {
        try {
            return !Double.isNaN(Double.parseDouble(value));
        } catch (NumberFormatException e) {
            return false;
        }
    }

    private static String format(double value) {
        return Double.isNaN(value) ? NA : Double.toString(value);
    }

    private static List<String[]> read(File file) throws IOException {
        List<String[]> rows = new ArrayList<String[]>();
        for (String line : Files.readAllLines(file.toPath())) {
            if (!line.isEmpty()) {
                rows.add(line.split("\t", -1));
            }
        }
        if (rows.isEmpty()) {
            throw new IOException(file + " is empty");
        }
        return rows;
    }

    private static void write(File file, List<String[]> rows) throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(file));
        for (String[] row : rows) {
            StringBuilder line = new StringBuilder();
            for (int i = 0; i < row.length; ++i) {
                line.append(i > 0 ? "\t" : "").append(row[i]);
            }
            writer.println(line);
        }
        writer.close();
    }
}
