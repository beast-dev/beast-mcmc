/*
 * TopologySearchStudy.java
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
import dr.app.treelime.options.PartitionClockModel;
import dr.app.treelime.options.PartitionSubstitutionModel;
import dr.app.treelime.options.TopologySearchOptions;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.types.ClockType;
import dr.app.treelime.types.StartingTreeType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.tree.treemetrics.SplitDistance;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.search.TimeTreeSearch;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Likelihood;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import dr.xml.XMLObject;

import java.io.File;
import java.io.PrintWriter;
import java.io.StringReader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

/**
 * The M7 simulation study of plan 7 (§13.2), one dataset per call: a constant-size coalescent tree (contemporaneous
 * tips, or tips dated uniformly over the first 20% of the expected root age), a strict or uncorrelated lognormal
 * (σ = 0.3) clock, HKY (κ = 4, frequencies 0.3/0.2/0.2/0.3) with four gamma categories (α = 0.5); then TreeLime fits it
 * (HKY+Γ4, constant coalescent, the strict clock or the HMC relaxed clock) through each named pipeline, from UPGMA:
 * <ul>
 * <li>A: the fixed UPGMA topology, fitted (TreeLime without the search);</li>
 * <li>fast: NNI rounds (and root scans); standard: NNI, SPR and root-changing SPR; thorough: standard with
 * perturbation (nStop 10); standard_nj: standard from the NJ start; truth: standard from the true tree;</li>
 * <li>any label=attributes pair given in -Dvariants (semicolon-separated): the search with those attributes.</li>
 * <li>label@r<k>: that pipeline from random starting tree k (BEAST's coalescent tree, drawn with seed k + 1) with
 * search seed k, one start of a process-level multi-start.</li>
 * </ul>
 * One TSV row per pipeline: F, seconds, rooted and unrooted RF (normalised) to the truth, whether the root split is
 * right, the Kuhner-Felsenstein distance over the true tree length, the root-relative node-age RMSE on shared clades,
 * and (dated tips) the relative errors of the root age and of the mean rate (substitutions over time).
 * <p>
 * Usage: TopologySearchStudy taxa sites strict|ucln contemporaneous|dated rep pipelines [seed], with -Dout=file to
 * append the rows to a file and -DmaxRunTime=seconds to cap each search; not a unit test.
 */
public class TopologySearchStudy {

    private static final double N0 = 10.0; // time units
    private static final double RATE = 0.005; // substitutions per site per time unit
    private static final double SIGMA = 0.3;
    private static final double KAPPA = 4.0;
    private static final double[] FREQUENCIES = {0.3, 0.2, 0.2, 0.3};
    private static final double ALPHA = 0.5;

    public static void main(String[] args) throws Exception {
        if (args[0].endsWith(".xml") || args[0].endsWith(".fasta") || args[0].endsWith(".nex")) {
            realData(args);
            return;
        }
        final int taxa = Integer.parseInt(args[0]);
        final int sites = Integer.parseInt(args[1]);
        final boolean ucln = args[2].equals("ucln");
        final boolean dated = args[3].equals("dated");
        final int rep = Integer.parseInt(args[4]);
        final String[] pipelines = args[5].split(",");
        final long seed = args.length > 6 ? Long.parseLong(args[6]) : 1L;
        java.util.logging.Logger.getLogger("dr").setLevel(java.util.logging.Level.SEVERE);
        final Map<String, String> variants = variants(System.getProperty("variants", ""));
        final String cell = taxa + "\t" + sites + "\t" + args[2] + "\t" + args[3] + "\t" + rep;

        MathUtils.setSeed(seed * 1000003L + 7919L * taxa + 104729L * sites + (ucln ? 13 : 0) + (dated ? 17 : 0)
                + 1000L * rep);
        final Map<String, Double> tipHeights = new LinkedHashMap<String, Double>();
        final double expectedRoot = 2.0 * N0 * (1.0 - 1.0 / taxa);
        for (int i = 0; i < taxa; ++i) {
            tipHeights.put("t" + (i + 1), dated ? 0.2 * expectedRoot * MathUtils.nextDouble() : 0.0);
        }
        if (dated) {
            tipHeights.put("t1", 0.0); // the most recent tip at 0
        }
        final SimpleTree truth = simulateTree(tipHeights);
        final Map<SimpleNode, Double> rates = new HashMap<SimpleNode, Double>();
        drawRates((SimpleNode) truth.getRoot(), ucln, rates);
        final File directory = Files.createTempDirectory("m7study").toFile();
        final File fasta = new File(directory, "data.fasta");
        writeFasta(simulate(truth, rates, sites), truth, fasta);
        final double trueRate = meanRate(truth, rates);

        final String outName = System.getProperty("out");
        final java.io.PrintStream out = outName != null ?
                new java.io.PrintStream(new java.io.FileOutputStream(outName, true), true) : System.out;
        out.println("#taxa\tsites\tclock\ttips\trep\tpipeline\tF\tseconds\trf_rooted\trf_unrooted\troot_ok\tkf_rel" +
                "\tage_rmse_rel\troot_age_err\trate_err\trounds\taccepted\tstop");
        for (String pipeline : pipelines) {
            final long started = System.nanoTime();
            final Result result = run(pipeline, variants, fasta, tipHeights, dated, ucln, truth, directory);
            final double seconds = (System.nanoTime() - started) / 1E9;
            final Tree found = result.tree;
            final double rootTruth = truth.getNodeHeight(truth.getRoot());
            final double rootFound = found.getNodeHeight(found.getRoot());
            out.printf("%s\t%s\t%.4f\t%.1f\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t%s\t%s\t%d\t%d\t%s%n", cell, pipeline,
                    result.f, seconds, SplitDistance.normalizedRooted(truth, found),
                    SplitDistance.normalizedUnrooted(truth, found), rootSplit(truth).equals(rootSplit(found)) ? 1 : 0,
                    kuhnerFelsenstein(truth, found) / treeLength(truth), ageRmse(truth, found),
                    dated ? String.format("%.4f", rootFound / rootTruth - 1.0) : "NA",
                    dated ? String.format("%.4f", result.rate / trueRate - 1.0) : "NA",
                    result.rounds, result.accepted, result.stop);
        }
        System.exit(0);
    }

    private static Map<String, String> variants(String text) {
        Map<String, String> map = new LinkedHashMap<String, String>();
        for (String part : text.split(";")) {
            final int at = part.indexOf('=');
            if (at > 0) {
                map.put(part.substring(0, at).trim(), part.substring(at + 1).trim());
            }
        }
        return map;
    }

    private static final class Result {
        double f;
        Tree tree;
        double rate;
        int rounds;
        int accepted;
        String stop = "NA";
    }

    // real data (§13.3): path strict|hmc pipelines; one row per pipeline with F, time, the unrooted and rooted RF to
    // pipeline A's tree (the fixed starting topology), the root age and the mean rate
    private static void realData(String[] args) throws Exception {
        java.util.logging.Logger.getLogger("dr").setLevel(java.util.logging.Level.SEVERE);
        final File data = new File(args[0]);
        final boolean hmc = args[1].equals("hmc");
        final Map<String, String> variants = variants(System.getProperty("variants", ""));
        final String outName = System.getProperty("out");
        final java.io.PrintStream out = outName != null ?
                new java.io.PrintStream(new java.io.FileOutputStream(outName, true), true) : System.out;
        final File directory = Files.createTempDirectory("m7real").toFile();
        out.println("#data\tclock\tpipeline\tF\tseconds\trf_rooted_to_A\trf_unrooted_to_A\troot_age\tmean_rate" +
                "\trounds\taccepted\tstop");
        Tree fixed = null;
        for (String pipeline : args[2].split(",")) {
            TreeLimeOptions options = TreeLimeHarness.importData(data.getPath());
            for (PartitionClockModel clock : options.getPartitionClockModels()) {
                clock.setClockType(hmc ? ClockType.HMC_CLOCK : ClockType.STRICT_CLOCK);
            }
            final long started = System.nanoTime();
            final Result result = run(pipeline, variants, options, null, directory);
            final double seconds = (System.nanoTime() - started) / 1E9;
            if (pipeline.equals("A")) {
                fixed = result.tree;
            }
            out.printf("%s\t%s\t%s\t%.4f\t%.1f\t%s\t%s\t%.4f\t%.6g\t%d\t%d\t%s%n", data.getName(), args[1], pipeline,
                    result.f, seconds,
                    fixed != null ? String.format("%.4f", SplitDistance.normalizedRooted(fixed, result.tree)) : "NA",
                    fixed != null ? String.format("%.4f", SplitDistance.normalizedUnrooted(fixed, result.tree)) : "NA",
                    result.tree.getNodeHeight(result.tree.getRoot()), result.rate, result.rounds, result.accepted,
                    result.stop);
        }
        System.exit(0);
    }

    private static Result run(String pipeline, Map<String, String> variants, File fasta, Map<String, Double> heights,
                              boolean dated, boolean ucln, Tree truth, File directory) throws Exception {
        if (pipeline.equals("E") || pipeline.equals("F")) {
            return external(pipeline, fasta, heights, dated, sitesOf(fasta), directory);
        }
        TreeLimeOptions options = TreeLimeHarness.importData(fasta.getPath());
        if (dated) {
            for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
                Taxon taxon = options.taxonList.getTaxon(i);
                final double height = heights.get(taxon.getId());
                taxon.setDate(new Date(100.0 - height, Units.Type.YEARS, false));
                taxon.setAttribute("_height", height);
            }
            options.useTipDates = true;
            double maximum = 0.0;
            for (double h : heights.values()) {
                maximum = Math.max(maximum, h);
            }
            options.maximumTipHeight = maximum;
        }
        for (PartitionClockModel clock : options.getPartitionClockModels()) {
            clock.setClockType(ucln ? ClockType.HMC_CLOCK : ClockType.STRICT_CLOCK);
        }
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
            model.setGammaHetero(true);
            model.setGammaCategories(4);
        }
        return run(pipeline, variants, options, truth, directory);
    }

    // the pipeline on prepared options (data, dates, models)
    private static Result run(String pipeline, Map<String, String> variants, TreeLimeOptions options, Tree truth,
                              File directory) throws Exception {
        String attributes;
        boolean search = true;
        boolean fromTruth = false;
        // label@r<k>: that pipeline from random starting tree k (a coalescent tree drawn with seed k), one start of a
        // process-level multi-start (H3)
        int randomStart = -1;
        final int at = pipeline.indexOf("@r");
        if (at > 0) {
            randomStart = Integer.parseInt(pipeline.substring(at + 2));
            pipeline = pipeline.substring(0, at);
        }
        if (pipeline.equals("A")) {
            search = false;
            attributes = "";
        } else if (pipeline.equals("fast")) {
            attributes = "moves=\"nni\"";
        } else if (pipeline.equals("standard")) {
            attributes = "moves=\"nni,spr,root\"";
        } else if (pipeline.equals("thorough")) {
            attributes = "moves=\"nni,spr,root\" nStop=\"10\"";
        } else if (pipeline.equals("standard_nj")) {
            attributes = "moves=\"nni,spr,root\" startingTree=\"nj\"";
        } else if (pipeline.equals("truth") && truth != null) {
            attributes = "moves=\"nni,spr,root\"";
            fromTruth = true;
        } else if (variants.containsKey(pipeline)) {
            attributes = variants.get(pipeline);
        } else {
            throw new IllegalArgumentException("Unknown pipeline " + pipeline);
        }
        final String cap = System.getProperty("maxRunTime");
        if (search && cap != null) {
            attributes += " maxRunTime=\"" + cap + "\"";
        }

        options.maximizationOptions.printToScreen = false;
        if (System.getProperty("objective") != null) { // -Dobjective=ML: the tree likelihoods only (H2's bias check)
            options.maximizationOptions.objective =
                    dr.app.treelime.options.MaximizationOptions.Objective.valueOf(System.getProperty("objective"));
        }
        options.maximizationOptions.getTopologySearch().mode = search ? TopologySearchOptions.Mode.NNI :
                TopologySearchOptions.Mode.OFF;
        if (fromTruth) {
            options.getPartitionTreeModels().get(0).setStartingTreeType(StartingTreeType.USER);
            options.getPartitionTreeModels().get(0).setUserStartingTree(new FlexibleTree(truth));
        }
        if (randomStart >= 0) {
            options.getPartitionTreeModels().get(0).setStartingTreeType(StartingTreeType.RANDOM);
            options.maximizationOptions.getTopologySearch().seed = randomStart;
            MathUtils.setSeed(randomStart + 1L);
            pipeline += "_r" + randomStart;
        }
        options.fileNameStem = new File(directory, pipeline).getAbsolutePath();
        String xml = TreeLimeSearchGeneratorTest.generateTo(options, new File(directory, pipeline + ".xml"));
        if (search) {
            xml = xml.replace("<timeTreeSearch id=\"treeSearch\" ", "<timeTreeSearch id=\"treeSearch\" " + attributes +
                    " ");
        }
        xml = xml.substring(0, xml.indexOf("<!-- The optimum")) + "</beast>\n";
        Map<String, XMLObject> objects = new BeastParser(new String[]{}, null, false, false, false,
                new BeastVersion()).parse(new StringReader(xml), false);

        Result result = new Result();
        TreeModel tree = (TreeModel) objects.get("treeModel").getNativeObject();
        if (search) {
            TimeTreeSearch timeTreeSearch = (TimeTreeSearch) objects.get("treeSearch").getNativeObject();
            result.f = timeTreeSearch.getFinalValue();
            result.rounds = timeTreeSearch.getRounds();
            result.accepted = timeTreeSearch.getAcceptedMoves();
            result.stop = timeTreeSearch.getStopReason().name();
        } else {
            result.f = ((Likelihood) objects.get("joint").getNativeObject()).getLogLikelihood();
        }
        result.tree = new FlexibleTree(tree);
        BranchRateModel branchRates = ((TreeDataLikelihood) objects.get("treeLikelihood").getNativeObject())
                .getBranchRateModel();
        double substitutions = 0.0, time = 0.0;
        for (int n = 0; n < tree.getNodeCount(); ++n) {
            NodeRef node = tree.getNode(n);
            if (!tree.isRoot(node)) {
                final double t = tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node);
                substitutions += branchRates.getBranchRate(tree, node) * t;
                time += t;
            }
        }
        result.rate = substitutions / time;
        return result;
    }

    // -------------------------------------------------------------------------------------------------------------
    // H2's comparators (§13.2), run from -Dtools=<dir> holding iqtree3, lsd2 and tt-venv/bin/treetime: IQ-TREE's ML
    // tree under the true model family (HKY+G4), then dated by LSD2 (E: -r a, all branches, no collapsing; without tip
    // dates the root at 0 and the tips at 1) or by TreeTime (F: --coalescent opt --reroot best; dated tips only). F is
    // NaN: the comparators do not report TreeLime's objective

    private static Result external(String pipeline, File fasta, Map<String, Double> heights, boolean dated, int sites,
                                   File directory) throws Exception {
        final String tools = System.getProperty("tools");
        if (tools == null) {
            throw new IllegalArgumentException("pipeline " + pipeline + " needs -Dtools=<dir>");
        }
        if (pipeline.equals("F") && !dated) {
            throw new IllegalArgumentException("TreeTime (F) needs tip dates");
        }
        final File iqTree = new File(directory, "iq.treefile");
        if (!iqTree.exists()) {
            execute(directory, tools + "/iqtree3", "-s", fasta.getPath(), "-m", "HKY+G4", "-T", "1", "-seed", "1",
                    "--prefix", new File(directory, "iq").getPath(), "-quiet", "-redo");
        }
        // dates forward in time, the most recent tip at 100 (as the TreeLime pipelines' taxon dates)
        PrintWriter lsdDates = new PrintWriter(new File(directory, "dates.txt"));
        PrintWriter csvDates = new PrintWriter(new File(directory, "dates.csv"));
        lsdDates.println(heights.size());
        csvDates.println("name,date");
        for (Map.Entry<String, Double> tip : heights.entrySet()) {
            lsdDates.println(tip.getKey() + " " + (100.0 - tip.getValue()));
            csvDates.println(tip.getKey() + "," + (100.0 - tip.getValue()));
        }
        lsdDates.close();
        csvDates.close();

        Result result = new Result();
        result.f = Double.NaN;
        File treeFile;
        if (pipeline.equals("E")) {
            List<String> command = new ArrayList<String>(java.util.Arrays.asList(tools + "/lsd2", "-i",
                    iqTree.getPath(), "-s", Integer.toString(sites), "-r", "a", "-l", "-1",
                    "-o", new File(directory, "lsd").getPath()));
            if (dated) {
                command.add("-d");
                command.add(new File(directory, "dates.txt").getPath());
            }
            execute(directory, command.toArray(new String[0]));
            treeFile = new File(directory, "lsd.date.nexus"); // lsd.nwk is the rooted input tree, in substitutions
            result.rate = firstNumberAfter(new File(directory, "lsd"), "rate "); // " rate 0.0042, tMRCA 82.4 , ..."
        } else {
            final File out = new File(directory, "tt");
            execute(directory, tools + "/tt-venv/bin/treetime", "--tree", iqTree.getPath(), "--aln", fasta.getPath(),
                    "--dates", new File(directory, "dates.csv").getPath(), "--coalescent", "opt", "--reroot", "best",
                    "--outdir", out.getPath());
            treeFile = new File(out, "timetree.nexus");
            result.rate = firstNumberAfter(new File(out, "molecular_clock.txt"), "--rate:");
        }
        java.io.FileReader reader = new java.io.FileReader(treeFile);
        result.tree = new FlexibleTree(treeFile.getName().endsWith(".nexus") ?
                new dr.evolution.io.NexusImporter(reader).importNextTree() :
                new dr.evolution.io.NewickImporter(reader).importNextTree());
        reader.close();
        return result;
    }

    private static int sitesOf(File fasta) throws Exception {
        final List<String> lines = Files.readAllLines(fasta.toPath());
        return lines.size() > 1 ? lines.get(1).trim().length() : 0;
    }

    private static void execute(File directory, String... command) throws Exception {
        ProcessBuilder builder = new ProcessBuilder(command).directory(directory).redirectErrorStream(true)
                .redirectOutput(new File(directory, new File(command[0]).getName() + ".log"));
        builder.environment().put("OMP_TOOL", "disabled"); // IQ-TREE's OpenMP otherwise looks for libarcher.so
        final int status = builder.start().waitFor();
        if (status != 0) {
            throw new RuntimeException(String.join(" ", command) + " exited with " + status + "; see " +
                    new File(directory, new File(command[0]).getName() + ".log"));
        }
    }

    // the number after the key on the first line of a tool's text output that starts with it
    private static double firstNumberAfter(File file, String key) throws Exception {
        for (String line : Files.readAllLines(file.toPath())) {
            final int at = line.indexOf(key);
            if (line.trim().startsWith(key)) {
                java.util.regex.Matcher m = java.util.regex.Pattern.compile("[-+]?[0-9.]+(?:[eE][-+]?[0-9]+)?")
                        .matcher(line.substring(at + key.length()));
                if (m.find()) {
                    return Double.parseDouble(m.group());
                }
            }
        }
        return Double.NaN;
    }

    // a constant-size coalescent with tips at the given heights
    static SimpleTree simulateTree(Map<String, Double> tipHeights) {
        List<SimpleNode> samples = new ArrayList<SimpleNode>();
        for (Map.Entry<String, Double> tip : tipHeights.entrySet()) {
            SimpleNode node = new SimpleNode();
            node.setTaxon(new Taxon(tip.getKey()));
            node.setHeight(tip.getValue());
            samples.add(node);
        }
        samples.sort((a, b) -> Double.compare(a.getHeight(), b.getHeight()));
        List<SimpleNode> active = new ArrayList<SimpleNode>();
        int next = 0;
        double t = 0.0;
        while (next < samples.size() || active.size() > 1) {
            while (next < samples.size() && samples.get(next).getHeight() <= t) {
                active.add(samples.get(next++));
            }
            final int k = active.size();
            final double wait = k >= 2 ? MathUtils.nextExponential(k * (k - 1) / (2.0 * N0)) : Double.POSITIVE_INFINITY;
            if (next < samples.size() && t + wait > samples.get(next).getHeight()) {
                t = samples.get(next).getHeight();
                continue;
            }
            t += wait;
            SimpleNode a = active.remove(MathUtils.nextInt(active.size()));
            SimpleNode b = active.remove(MathUtils.nextInt(active.size()));
            SimpleNode parent = new SimpleNode();
            parent.setHeight(t);
            parent.addChild(a);
            parent.addChild(b);
            active.add(parent);
        }
        return new SimpleTree(active.get(0));
    }

    private static void drawRates(SimpleNode node, boolean ucln, Map<SimpleNode, Double> rates) {
        for (int c = 0; c < node.getChildCount(); ++c) {
            SimpleNode child = node.getChild(c);
            rates.put(child, ucln ? RATE * Math.exp(SIGMA * MathUtils.nextGaussian() - 0.5 * SIGMA * SIGMA) : RATE);
            drawRates(child, ucln, rates);
        }
    }

    private static double meanRate(SimpleTree tree, Map<SimpleNode, Double> rates) {
        double substitutions = 0.0, time = 0.0;
        for (Map.Entry<SimpleNode, Double> entry : rates.entrySet()) {
            final double t = entry.getKey().getParent().getHeight() - entry.getKey().getHeight();
            substitutions += entry.getValue() * t;
            time += t;
        }
        return substitutions / time;
    }

    // HKY + Γ4 states down the tree
    private static int[][] simulate(SimpleTree tree, Map<SimpleNode, Double> rates, int sites) {
        final HKY hky = new HKY(KAPPA, new FrequencyModel(Nucleotides.INSTANCE, FREQUENCIES));
        final double[] categories = new double[4];
        final GammaDistribution gamma = new GammaDistribution(ALPHA, 1.0 / ALPHA);
        double sum = 0.0;
        for (int c = 0; c < 4; ++c) {
            categories[c] = gamma.quantile((2.0 * c + 1.0) / 8.0);
            sum += categories[c];
        }
        for (int c = 0; c < 4; ++c) {
            categories[c] *= 4.0 / sum;
        }
        final int[] category = new int[sites];
        final int[] rootStates = new int[sites];
        for (int s = 0; s < sites; ++s) {
            category[s] = MathUtils.nextInt(4);
            rootStates[s] = draw(FREQUENCIES, 0);
        }
        final int[][] states = new int[tree.getNodeCount()][];
        final Map<SimpleNode, int[]> byNode = new HashMap<SimpleNode, int[]>();
        byNode.put((SimpleNode) tree.getRoot(), rootStates);
        descend(hky, (SimpleNode) tree.getRoot(), rates, categories, category, byNode);
        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            states[i] = byNode.get((SimpleNode) tree.getExternalNode(i));
        }
        return states;
    }

    private static void descend(HKY hky, SimpleNode node, Map<SimpleNode, Double> rates, double[] categories,
                                int[] category, Map<SimpleNode, int[]> byNode) {
        final int[] parent = byNode.get(node);
        for (int c = 0; c < node.getChildCount(); ++c) {
            final SimpleNode child = node.getChild(c);
            final double distance = rates.get(child) * (node.getHeight() - child.getHeight());
            final double[][] matrices = new double[4][16];
            for (int k = 0; k < 4; ++k) {
                hky.getTransitionProbabilities(distance * categories[k], matrices[k]);
            }
            final int[] states = new int[parent.length];
            for (int s = 0; s < parent.length; ++s) {
                states[s] = draw(matrices[category[s]], 4 * parent[s]);
            }
            byNode.put(child, states);
            descend(hky, child, rates, categories, category, byNode);
        }
    }

    private static int draw(double[] probabilities, int offset) {
        double u = MathUtils.nextDouble();
        for (int k = 0; k < 3; ++k) {
            u -= probabilities[offset + k];
            if (u < 0) {
                return k;
            }
        }
        return 3;
    }

    private static void writeFasta(int[][] states, Tree tree, File file) throws Exception {
        PrintWriter writer = new PrintWriter(file);
        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            writer.println(">" + tree.getNodeTaxon(tree.getExternalNode(i)).getId());
            StringBuilder sb = new StringBuilder();
            for (int state : states[i]) {
                sb.append("ACGT".charAt(state));
            }
            writer.println(sb);
        }
        writer.close();
    }

    // -------------------------------------------------------------------------------------------------------------
    // metrics on rooted clades, keyed by their sorted taxa

    private static Map<String, double[]> clades(Tree tree) {
        Map<String, double[]> map = new HashMap<String, double[]>(); // {height, branch length above (0 at the root)}
        fill(tree, tree.getRoot(), map);
        return map;
    }

    private static TreeSet<String> fill(Tree tree, NodeRef node, Map<String, double[]> map) {
        TreeSet<String> taxa = new TreeSet<String>();
        if (tree.isExternal(node)) {
            taxa.add(tree.getNodeTaxon(node).getId());
        } else {
            for (int c = 0; c < tree.getChildCount(node); ++c) {
                taxa.addAll(fill(tree, tree.getChild(node, c), map));
            }
        }
        final double branch = tree.isRoot(node) ? 0.0 :
                tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node);
        map.put(taxa.toString(), new double[]{tree.getNodeHeight(node), branch, tree.isExternal(node) ? 1 : 0});
        return taxa;
    }

    private static String rootSplit(Tree tree) {
        TreeSet<String> a = new TreeSet<String>(), b = new TreeSet<String>();
        fill(tree, tree.getChild(tree.getRoot(), 0), new HashMap<String, double[]>()).forEach(a::add);
        fill(tree, tree.getChild(tree.getRoot(), 1), new HashMap<String, double[]>()).forEach(b::add);
        return a.first().compareTo(b.first()) < 0 ? a + "|" + b : b + "|" + a;
    }

    // over the union of rooted clades (tips included), the branch lengths above them, absent clades counting 0
    private static double kuhnerFelsenstein(Tree truth, Tree found) {
        Map<String, double[]> x = clades(truth), y = clades(found);
        double sum = 0.0;
        for (Map.Entry<String, double[]> e : x.entrySet()) {
            final double other = y.containsKey(e.getKey()) ? y.get(e.getKey())[1] : 0.0;
            sum += (e.getValue()[1] - other) * (e.getValue()[1] - other);
        }
        for (Map.Entry<String, double[]> e : y.entrySet()) {
            if (!x.containsKey(e.getKey())) {
                sum += e.getValue()[1] * e.getValue()[1];
            }
        }
        return Math.sqrt(sum);
    }

    private static double treeLength(Tree tree) {
        double length = 0.0;
        for (double[] clade : clades(tree).values()) {
            length += clade[1];
        }
        return length;
    }

    // internal node ages over each tree's own root age, on the clades both trees have
    private static double ageRmse(Tree truth, Tree found) {
        Map<String, double[]> x = clades(truth), y = clades(found);
        final double rx = truth.getNodeHeight(truth.getRoot()), ry = found.getNodeHeight(found.getRoot());
        double sum = 0.0;
        int count = 0;
        for (Map.Entry<String, double[]> e : x.entrySet()) {
            if (e.getValue()[2] == 0 && y.containsKey(e.getKey())) {
                final double d = e.getValue()[0] / rx - y.get(e.getKey())[0] / ry;
                sum += d * d;
                count++;
            }
        }
        return count > 0 ? Math.sqrt(sum / count) : Double.NaN;
    }
}
