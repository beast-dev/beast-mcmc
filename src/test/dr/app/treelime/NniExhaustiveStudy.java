/*
 * NniExhaustiveStudy.java
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
import dr.app.treelime.options.TopologySearchOptions;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.types.StartingTreeType;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.treemetrics.SplitDistance;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.search.TimeTreeSearch;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.model.Likelihood;
import dr.math.MathUtils;
import dr.xml.XMLObject;

import java.io.File;
import java.io.StringReader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * V10 of plan 7, NNI only (M2): the exhaustive gold standard on a small tree. Sequences are simulated on a coalescent
 * tree of n taxa; every rooted topology ((2n - 3)!! of them) is fitted with the TreeLime maximizer from several height
 * starts (the rescaled starting tree and random feasible heights), keeping the best F*; the spread between starts is
 * the gold standard's own error. Then the search runs from UPGMA and from random starting trees, and fails only when
 * it ends below the global optimum by more than that error.
 * <p>
 * Usage: NniExhaustiveStudy taxa sites theta heightStarts randomStarts seed; run from the repository root with BEAGLE
 * on java.library.path. JC69 sequences, contemporaneous tips, the default TreeLime model. Not a unit test.
 */
public class NniExhaustiveStudy {

    public static void main(String[] args) throws Exception {
        final int taxa = Integer.parseInt(args[0]);
        final int sites = Integer.parseInt(args[1]);
        final double theta = Double.parseDouble(args[2]);
        final int heightStarts = Integer.parseInt(args[3]);
        final int randomStarts = Integer.parseInt(args[4]);
        final long seed = Long.parseLong(args[5]);
        java.util.logging.Logger.getLogger("dr").setLevel(java.util.logging.Level.SEVERE);
        MathUtils.setSeed(seed);

        Taxa list = new Taxa();
        List<String> names = new ArrayList<String>();
        for (int i = 0; i < taxa; ++i) {
            list.addTaxon(new Taxon("t" + (i + 1)));
            names.add("t" + (i + 1));
        }
        ConstantPopulation population = new ConstantPopulation(Units.Type.SUBSTITUTIONS);
        population.setN0(theta);
        Tree truth = new CoalescentSimulator().simulateTree(list, population);
        File directory = Files.createTempDirectory("nniExhaustive").toFile();
        File fasta = new File(directory, "data.fasta");
        NniRecoveryStudy.writeFasta(NniRecoveryStudy.simulateJC(truth, sites), truth, fasta);

        // the gold standard
        List<String> topologies = topologies(names);
        System.out.println(topologies.size() + " rooted topologies of " + taxa + " taxa, " + sites + " sites");
        long started = System.nanoTime();
        double best = Double.NEGATIVE_INFINITY;
        double secondBest = Double.NEGATIVE_INFINITY;
        String bestTopology = null;
        double bestSpread = 0.0;
        double maxSpread = 0.0;
        java.util.Map<String, Double> gold = new java.util.HashMap<String, Double>();
        for (String topology : topologies) {
            double[] f = fit(fasta, topology, directory, heightStarts);
            double max = Double.NEGATIVE_INFINITY, min = Double.POSITIVE_INFINITY;
            for (double x : f) {
                max = Math.max(max, x);
                min = Math.min(min, x);
            }
            maxSpread = Math.max(maxSpread, max - min);
            Tree t = tree(topology);
            gold.put(dr.evolution.tree.TreeUtils.uniqueNewick(t, t.getRoot()), max);
            if (max > best) {
                secondBest = best;
                best = max;
                bestTopology = topology;
                bestSpread = max - min;
            } else if (max > secondBest) {
                secondBest = max;
            }
        }
        Tree bestTree = tree(bestTopology);
        System.out.printf("gold standard: F* = %.4f at %s (spread between its starts %.2e; largest spread %.2e; " +
                        "next topology %.4f), rooted RF to the true tree %d (%.0fs)%n", best, bestTopology, bestSpread,
                maxSpread, secondBest, SplitDistance.rooted(truth, bestTree), (System.nanoTime() - started) / 1E9);

        // the search from UPGMA and from random starting trees
        final double error = Math.max(0.01, bestSpread);
        int failures = 0;
        for (int start = 0; start <= randomStarts; ++start) {
            TreeLimeOptions options = TreeLimeHarness.importData(fasta.getPath());
            options.maximizationOptions.printToScreen = false;
            options.maximizationOptions.getTopologySearch().mode = TopologySearchOptions.Mode.NNI;
            options.getPartitionTreeModels().get(0).setStartingTreeType(start == 0 ? StartingTreeType.UPGMA :
                    StartingTreeType.RANDOM);
            MathUtils.setSeed(seed + start);
            TimeTreeSearch search = (TimeTreeSearch) parse(options, directory, "start" + start).get("treeSearch")
                    .getNativeObject();
            Tree startTree = new FlexibleTree(search.getTree());
            search.run();
            double gap = best - search.getFinalValue();
            boolean failed = gap > error;
            failures += failed ? 1 : 0;
            System.out.printf("%s start: F %.4f -> %.4f, gap to F* %.4f%s; rooted RF to the optimum %d -> %d; " +
                            "%d accepted; %d zero-length internal branches%n",
                    start == 0 ? "UPGMA" : "random " + start, search.getStartValue(), search.getFinalValue(), gap,
                    failed ? " FAILED" : "", SplitDistance.rooted(bestTree, startTree),
                    SplitDistance.rooted(bestTree, new FlexibleTree(search.getTree())), search.getAcceptedMoves(),
                    zeroBranches(search.getTree()));
            if (failed) {
                TreeModel found = search.getTree();
                System.out.printf("    F* of the topology it ended on: %.4f%n",
                        gold.get(dr.evolution.tree.TreeUtils.uniqueNewick(found, found.getRoot())));
                for (String row : search.getLog().getRows()) {
                    String[] fields = row.split("\t", -1);
                    System.out.println("    " + fields[1] + " " + fields[2] + " nT0 " + fields[3] + " F " + fields[7] +
                            " -> " + fields[8] + " " + fields[17]);
                }
            }
        }
        System.out.println(failures == 0 ? "all starts reached the global optimum" :
                failures + " of " + (randomStarts + 1) + " starts ended below the global optimum");
        System.exit(0);
    }

    // internal branches shorter than 1e-6 of the root height: ratios driven to 1, where NNIs change nothing
    static int zeroBranches(TreeModel tree) {
        int count = 0;
        final double root = tree.getNodeHeight(tree.getRoot());
        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
            NodeRef node = tree.getInternalNode(i);
            if (!tree.isRoot(node) &&
                    tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node) < 1E-6 * root) {
                count++;
            }
        }
        return count;
    }

    // F after the maximizer from the rescaled starting tree and from random feasible heights
    static double[] fit(File fasta, String topology, File directory, int heightStarts) throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(fasta.getPath());
        options.maximizationOptions.printToScreen = false;
        options.getPartitionTreeModels().get(0).setStartingTreeType(StartingTreeType.USER);
        options.getPartitionTreeModels().get(0).setUserStartingTree(tree(topology));
        Map<String, XMLObject> objects = parse(options, directory, "fit");
        MaximizerWrtParameter maximizer = (MaximizerWrtParameter) objects.get("maximizer").getNativeObject();
        TreeModel model = (TreeModel) objects.get("treeModel").getNativeObject();
        Likelihood joint = (Likelihood) objects.get("joint").getNativeObject();
        double[] f = new double[heightStarts];
        final double rootHeight = model.getNodeHeight(model.getRoot());
        for (int s = 0; s < heightStarts; ++s) {
            if (s > 0) {
                randomHeights(model, model.getRoot(), rootHeight);
            }
            maximizer.maximize();
            f[s] = joint.getLogLikelihood();
        }
        return f;
    }

    // a feasible height for every internal node: a uniform fraction of the way from its oldest child to its parent
    private static void randomHeights(TreeModel model, NodeRef node, double rootHeight) {
        if (model.isExternal(node)) {
            return;
        }
        model.beginTreeEdit();
        assign(model, node, rootHeight * (0.5 + MathUtils.nextDouble()));
        model.endTreeEdit();
    }

    // top-down: the node at the given height, each internal child between its subtree's oldest tip and the node
    private static void assign(TreeModel model, NodeRef node, double height) {
        model.setNodeHeight(node, height);
        for (int c = 0; c < model.getChildCount(node); ++c) {
            NodeRef child = model.getChild(node, c);
            if (!model.isExternal(child)) {
                double lower = oldestTip(model, child);
                assign(model, child, lower + (height - lower) * (0.05 + 0.9 * MathUtils.nextDouble()));
            }
        }
    }

    private static double oldestTip(TreeModel model, NodeRef node) {
        if (model.isExternal(node)) {
            return model.getNodeHeight(node);
        }
        double h = 0.0;
        for (int c = 0; c < model.getChildCount(node); ++c) {
            h = Math.max(h, oldestTip(model, model.getChild(node, c)));
        }
        return h;
    }

    static Map<String, XMLObject> parse(TreeLimeOptions options, File directory, String name) throws Exception {
        options.fileNameStem = new File(directory, name).getAbsolutePath();
        String xml = TreeLimeSearchGeneratorTest.generateTo(options, new File(directory, name + ".xml"));
        xml = xml.replace("<timeTreeSearch id=\"treeSearch\"", "<timeTreeSearch id=\"treeSearch\" searchOnParse=\"false\"");
        xml = xml.replace("<maximizeWrtParameter id=\"maximizer\" nIterations",
                "<maximizeWrtParameter id=\"maximizer\" maximizeOnParse=\"false\" nIterations");
        xml = xml.substring(0, xml.indexOf("<!-- The optimum")) + "</beast>\n";
        return new BeastParser(new String[]{}, null, false, false, false, new BeastVersion())
                .parse(new StringReader(xml), false);
    }

    // a rooted topology as a newick tree with heights equal to the number of tips below each internal node
    static Tree tree(String topology) throws Exception {
        return new NewickImporter(new StringReader(withLengths(topology) + ";")).importTree(null);
    }

    private static String withLengths(String topology) {
        return lengths(topology, new int[]{0}, Integer.MAX_VALUE).newick;
    }

    private static final class Clade {
        final String newick;
        final int tips;

        Clade(String newick, int tips) {
            this.newick = newick;
            this.tips = tips;
        }
    }

    // parses "(x,y)" recursively; each child's branch is the parent's height minus its own
    private static Clade lengths(String s, int[] at, int unused) {
        if (s.charAt(at[0]) == '(') {
            at[0]++;
            Clade left = lengths(s, at, unused);
            at[0]++; // ','
            Clade right = lengths(s, at, unused);
            at[0]++; // ')'
            int tips = left.tips + right.tips;
            double height = tips - 1;
            return new Clade("(" + left.newick + ":" + (height - heightOf(left)) + "," + right.newick + ":" +
                    (height - heightOf(right)) + ")", tips);
        }
        int end = at[0];
        while (end < s.length() && s.charAt(end) != ',' && s.charAt(end) != ')') {
            end++;
        }
        String name = s.substring(at[0], end);
        at[0] = end;
        return new Clade(name, 1);
    }

    private static double heightOf(Clade clade) {
        return clade.tips - 1;
    }

    // every rooted binary topology: each taxon in turn joins any branch, the root's included
    static List<String> topologies(List<String> names) {
        List<String> trees = new ArrayList<String>();
        trees.add("(" + names.get(0) + "," + names.get(1) + ")");
        for (int i = 2; i < names.size(); ++i) {
            List<String> next = new ArrayList<String>();
            for (String tree : trees) {
                next.addAll(insertions(tree, names.get(i)));
            }
            trees = next;
        }
        return trees;
    }

    // the tree with the taxon joined to each branch in turn: above every clade, the whole tree included
    private static List<String> insertions(String tree, String taxon) {
        List<String> result = new ArrayList<String>();
        result.add("(" + tree + "," + taxon + ")");
        if (tree.charAt(0) == '(') {
            int depth = 0;
            int comma = -1;
            for (int i = 1; i < tree.length() - 1; ++i) {
                char c = tree.charAt(i);
                if (c == '(') {
                    depth++;
                } else if (c == ')') {
                    depth--;
                } else if (c == ',' && depth == 0) {
                    comma = i;
                }
            }
            String left = tree.substring(1, comma);
            String right = tree.substring(comma + 1, tree.length() - 1);
            for (String l : insertions(left, taxon)) {
                result.add("(" + l + "," + right + ")");
            }
            for (String r : insertions(right, taxon)) {
                result.add("(" + left + "," + r + ")");
            }
        }
        return result;
    }
}
