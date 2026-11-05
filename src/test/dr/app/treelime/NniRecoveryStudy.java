/*
 * NniRecoveryStudy.java
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
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.treemetrics.SplitDistance;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.search.NniNeighbourhood;
import dr.evomodel.tree.search.RootedNniMove;
import dr.evomodel.tree.search.TimeTreeSearch;
import dr.math.MathUtils;
import dr.xml.XMLObject;

import java.io.File;
import java.io.PrintWriter;
import java.io.StringReader;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

/**
 * V11 of plan 7, NNI only (M2): sequences simulated on a coalescent tree; TreeLime starts from that tree perturbed by
 * k random rooted NNIs and searches. A replicate recovers the tree when the rooted RF distance to the truth is 0; it
 * reaches the truth's level when the final F is at least the F of the same search started from the true tree (the
 * data may not resolve every short branch, so the true topology need not be the best). Gate: RF = 0 for k <= 3 in at
 * least 95% of replicates.
 * <p>
 * Usage: NniRecoveryStudy taxa sites replicates kMax theta seed [firstReplicate ["search attributes"]]; run from the repository root with BEAGLE on
 * java.library.path. JC69 sequences, contemporaneous tips, the default TreeLime model (HKY, strict clock, constant
 * coalescent). Not a unit test.
 */
public class NniRecoveryStudy {

    public static void main(String[] args) throws Exception {
        final int taxa = Integer.parseInt(args[0]);
        final int sites = Integer.parseInt(args[1]);
        final int replicates = Integer.parseInt(args[2]);
        final int kMax = Integer.parseInt(args[3]);
        final double theta = Double.parseDouble(args[4]);
        final long seed = Long.parseLong(args[5]);
        final int firstRep = args.length > 6 ? Integer.parseInt(args[6]) : 0;
        searchAttributes = args.length > 7 ? args[7] : "";
        java.util.logging.Logger.getLogger("dr").setLevel(java.util.logging.Level.SEVERE);
        MathUtils.setSeed(seed);

        int[] recovered = new int[kMax + 1];
        int[] asGood = new int[kMax + 1];
        double[] startRf = new double[kMax + 1];
        double[] finalRf = new double[kMax + 1];
        System.out.println("rep\tk\trf_start\trf_final\trf_final_unrooted\tF_truth_start\tF_final\tseconds\tunderflows");
        for (int rep = firstRep; rep < replicates; ++rep) {
            // the search reseeds MathUtils (its seed attribute), so each replicate draws from its own seed
            MathUtils.setSeed(seed + 1000L * rep);
            Taxa list = new Taxa();
            for (int i = 0; i < taxa; ++i) {
                list.addTaxon(new Taxon("t" + (i + 1)));
            }
            ConstantPopulation population = new ConstantPopulation(Units.Type.SUBSTITUTIONS);
            population.setN0(theta);
            Tree truth = new CoalescentSimulator().simulateTree(list, population);
            File directory = Files.createTempDirectory("nniRecovery").toFile();
            File fasta = new File(directory, "data.fasta");
            writeFasta(simulateJC(truth, sites), truth, fasta);

            final double fTruth = search(fasta, truth, directory, "truth").getFinalValue();
            for (int k = 1; k <= kMax; ++k) {
                MathUtils.setSeed(seed + 1000L * rep + k);
                TreeModel model = new DefaultTreeModel("perturbed", truth);
                for (int m = 0; m < k; ++m) {
                    List<RootedNniMove> moves = NniNeighbourhood.enumerate(model, null);
                    moves.get(MathUtils.nextInt(moves.size())).apply(model, null);
                }
                Tree start = new FlexibleTree(model);
                long t = System.nanoTime();
                TimeTreeSearch search = search(fasta, start, directory, "k" + k);
                Tree found = new FlexibleTree(search.getTree());
                int rf0 = SplitDistance.rooted(truth, start);
                int rf = SplitDistance.rooted(truth, found);
                recovered[k] += rf == 0 ? 1 : 0;
                asGood[k] += search.getFinalValue() >= fTruth - 0.01 ? 1 : 0;
                startRf[k] += rf0;
                finalRf[k] += rf;
                System.out.printf("%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.1f\t%d%n", rep, k, rf0, rf,
                        SplitDistance.unrooted(truth, found), fTruth, search.getFinalValue(), (System.nanoTime() - t) / 1E9,
                        search.getEvaluator().getUnderflowRecomputes());
            }
        }
        System.out.println();
        System.out.println("k\tRF=0\tF>=F(truth start)\tmean RF start\tmean RF final");
        for (int k = 1; k <= kMax; ++k) {
            System.out.printf("%d\t%d/%d\t%d/%d\t%.2f\t%.2f%n", k, recovered[k], replicates, asGood[k], replicates,
                    startRf[k] / replicates, finalRf[k] / replicates);
        }
        System.exit(0);
    }

    private static String searchAttributes = "";

    // the search from the given starting tree, through TreeLime's generator
    static TimeTreeSearch search(File fasta, Tree start, File directory, String name) throws Exception {
        return search(fasta, start, directory, name, searchAttributes);
    }

    static TimeTreeSearch search(File fasta, Tree start, File directory, String name, String searchAttributes)
            throws Exception {
        TimeTreeSearch search = (TimeTreeSearch) parse(fasta, start, directory, name, searchAttributes)
                .get("treeSearch").getNativeObject();
        search.run();
        return search;
    }

    // the TreeLime XML with the search, parsed without running it or the maximizer
    static Map<String, XMLObject> parse(File fasta, Tree start, File directory, String name, String searchAttributes)
            throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(fasta.getPath());
        options.maximizationOptions.printToScreen = false;
        options.maximizationOptions.getTopologySearch().mode = TopologySearchOptions.Mode.NNI;
        options.getPartitionTreeModels().get(0).setStartingTreeType(StartingTreeType.USER);
        options.getPartitionTreeModels().get(0).setUserStartingTree(start);
        options.fileNameStem = new File(directory, name).getAbsolutePath();
        String xml = TreeLimeSearchGeneratorTest.generateTo(options, new File(directory, name + ".xml"));
        xml = xml.replace("<timeTreeSearch id=\"treeSearch\"", "<timeTreeSearch id=\"treeSearch\" searchOnParse=\"false\" " + searchAttributes);
        xml = xml.substring(0, xml.indexOf("<!-- The optimum")) + "</beast>\n";
        return new BeastParser(new String[]{}, null, false, false, false, new BeastVersion())
                .parse(new StringReader(xml), false);
    }

    // Jukes-Cantor states down the tree, one site at a time
    static int[][] simulateJC(Tree tree, int sites) {
        int[][] states = new int[tree.getNodeCount()][sites];
        for (int s = 0; s < sites; ++s) {
            states[tree.getRoot().getNumber()][s] = MathUtils.nextInt(4);
        }
        descend(tree, tree.getRoot(), states, sites);
        return states;
    }

    private static void descend(Tree tree, NodeRef node, int[][] states, int sites) {
        for (int c = 0; c < tree.getChildCount(node); ++c) {
            NodeRef child = tree.getChild(node, c);
            double t = tree.getNodeHeight(node) - tree.getNodeHeight(child);
            double change = 0.75 * (1.0 - Math.exp(-4.0 * t / 3.0)); // to any other state
            for (int s = 0; s < sites; ++s) {
                int parent = states[node.getNumber()][s];
                if (MathUtils.nextDouble() < change) {
                    int other = MathUtils.nextInt(3);
                    states[child.getNumber()][s] = other >= parent ? other + 1 : other;
                } else {
                    states[child.getNumber()][s] = parent;
                }
            }
            descend(tree, child, states, sites);
        }
    }

    static void writeFasta(int[][] states, Tree tree, File file) throws Exception {
        PrintWriter writer = new PrintWriter(file);
        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            NodeRef tip = tree.getExternalNode(i);
            writer.println(">" + tree.getNodeTaxon(tip).getId());
            StringBuilder sb = new StringBuilder();
            for (int state : states[tip.getNumber()]) {
                sb.append("ACGT".charAt(state));
            }
            writer.println(sb);
        }
        writer.close();
    }
}
