/*
 * TreeSummarizer.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
 */

package dr.app.tools;

import dr.util.Version;
import dr.util.HeapSort;
import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.Importer;
import dr.evolution.io.TreeImporter;
import dr.evolution.util.TaxonList;
import dr.stats.DiscreteStatistics;
import dr.evomodel.tree.TreeTraceAnalysis;

import java.io.*;
import java.util.*;

import org.virion.jam.console.ConsoleApplication;

import javax.swing.*;

public class TreeSummarizer {

    private final static Version version = new BeastVersion();

    public final static int KEEP_HEIGHTS = 0;
    public final static int MEAN_HEIGHTS = 1;
    public final static int MEDIAN_HEIGHTS = 2;

    public TreeSummarizer(int burnin, String inputFileName, String outputFileName) throws IOException {

        System.out.println("Reading trees...");

        CladeSystem cladeSystem = new CladeSystem();

        int fileTotalTrees = 0;
        TreeImporter importer = new NexusImporter(new FileReader(inputFileName));
        try {
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (fileTotalTrees >= burnin) {
                    cladeSystem.add(tree);

                    totalTreesUsed += 1;
                }
                totalTrees += 1;
            }
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }

        System.out.println("Total trees read: " + totalTrees);
        System.out.println("Total trees summarized: " + totalTreesUsed);

        System.out.println("Finding summary tree...");
        Tree summaryTree = summarizeTrees(burnin, cladeSystem, inputFileName);

        System.out.println("Writing summary tree...");
        if (outputFileName != null) {
            NexusExporter exporter = new NexusExporter(new PrintStream(new FileOutputStream(outputFileName)));
            exporter.exportTree(summaryTree);
        } else {
            NexusExporter exporter = new NexusExporter(System.out);
            exporter.exportTree(summaryTree);
        }
        System.out.println("Finished.");
    }

    private Tree summarizeTrees(int burnin, CladeSystem cladeSystem, String inputFileName) throws IOException {

        Tree bestTree = null;
        double bestScore = 0.0;

        System.out.println("Analyzing " + totalTreesUsed + " trees...");
        System.out.println("0              25             50             75            100");
        System.out.println("|--------------|--------------|--------------|--------------|");
        System.out.print(  "*");

        double stepSize = totalTreesUsed/60.0;
        int counter = 0;

        int fileTotalTrees = 0;
        TreeImporter importer = new NexusImporter(new FileReader(inputFileName));
        try {
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (fileTotalTrees >= burnin) {
                    double score = scoreTree(tree, cladeSystem);
//                    System.out.println(score);
                    if (score > bestScore) {
                        bestTree = tree;
                        bestScore = score;
                    }
                }
                if (counter > 0 && counter % stepSize == 0) {
                    System.out.print("*");
                    System.out.flush();
                }
            }
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return null;
        }
        System.out.println();
        System.out.println();
        System.out.println("Best Sum Clade Support: " + bestScore);

        return bestTree;
    }

    private double scoreTree(Tree tree, CladeSystem cladeSystem) {
        return cladeSystem.getSumCladeFrequency(tree, tree.getRoot(), null);
    }

    public static void printTitle() {

        System.out.println("+-----------------------------------------------\\");
        System.out.println("|            TreeSummarizer v1.4 2006           |\\");
        System.out.println("|              MCMC Output analysis             ||");

        String versionString = "BEAST Library: " + version.getVersionString();
        System.out.print("|");
        int n = 47 - versionString.length();
        int n1 = n / 2;
        int n2 = n1 + (n % 2);
        for (int i = 0; i < n1; i++) { System.out.print(" "); }
        System.out.print(versionString);
        for (int i = 0; i < n2; i++) { System.out.print(" "); }
        System.out.println("||");

        System.out.println("|       Andrew Rambaut and Alexei Drummond      ||");
        System.out.println("|              University of Oxford             ||");
        System.out.println("|      http://evolve.zoo.ox.ac.uk/beast/        ||");
        System.out.println("\\-----------------------------------------------\\|");
        System.out.println(" \\-----------------------------------------------\\");
        System.out.println();
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("treesummarizer", "[-burnin <burnin>] <input-file-name> [<output-file-name>]");
        System.out.println();
        System.out.println("  Example: treesummarizer test.trees out.txt");
        System.out.println("  Example: treesummarizer -burnin 100 test.trees out.txt");
        System.out.println();
    }

    private class CladeSystem
    {
        //
        // Public stuff
        //

        /**
         */
        public CladeSystem()
        {
        }

        /** adds all the clades in the tree */
        public void add(Tree tree)
        {
            if (taxonList == null) {
                taxonList = tree;
            }

            // Recurse over the tree and add all the clades (or increment their
            // frequency if already present). The root clade is not added.
            addClades(tree, tree.getRoot(), null);
        }

        public Map getCladeMap() {
            return cladeMap;
        }

        private void addClades(Tree tree, NodeRef node, BitSet bits) {

            if (tree.isExternal(node)) {

                if (taxonList != null) {
                    int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                    bits.set(index);
                } else {
                    bits.set(node.getNumber());
                }
            } else {

                BitSet bits2 = new BitSet();
                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    addClades(tree, node1, bits2);
                }

                addClade(bits2);

                if (bits != null) {
                    bits.or(bits2);
                }
            }
        }

        private void addClade(BitSet bits) {
            Clade clade = (Clade)cladeMap.get(bits);
            if (clade == null) {
                clade = new Clade(bits);
                cladeMap.put(bits, clade);
            }
            clade.setFrequency(clade.getFrequency() + 1);
        }

        public double getSumCladeFrequency(Tree tree, NodeRef node, BitSet bits) {

            double sumCladeFrequency = 0.0;

            if (tree.isExternal(node)) {

                if (taxonList != null) {
                    int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                    bits.set(index);
                } else {
                    bits.set(node.getNumber());
                }
            } else {

                BitSet bits2 = new BitSet();
                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    sumCladeFrequency += getSumCladeFrequency(tree, node1, bits2);
                }

                sumCladeFrequency += getCladeFrequency(bits2);

                if (bits != null) {
                    bits.or(bits2);
                }
            }

            return sumCladeFrequency;
        }

        private double getCladeFrequency(BitSet bits) {
            Clade clade = (Clade)cladeMap.get(bits);
            if (clade == null) {
                return 0.0;
            }
            return clade.getFrequency();
        }

        class Clade {
            public Clade(BitSet bits) {
                this.bits = bits;
                frequency = 0;
            }

            public int getFrequency() {
                return frequency;
            }

            public void setFrequency(int frequency) {
                this.frequency = frequency;
            }

            public boolean equals(Object o) {
                if (this == o) return true;
                if (o == null || getClass() != o.getClass()) return false;

                final Clade clade = (Clade) o;

                if (bits != null ? !bits.equals(clade.bits) : clade.bits != null) return false;

                return true;
            }

            public int hashCode() {
                return (bits != null ? bits.hashCode() : 0);
            }

            int frequency;
            BitSet bits;
        }

        //
        // Private stuff
        //
        TaxonList taxonList = null;

        Map cladeMap = new HashMap();
    }

    int totalTrees = 0;
    int totalTreesUsed = 0;

//Main method
    public static void main(String[] args) throws IOException {

        String inputFileName = null;
        String outputFileName = null;
        String outputTreeFileName = null;

        if (args.length == 0) {
            System.setProperty("com.apple.macos.useScreenMenuBar","true");
            System.setProperty("apple.laf.useScreenMenuBar","true");
            System.setProperty("apple.awt.showGrowBox","true");

            java.net.URL url = LogCombiner.class.getResource("/images/utility.png");
            javax.swing.Icon icon = null;

            if (url != null) {
                icon = new javax.swing.ImageIcon(url);
            }

            String nameString = "TreeSummarizer v1.4";
            String aboutString = "TreeSummarizer v1.4\n" +
                    "©2006 Andrew Rambaut & Alexei Drummond\n" +
                    "University of Oxford";

            ConsoleApplication consoleApp = new ConsoleApplication(nameString, aboutString, icon);

            printTitle();

            TreeSummarizerDialog dialog = new TreeSummarizerDialog(new JFrame());

            if (!dialog.showDialog("TreeSummarizer v1.4")) {
                return;
            }

            int burnin = dialog.getBurnin();

            inputFileName = dialog.getInputFileName();
            if (inputFileName == null) {
                System.err.println("No input file specified");
                return;
            }

            outputFileName = dialog.getOutputFileName();
//            if (outputFileName == null) {
//                System.err.println("No output file specified");
//                return;
//            }

//            outputTreeFileName = dialog.getOutputTreeFileName();
//            if (outputFileName == null) {
//                System.err.println("No output file specified");
//                return;
//            }

            try {
                new TreeSummarizer(burnin, inputFileName, outputFileName);

            } catch (Exception ex) {
                System.err.println("Exception: " + ex.getMessage());
            }
            while (true) {
                Thread.yield();
            }
        }

        TreeSummarizer.printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[] {
                        new Arguments.StringOption("heights", new String[] { "keep", "median", "mean" }, false, "an option of 'keep', 'median' or 'mean'"),
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in'"),
                        new Arguments.Option("help", "option to print this message")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println(ae);
            TreeSummarizer.printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("help")) {
            TreeSummarizer.printUsage(arguments);
            System.exit(0);
        }

        int heights = TreeSummarizer.KEEP_HEIGHTS;
        if (arguments.hasOption("heights")) {
            String value = arguments.getStringOption("heights");
            if (value.equalsIgnoreCase("mean")) {
                heights = TreeSummarizer.MEAN_HEIGHTS;
            } else if (value.equalsIgnoreCase("median")) {
                heights = TreeSummarizer.MEDIAN_HEIGHTS;
            }
        }

        int burnin = -1;
        if (arguments.hasOption("burnin")) {
            burnin = arguments.getIntegerOption("burnin");
        }

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 2) {
            System.err.println("Unknown option: " + args2[2]);
            System.err.println();
            TreeSummarizer.printUsage(arguments);
            System.exit(1);
        }

        if (args2.length == 2) {
            inputFileName = args2[0];
            outputFileName = args2[1];
        } else {
            TreeSummarizer.printUsage(arguments);
            System.exit(1);
        }

        new TreeSummarizer(burnin, inputFileName, outputFileName);

        System.exit(0);
    }
}
