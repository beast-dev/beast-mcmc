/*
 * TreeAnnotator.java
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

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.evolution.io.*;
import dr.evolution.tree.*;
import dr.evolution.util.TaxonList;
import dr.stats.DiscreteStatistics;
import dr.util.Version;
import dr.util.HeapSort;
import org.virion.jam.console.ConsoleApplication;

import javax.swing.*;
import java.io.*;
import java.util.*;

public class TreeAnnotator {

    private final static Version version = new BeastVersion();

    public final static int KEEP_HEIGHTS = 0;
    public final static int MEAN_HEIGHTS = 1;
    public final static int MEDIAN_HEIGHTS = 2;

    public TreeAnnotator(int burnin, int heightsOption, String targetTreeFileName, String inputFileName, String outputFileName) throws IOException {

        MutableTree targetTree = null;
        if (targetTreeFileName != null) {
            NexusImporter importer = new NexusImporter(new FileReader(targetTreeFileName));
            try {
                targetTree = new FlexibleTree(importer.importNextTree());
            } catch (Importer.ImportException e) {
                System.err.println("Error Parsing Target Tree: " + e.getMessage());
                return;
            }
        }

        CladeSystem targetClades = new CladeSystem(targetTree);
        targetClades.setAttributeNames(attributeNames);

        System.out.println("Reading trees...");

        int fileTotalTrees = 0;
        TreeImporter importer = new NexusImporter(new FileReader(inputFileName));
        try {
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (fileTotalTrees >= burnin) {
                    targetClades.add(tree);

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

        System.out.println("Annotating target tree...");
        annotateTree(targetClades, targetTree, heightsOption);

        System.out.println("Writing annotated tree....");
        if (outputFileName != null) {
            NexusExporter exporter = new NexusExporter(new PrintStream(new FileOutputStream(outputFileName)));
            exporter.exportTree(targetTree);
        } else {
            NexusExporter exporter = new NexusExporter(System.out);
            exporter.exportTree(targetTree);
        }
        System.out.println("Finished.");
    }

    private void annotateTree(CladeSystem cladeSystem, MutableTree targetTree, int heightsOption) {
        Map clades = cladeSystem.getCladeMap();
        Iterator iter = clades.keySet().iterator();
        while (iter.hasNext()) {
            Object key = iter.next();
            CladeSystem.Clade clade = (CladeSystem.Clade)clades.get(key);

            annotateNode(targetTree, clade, heightsOption);
        }
    }

    private void annotateNode(MutableTree tree, CladeSystem.Clade clade, int heightsOption) {
        NodeRef node = clade.node;
        double posterior = ((double)clade.frequency) / totalTreesUsed;
        tree.setNodeAttribute(node, "posterior", new Double(posterior));

        if (posterior >= 0.5) {
            for (int i = 0; i < clade.attributeLists.length; i++) {
                double[] values = new double[clade.attributeLists[i].size()];
                if (values.length > 0) {
                    for (int j = 0; j < clade.attributeLists[i].size(); j++) {
                        values[j] = ((Double)clade.attributeLists[i].get(j)).doubleValue();
                    }
                    if (attributeNames[i].equals("height")) {
                        if (heightsOption == MEAN_HEIGHTS) {
                            double mean = DiscreteStatistics.mean(values);
                            tree.setNodeHeight(node, mean);
                        } else if (heightsOption == MEDIAN_HEIGHTS) {
                            double median = DiscreteStatistics.median(values);
                            tree.setNodeHeight(node, median);
                        } else {
                            // keep the existing height
                        }
                    }
                    annotateMeanAttribute(tree, node, attributeNames[i] + "_mean", values);
                    annotateMedianAttribute(tree, node, attributeNames[i] + "_median", values);
                    annotateHPDAttribute(tree, node, attributeNames[i] + "_95%_HPD", 0.95, values);
                    annotateQuantileAttribute(tree, node, attributeNames[i] + "_95%_quantiles", 0.95, values);
                    annotateRangeAttribute(tree, node, attributeNames[i] + "_range", values);
                }
            }
        }
    }

    private void annotateMeanAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
        double mean = DiscreteStatistics.mean(values);
        tree.setNodeAttribute(node, label, new Double(mean));
    }

    private void annotateMedianAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
        double median = DiscreteStatistics.median(values);
        tree.setNodeAttribute(node, label, new Double(median));

    }

    private void annotateRangeAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
        double min = DiscreteStatistics.min(values);
        double max = DiscreteStatistics.max(values);
        tree.setNodeAttribute(node, label, new Object[] { new Double(min), new Double(max) });
    }

    private void annotateQuantileAttribute(MutableTree tree, NodeRef node, String label, double quantile, double[] values) {
        double lower = DiscreteStatistics.quantile(1.0 - (quantile / 2), values);
        double upper = DiscreteStatistics.quantile(quantile / 2, values);
        tree.setNodeAttribute(node, label, new Object[] { new Double(lower), new Double(upper) });
    }

    private void annotateHPDAttribute(MutableTree tree, NodeRef node, String label, double hpd, double[] values) {
        int[] indices = new int[values.length];
        HeapSort.sort(values, indices);

        double minRange = Double.MAX_VALUE;
        int hpdIndex = 0;

        int diff = (int)Math.round(hpd * (double)values.length);
        for (int i =0; i <= (values.length - diff); i++) {
            double minValue = values[indices[i]];
            double maxValue = values[indices[i+diff-1]];
            double range = Math.abs(maxValue - minValue);
            if (range < minRange) {
                minRange = range;
                hpdIndex = i;
            }
        }
        double lower = values[indices[hpdIndex]];
        double upper = values[indices[hpdIndex+diff-1]];
        tree.setNodeAttribute(node, label, new Object[] { new Double(lower), new Double(upper) });
    }

    private class CladeSystem
    {
        //
        // Public stuff
        //

        /**
         * @param tree
         */
        public CladeSystem(Tree tree)
        {
            this.taxonList = tree;
            addClades(tree, tree.getRoot(), null, false);
        }

        public void setAttributeNames(String[] attributeNames) {
            this.attributeNames = attributeNames;
        }

        /** adds all the clades in the tree */
        public void add(Tree tree)
        {
            if (taxonList == null) {
                taxonList = tree;
            }

            // Recurse over the tree and add all the clades (or increment their
            // frequency if already present). The root clade is not added.
            addClades(tree, tree.getRoot(), null, true);
        }

        public Map getCladeMap() {
            return cladeMap;
        }

        private void addClades(Tree tree, NodeRef node, BitSet bits, boolean collectAttributes) {

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

                    addClades(tree, node1, bits2, collectAttributes);
                }

                if (collectAttributes) {
                    addClade(bits2, tree, node);
                } else {
                    putClade(bits2, tree, node);
                }

                if (bits != null) {
                    bits.or(bits2);
                }
            }
        }

        private void putClade(BitSet bits, Tree tree, NodeRef node) {
            Clade clade = (Clade)cladeMap.get(bits);
            if (clade == null) {
                clade = new Clade(node, bits);
                cladeMap.put(bits, clade);
            } else {
                throw new IllegalArgumentException("Clade already put in CladeSystem");
            }
        }

        private void addClade(BitSet bits, Tree tree, NodeRef node) {
            Clade clade = (Clade)cladeMap.get(bits);
            if (clade == null) {
                return;
            }
            clade.setFrequency(clade.getFrequency() + 1);

            if (clade.attributeLists == null) {
                clade.attributeLists = new List[attributeNames.length];
                for (int i = 0; i < attributeNames.length; i++) {
                    clade.attributeLists[i] = new ArrayList();
                }
            }

            clade.attributeLists[0].add(new Double(tree.getNodeHeight(node)));
            for (int i = 0; i < attributeNames.length; i++) {
                Object value;
                if (attributeNames[i].equals("height")) {
                    value = new Double(tree.getNodeHeight(node));
                } else if (attributeNames[i].equals("length")) {
                    value = new Double(tree.getBranchLength(node));
                } else {
                    value = tree.getAttribute(attributeNames[i]);
                }
                if (value != null) {
                    clade.attributeLists[i].add(value);
                }
            }
        }

        class Clade {
            public Clade(NodeRef node, BitSet bits) {
                this.node = node;
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
            NodeRef node;
            BitSet bits;
            List[] attributeLists = null;
        }

        //
        // Private stuff
        //
        TaxonList taxonList = null;

        Map cladeMap = new HashMap();
        String[] attributeNames = new String[0];
    }

    int totalTrees = 0;
    int totalTreesUsed = 0;
    String[] attributeNames = new String[] { "height", "length", "rate" };
    TaxonList taxa = null;

    public static void printTitle() {

        System.out.println("+-----------------------------------------------\\");
        System.out.println("|            TreeAnnotator v1.4 2006            |\\");
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

        arguments.printUsage("treeannotator", "[-burnin <burnin>] <target-tree-file-name> <input-file-name> [<output-file-name>]");
        System.out.println();
        System.out.println("  Example: treeannotator test.trees out.txt");
        System.out.println("  Example: treeannotator -burnin 10000 -target map.tree test.trees out.txt");
        System.out.println();
    }

//Main method
    public static void main(String[] args) throws IOException {

        String targetTreeFileName = null;
        String inputFileName = null;
        String outputFileName = null;

        if (args.length == 0) {
            System.setProperty("com.apple.macos.useScreenMenuBar","true");
            System.setProperty("apple.laf.useScreenMenuBar","true");
            System.setProperty("apple.awt.showGrowBox","true");

            java.net.URL url = LogCombiner.class.getResource("/images/utility.png");
            javax.swing.Icon icon = null;

            if (url != null) {
                icon = new javax.swing.ImageIcon(url);
            }

            String nameString = "TreeAnnotator v1.4";
            String aboutString = "TreeAnnotator v1.4\n" +
                    "©2006 Andrew Rambaut & Alexei Drummond\n" +
                    "University of Oxford";

            ConsoleApplication consoleApp = new ConsoleApplication(nameString, aboutString, icon);

            printTitle();

            TreeAnnotatorDialog dialog = new TreeAnnotatorDialog(new JFrame());

            if (!dialog.showDialog("TreeAnnotator v1.4")) {
                return;
            }

            int burnin = dialog.getBurnin();
            int heightsOption = dialog.getHeightsOption();

            targetTreeFileName = dialog.getTargetFileName();
            if (targetTreeFileName == null) {
                System.err.println("No target file specified");
                return;
            }

            inputFileName = dialog.getInputFileName();
            if (inputFileName == null) {
                System.err.println("No input file specified");
                return;
            }

            outputFileName = dialog.getOutputFileName();
            if (outputFileName == null) {
                System.err.println("No output file specified");
                return;
            }

            try {
                new TreeAnnotator(burnin, heightsOption, targetTreeFileName, inputFileName, outputFileName);

            } catch (Exception ex) {
                System.err.println("Exception: " + ex.getMessage());
            }
            while (true) {
                Thread.yield();
            }
        }

        printTitle();

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
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("help")) {
            printUsage(arguments);
            System.exit(0);
        }

        int heights = KEEP_HEIGHTS;
        if (arguments.hasOption("heights")) {
            String value = arguments.getStringOption("heights");
            if (value.equalsIgnoreCase("mean")) {
                heights = MEAN_HEIGHTS;
            } else if (value.equalsIgnoreCase("median")) {
                heights = MEDIAN_HEIGHTS;
            }
        }

        int burnin = -1;
        if (arguments.hasOption("burnin")) {
            burnin = arguments.getIntegerOption("burnin");
        }

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 3) {
            System.err.println("Unknown option: " + args2[3]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (args2.length == 3) {
            targetTreeFileName = args2[0];
            inputFileName = args2[1];
            outputFileName = args2[2];
        } else {
            printUsage(arguments);
            System.exit(1);
        }

        new TreeAnnotator(burnin, heights, targetTreeFileName, inputFileName, outputFileName);

        System.exit(0);
    }
}

