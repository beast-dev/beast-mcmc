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
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;
import dr.util.Version;
import org.virion.jam.console.ConsoleApplication;

import javax.swing.*;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class TreeAnnotator {

    private final static Version version = new BeastVersion();

    public final static int MAX_CLADE_CREDIBILITY = 0;
    public final static int USER_TARGET_TREE = 1;

    public final static int KEEP_HEIGHTS = 0;
    public final static int MEAN_HEIGHTS = 1;
    public final static int MEDIAN_HEIGHTS = 2;

    public TreeAnnotator(int burnin, int heightsOption, double posteriorLimit, int targetOption, String targetTreeFileName, String inputFileName, String outputFileName) throws IOException {

        this.posteriorLimit = posteriorLimit;

        System.out.println("Reading trees...");

        CladeSystem cladeSystem = new CladeSystem();
        cladeSystem.setAttributeNames(attributeNames);

        TreeImporter importer = new NexusImporter(new FileReader(inputFileName));
        try {
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (totalTrees >= burnin) {
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
        if (burnin > 0) {
            System.out.println("Ignoring first" + burnin + " trees.");
        }

        MutableTree targetTree;

        if (targetOption == USER_TARGET_TREE) {
            if (targetTreeFileName != null) {
                System.out.println("Reading user specified target tree, " + targetTreeFileName);

                importer = new NexusImporter(new FileReader(targetTreeFileName));
                try {
                    targetTree = new FlexibleTree(importer.importNextTree());
                } catch (Importer.ImportException e) {
                    System.err.println("Error Parsing Target Tree: " + e.getMessage());
                    return;
                }
            } else {
                System.err.println("No user target tree specified.");
                return;
            }
        } else {
            System.out.println("Finding maximum clade credibility tree...");
            targetTree = new FlexibleTree(summarizeTrees( burnin, cladeSystem, inputFileName));
        }

        System.out.println("Annotating target tree...");
        cladeSystem.annotateTree(targetTree, targetTree.getRoot(), null, heightsOption);

        System.out.println("Writing annotated tree....");
        if (outputFileName != null) {
            NexusExporter exporter = new NexusExporter(new PrintStream(new FileOutputStream(outputFileName)));
            exporter.exportTree(targetTree);
        } else {
            NexusExporter exporter = new NexusExporter(System.out);
            exporter.exportTree(targetTree);
        }
    }

    private Tree summarizeTrees(int burnin, CladeSystem cladeSystem, String inputFileName) throws IOException {

        Tree bestTree = null;
        double bestScore = 0.0;

        System.out.println("Analyzing " + totalTreesUsed + " trees...");
        System.out.println("0              25             50             75            100");
        System.out.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        int counter = 0;
        TreeImporter importer = new NexusImporter(new FileReader(inputFileName));
        try {
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (counter >= burnin) {
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
                counter ++;
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
            // frequency if already present). The root clade is added too (for
            // annotation purposes).
            addClades(tree, tree.getRoot(), null);
        }

        public Map getCladeMap() {
            return cladeMap;
        }

        public Clade getClade(NodeRef node) {
            return null;
        }

        private void addClades(Tree tree, NodeRef node, BitSet bits) {

            BitSet bits2 = new BitSet();

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits2.set(index);

            } else {

                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    addClades(tree, node1, bits2);
                }
            }

            addClade(bits2, tree, node);

            if (bits != null) {
                bits.or(bits2);
            }
        }

        private void addClade(BitSet bits, Tree tree, NodeRef node) {
            Clade clade = (Clade)cladeMap.get(bits);
            if (clade == null) {
                clade = new Clade(bits);
                cladeMap.put(bits, clade);
            }
            clade.setFrequency(clade.getFrequency() + 1);

            if (attributeNames != null) {
                if (clade.attributeLists == null) {
                    clade.attributeLists = new List[attributeNames.length];
                    for (int i = 0; i < attributeNames.length; i++) {
                        clade.attributeLists[i] = new ArrayList();
                    }
                }

                for (int i = 0; i < attributeNames.length; i++) {
                    Object value;
                    if (attributeNames[i].equals("height")) {
                        value = new Double(tree.getNodeHeight(node));
                    } else if (attributeNames[i].equals("length")) {
                        value = new Double(tree.getBranchLength(node));
                    } else {
                        value = tree.getNodeAttribute(node, attributeNames[i]);
                    }

                    //if (value == null) {
                    //    System.out.println("attribute " + attributeNames[i] + " is null.");
                    //}

                    if (value != null) {
                        clade.attributeLists[i].add(value);
                    }
                }
            }
        }

        public double getSumCladeFrequency(Tree tree, NodeRef node, BitSet bits) {

            double sumCladeFrequency = 0.0;

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits.set(index);
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

        public void annotateTree(MutableTree tree, NodeRef node, BitSet bits, int heightsOption) {

            BitSet bits2 = new BitSet();

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits2.set(index);

                annotateNode(tree, node, bits2, true, heightsOption);
            } else {

                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    annotateTree(tree, node1, bits2, heightsOption);
                }

                annotateNode(tree, node, bits2, false, heightsOption);
            }

            if (bits != null) {
                bits.or(bits2);
            }
        }

        private void annotateNode(MutableTree tree, NodeRef node, BitSet bits, boolean isTip, int heightsOption) {
            Clade clade = (Clade)cladeMap.get(bits);
            if (clade == null) {
                throw new RuntimeException("Clade missing");
            }

            boolean filter = false;
            if (!isTip) {
                double posterior = ((double)clade.frequency) / totalTreesUsed;
                tree.setNodeAttribute(node, "posterior", new Double(posterior));
                if (posterior < posteriorLimit) {
                    filter = true;
                }
            }

            for (int i = 0; i < clade.attributeLists.length; i++) {
                boolean isHeight = attributeNames[i].equals("height");
                boolean isBoolean = attributeTypes[i].equals(Boolean.class);

                double[] values = new double[clade.attributeLists[i].size()];
                if (values.length > 0) {
                    double minValue = Double.MAX_VALUE;
                    double maxValue = -Double.MAX_VALUE;
                    for (int j = 0; j < clade.attributeLists[i].size(); j++) {
                        values[j] = ((Number)clade.attributeLists[i].get(j)).doubleValue();
                        if (values[j] < minValue) minValue = values[j];
                        if (values[j] > maxValue) maxValue = values[j];
                    }
                    if (isHeight) {
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

                    if (!filter) {
                        annotateMeanAttribute(tree, node, attributeNames[i], values);
                        if (!isBoolean && minValue < maxValue) {
                            // Basically, if it is a boolean (0, 1) then we don't need the distribution information
                            // Likewise if it doesn't vary.
                            annotateMedianAttribute(tree, node, attributeNames[i] + "_median", values);
                            annotateHPDAttribute(tree, node, attributeNames[i] + "_95%_HPD", 0.95, values);
                            annotateRangeAttribute(tree, node, attributeNames[i] + "_range", values);
                        }
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
    double posteriorLimit = 0.0;
    String[] attributeNames = new String[] { "height", "rate", "changed" };
    Class[] attributeTypes = new Class[] { Double.class, Double.class, Boolean.class};
    TaxonList taxa = null;

    public static void printTitle() {
        System.out.println();
        centreLine("TreeAnnotator " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("MCMC Output analysis", 60);
        centreLine("by", 60);
        centreLine("Andrew Rambaut and Alexei J. Drummond", 60);
        System.out.println();
        centreLine("Institute of Evolutionary Biology", 60);
        centreLine("University of Edinburgh", 60);
        centreLine("a.rambaut@ed.ac.uk", 60);
        System.out.println();
        centreLine("Department of Computer Science", 60);
        centreLine("University of Auckland", 60);
        centreLine("alexei@cs.auckland.ac.nz", 60);
        System.out.println();
        System.out.println();
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) { System.out.print(" "); }
        System.out.println(line);
    }


    public static void printUsage(Arguments arguments) {

        arguments.printUsage("treeannotator", "[-burnin <burnin>] [-heights <height_option>] [-limit <posterior_limit>] [-target <target-tree-file-name>] <input-file-name> [<output-file-name>]");
        System.out.println();
        System.out.println("  Example: treeannotator test.trees out.txt");
        System.out.println("  Example: treeannotator -burnin 100 -heights mean test.trees out.txt");
        System.out.println("  Example: treeannotator -burnin 100 -target map.tree test.trees out.txt");
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

	        final String versionString = version.getVersionString();
            String nameString = "TreeAnnotator " + versionString;
            String aboutString = "<html><center><p>" + versionString + ", " + version.getDateString() +"</p>" +
                    "<p>by<br>" +
                    "Andrew Rambaut and Alexei J. Drummond</p>" +
                    "<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
                    "<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
                    "<p>Department of Computer Science, University of Auckland<br>" +
                    "<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +
                    "<p>Part of the BEAST package:<br>" +
                    "<a href=\"http://beast.bio.ed.ac.uk/\">http://beast.bio.ed.ac.uk/</a></p>" +
                    "</center></html>";

            ConsoleApplication consoleApp = new ConsoleApplication(nameString, aboutString, icon, true);

            printTitle();

            TreeAnnotatorDialog dialog = new TreeAnnotatorDialog(new JFrame());

            if (!dialog.showDialog("TreeAnnotator " + versionString)) {
                return;
            }

            int burnin = dialog.getBurnin();
            double posteriorLimit = dialog.getPosteriorLimit();
            int targetOption = dialog.getTargetOption();
            int heightsOption = dialog.getHeightsOption();

            targetTreeFileName = dialog.getTargetFileName();
            if (targetOption == USER_TARGET_TREE && targetTreeFileName == null) {
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
                new TreeAnnotator(burnin, heightsOption, posteriorLimit, targetOption, targetTreeFileName, inputFileName, outputFileName);

            } catch (Exception ex) {
                System.err.println("Exception: " + ex.getMessage());
            }

            System.out.println("Finished - Quit program to exit.");
            while (true) {
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[] {
                        //new Arguments.StringOption("target", new String[] { "maxclade", "maxtree" }, false, "an option of 'maxclade' or 'maxtree'"),
                        new Arguments.StringOption("heights", new String[] { "keep", "median", "mean" }, false, "an option of 'keep', 'median' or 'mean'"),
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in'"),
                        new Arguments.RealOption("limit", "the minimum posterior probability for a node to be annoated"),
                        new Arguments.StringOption("target", "target_file_name", "specifies a user target tree to be annotated"),
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

        double posteriorLimit = 0.0;
        if (arguments.hasOption("limit")) {
            posteriorLimit = arguments.getRealOption("limit");
        }

        int target = MAX_CLADE_CREDIBILITY;
        if (arguments.hasOption("target")) {
            target = USER_TARGET_TREE;
            targetTreeFileName = arguments.getStringOption("target");
        }


        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 2) {
            System.err.println("Unknown option: " + args2[2]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (args2.length == 2) {
            targetTreeFileName = null;
            inputFileName = args2[0];
            outputFileName = args2[1];
        } else {
            printUsage(arguments);
            System.exit(1);
        }

        new TreeAnnotator(burnin, heights, posteriorLimit, target, targetTreeFileName, inputFileName, outputFileName);

        System.exit(0);
    }

}

