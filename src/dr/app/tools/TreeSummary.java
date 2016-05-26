/*
 * TreeSummary.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.util.Version;

import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * @author Andrew Rambaut
 */
public class TreeSummary {

    private final static Version version = new BeastVersion();

    // Messages to stderr, output to stdout
    private static PrintStream progressStream = System.err;

    /**
     * Burnin can be specified as the number of trees or the number of states
     * (one or other should be zero).
     * @param burninTrees
     * @param burninStates
     * @param posteriorLimit
     * @param inputFileName
     * @param outputFileName
     * @throws java.io.IOException
     */
    public TreeSummary(final int burninTrees,
                       final int burninStates,
                       double posteriorLimit,
                       boolean createSummaryTree,
                       boolean createCladeMap,
                       String inputFileName,
                       String outputFileName
    ) throws IOException {

        this.posteriorLimit = posteriorLimit;

        CladeSystem cladeSystem = new CladeSystem();

        int burnin = -1;

        totalTrees = 10000;
        totalTreesUsed = 0;

        progressStream.println("Reading trees (bar assumes 10,000 trees)...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        cladeSystem = new CladeSystem();
        FileReader fileReader = new FileReader(inputFileName);
        TreeImporter importer = new NexusImporter(fileReader);
        try {
            totalTrees = 0;
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                int state = Integer.MAX_VALUE;

                if (burninStates > 0) {
                    // if burnin has been specified in states, try to parse it out...
                    String name = tree.getId().trim();

                    if (name != null && name.length() > 0 && name.startsWith("STATE_")) {
                        state = Integer.parseInt(name.split("_")[1]);
                    }
                }

                if (totalTrees >= burninTrees && state >= burninStates) {
                    // if either of the two burnin thresholds have been reached...

                    if (burnin < 0) {
                        // if this is the first time this point has been reached,
                        // record the number of trees this represents for future use...
                        burnin = totalTrees;
                    }

                    cladeSystem.add(tree, true);

                    totalTreesUsed += 1;
                }

                if (totalTrees > 0 && totalTrees % stepSize == 0) {
                    progressStream.print("*");
                    progressStream.flush();
                }
                totalTrees++;
            }

        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }
        fileReader.close();
        progressStream.println();
        progressStream.println();

        if (totalTrees < 1) {
            System.err.println("No trees");
            return;
        }
        if (totalTreesUsed <= 1) {
            if (burnin > 0) {
                System.err.println("No trees to use: burnin too high");
                return;
            }
        }
        cladeSystem.calculateCladeCredibilities(totalTreesUsed);

        progressStream.println("Total trees read: " + totalTrees);
        if (burninTrees > 0) {
            progressStream.println("Ignoring first " + burninTrees + " trees" +
                    (burninStates > 0 ? " (" + burninStates + " states)." : "." ));
        } else if (burninStates > 0) {
            progressStream.println("Ignoring first " + burninStates + " states (" + burnin + " trees).");
        }

        progressStream.println("Total unique clades: " + cladeSystem.getCladeMap().keySet().size());
        progressStream.println();

        if (createCladeMap) {

            Map<BitSet, Integer> cladeCountMap = cladeSystem.getCladeCounts();

            System.out.println("No.\tSize\tCred\tMembers");
            int n = 1;
            for (BitSet bits : cladeCountMap.keySet()) {
                System.out.print(n);
                System.out.print("\t");
                System.out.print(bits.cardinality());
                System.out.print("\t");
                System.out.print(cladeSystem.getCladeCredibility(bits));
                System.out.print("\t");
                System.out.println(cladeSystem.getCladeString(bits));
                n++;
            }
            System.out.println();


            progressStream.println("Reading trees...");
            progressStream.println("0              25             50             75            100");
            progressStream.println("|--------------|--------------|--------------|--------------|");

            stepSize = totalTrees / 60;

            fileReader = new FileReader(inputFileName);
            importer = new NexusImporter(fileReader);
            final PrintStream stream = outputFileName != null ?
                    new PrintStream(new FileOutputStream(outputFileName)) :
                    System.out;

            stream.print("Clade");
             n = 1;
            for (BitSet bits : cladeCountMap.keySet()) {
                stream.print("\t");
                stream.print(n);
                n++;
            }
            stream.println();

            stream.print("State");
            for (BitSet bits : cladeCountMap.keySet()) {
                stream.print("\t");
                stream.print(cladeSystem.getCladeCredibility(bits));
            }
            stream.println();

            try {
                totalTrees = 0;
                while (importer.hasTree()) {
                    Tree tree = importer.importNextTree();

                    int state = totalTrees;

                    if (burninStates > 0) {
                        // if burnin has been specified in states, try to parse it out...
                        String name = tree.getId().trim();

                        if (name != null && name.length() > 0 && name.startsWith("STATE_")) {
                            state = Integer.parseInt(name.split("_")[1]);
                        }
                    }

                    if (totalTrees >= burninTrees && state >= burninStates) {
                        // if either of the two burnin thresholds have been reached...

                        if (burnin < 0) {
                            // if this is the first time this point has been reached,
                            // record the number of trees this represents for future use...
                            burnin = totalTrees;
                        }

                        stream.print(state);

                        Set<BitSet> cladeSet = cladeSystem.getCladeSet(tree);

                        for (BitSet bits : cladeCountMap.keySet()) {
                            stream.print("\t");
                            stream.print(cladeSet.contains(bits) ? "1" : "0");
                        }
                        stream.println();

                        totalTreesUsed += 1;
                    }

                    if (totalTrees > 0 && totalTrees % stepSize == 0) {
                        progressStream.print("*");
                        progressStream.flush();
                    }
                    totalTrees++;
                }

            } catch (Importer.ImportException e) {
                System.err.println("Error Parsing Input Tree: " + e.getMessage());
                return;
            }
            fileReader.close();
            stream.close();
            progressStream.println();
            progressStream.println();
        }
        if (createSummaryTree) {
            progressStream.println("Finding summary tree...");

//        MutableTree targetTree = new FlexibleTree(summarizeTrees(burnin, cladeSystem, inputFileName /*, false*/));

            Tree consensusTree = buildConsensusTree(cladeSystem);

            progressStream.println("Writing consensus tree....");

            try {
                final PrintStream stream = outputFileName != null ?
                        new PrintStream(new FileOutputStream(outputFileName)) :
                        System.out;

                new NexusExporter(stream).exportTree(consensusTree);
            } catch (Exception e) {
                System.err.println("Error writing consensus tree file: " + e.getMessage());
                return;
            }
        }

    }

    private Tree buildConsensusTree(CladeSystem cladeSystem) {
        List<BitSet> bitSets = new ArrayList<BitSet>(cladeSystem.getCladeMap().keySet());

        Collections.sort(bitSets, new Comparator<BitSet>() {
            @Override
            public int compare(BitSet b1, BitSet b2) {
                return b1.cardinality() - b2.cardinality();
            }
        });

        SimpleNode root = null;

        for (int i = 0; i < bitSets.size(); i++) {
            BitSet key1 = bitSets.get(i);
            CladeSystem.Clade clade1 = cladeSystem.getCladeMap().get(key1);

            if (key1.cardinality() == 1) {
                Taxon taxon = cladeSystem.getTaxon(key1.nextSetBit(0));
                clade1.node = new SimpleNode();
                clade1.node.setTaxon(taxon);

                clade1.node.setAttribute("clade", clade1);
            }

            if (clade1.credibility >= posteriorLimit) {

                for (int j = i + 1; j < bitSets.size(); j++) {
                    BitSet key2 = bitSets.get(j);
                    if (isSubSet(key1, key2) && cladeSystem.getCladeCredibility(key2) >= posteriorLimit) {
                        // the clades are ordered by size so this is the smallest clade for which
                        // clade1 is a subset and has high credibility.
                        CladeSystem.Clade clade2 = cladeSystem.getCladeMap().get(key2);
                        if (clade2.node == null) {
                            clade2.node = new SimpleNode();
                        }
                        if (clade1.node == null) {
                            throw new RuntimeException("null node");
                        }
                        clade2.node.addChild(clade1.node);

                        clade2.node.setAttribute("credibility", clade2.credibility);
                        clade2.node.setAttribute("clade", clade2);

                        if (key2.cardinality() == cladeSystem.taxonList.getTaxonCount()) {
                            root = clade2.node;
                        }
                        break;
                    }

                }
            }
        }
        SimpleTree tree = new SimpleTree(root);

        return tree;
    }

    private Tree summarizeTrees(int burnin, Tree consensusTree, CladeSystem cladeSystem, String inputFileName)
            throws IOException {

        progressStream.println("Analyzing " + totalTreesUsed + " trees...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        int counter = 0;
        int bestTreeNumber = 0;
        TreeImporter importer = new NexusImporter(new FileReader(inputFileName));

        try {
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (counter >= burnin) {
                    cladeSystem.addSubTrees(tree);
                }
                if (counter > 0 && counter % stepSize == 0) {
                    progressStream.print("*");
                    progressStream.flush();
                }
                counter++;
            }
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return null;
        }

        Tree bestTree = cladeSystem.getBestTree(consensusTree);

        return bestTree;
    }

    private class CladeSystem {

        public CladeSystem() {
        }

        public CladeSystem(TaxonList taxonList) {
            this.taxonList = taxonList;
        }

        /**
         * adds all the clades in the tree
         */
        public void add(Tree tree, boolean includeTips) {
            if (taxonList == null) {
                taxonList = tree;
            }

            // Recurse over the tree and add all the clades (or increment their
            // frequency if already present). The root clade is added too (for
            // annotation purposes).
            addClades(tree, tree.getRoot(), includeTips);
        }

        private BitSet addClades(Tree tree, NodeRef node, boolean includeTips) {

            BitSet bits = new BitSet();

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits.set(index);

                if (includeTips) {
                    addClade(bits);
                }

            } else {

                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    bits.or(addClades(tree, node1, includeTips));
                }

                addClade(bits);
            }

            return bits;
        }

        private void addClade(BitSet bits) {
            Clade clade = cladeMap.get(bits);
            if (clade == null) {
                clade = new Clade(bits);
                cladeMap.put(bits, clade);
            }
            clade.setCount(clade.getCount() + 1);
        }

        public void addSubTrees(Tree tree) {
            addSubTrees(tree, tree.getRoot());
        }

        private BitSet addSubTrees(Tree tree, NodeRef node) {

            BitSet bits = new BitSet();

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits.set(index);

            } else {

                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    BitSet bits2 = addSubTrees(tree, node1);

                    bits.or(bits2);
                }

                Clade clade = cladeMap.get(bits);

                if (clade.credibility >= posteriorLimit) {
                    if (clade.conditionalCladeSystem == null) {
                        clade.conditionalCladeSystem = new CladeSystem(taxonList);
                    }
                    clade.conditionalCladeSystem.addClades(tree, node, false);
                }
            }

            return bits;
        }

        public Map<BitSet, Clade> getCladeMap() {
            return cladeMap;
        }

        public void calculateCladeCredibilities(int totalTreesUsed) {
            for (Clade clade : cladeMap.values()) {

                if (clade.getCount() > totalTreesUsed) {

                    throw new AssertionError("clade.getCount=(" + clade.getCount() +
                            ") should be <= totalTreesUsed = (" + totalTreesUsed + ")");
                }

                clade.setCredibility(((double) clade.getCount()) / (double) totalTreesUsed);

                if (clade.conditionalCladeSystem != null) {
                    clade.conditionalCladeSystem.calculateCladeCredibilities(totalTreesUsed);
                }
            }
        }

        public double getLogCladeCredibility(Tree tree, NodeRef node, BitSet bits) {

            double logCladeCredibility = 0.0;

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits.set(index);
            } else {

                BitSet bits2 = new BitSet();
                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    logCladeCredibility += getLogCladeCredibility(tree, node1, bits2);
                }

                logCladeCredibility += Math.log(getCladeCredibility(bits2));

                if (bits != null) {
                    bits.or(bits2);
                }
            }

            return logCladeCredibility;
        }

        private double getCladeCredibility(BitSet bits) {
            Clade clade = cladeMap.get(bits);
            if (clade == null) {
                return 0.0;
            }
            return clade.getCredibility();
        }

        private int getCladeCount(BitSet bits) {
            Clade clade = cladeMap.get(bits);
            if (clade == null) {
                return 0;
            }
            return clade.getCount();
        }

        public Tree getBestTree(Tree consensusTree) {
            return new SimpleTree(getBestSubTree(consensusTree, consensusTree.getRoot()));
        }

        public SimpleNode getBestSubTree(Tree tree, NodeRef node) {
            SimpleNode subTree = null;

            if (tree.isExternal(node)) {
            } else {

                Clade clade = (Clade) tree.getNodeAttribute(node, "clade");
                if (clade.conditionalCladeSystem != null) {
                    subTree = clade.conditionalCladeSystem.getBestSubTree(tree, node);
                } else {

                    for (int i = 0; i < tree.getChildCount(node); i++) {

                        NodeRef child = tree.getChild(node, i);

                        SimpleNode newChild = getBestSubTree(tree, child);
                    }

                }

            }


            return subTree;
        }

        public Taxon getTaxon(int index) {
            return taxonList.getTaxon(index);
        }

        private BitSet getClades(Tree tree, NodeRef node, boolean includeTips, Set<BitSet> cladeSet) {

            BitSet bits = new BitSet();

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits.set(index);

                if (includeTips) {
                    cladeSet.add(bits);
                }

            } else {

                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    bits.or(getClades(tree, node1, includeTips, cladeSet));
                }

                cladeSet.add(bits);
            }

            return bits;
        }

        public Set<BitSet> getCladeSet(Tree tree) {
            Set<BitSet> cladeSet = new HashSet<BitSet>();
            getClades(tree, tree.getRoot(), false, cladeSet);
            return cladeSet;
        }

        public Map<BitSet, Integer> getCladeCounts() {
            Map<BitSet, Integer> countMap = new HashMap<BitSet, Integer>();

            for (BitSet bits : cladeMap.keySet()) {
                int count = getCladeCount(bits);
                if (count > 1) {
                    countMap.put(bits, count);
                }
            }

            return countMap;
        }

        public String getCladeString(BitSet bits) {
            StringBuilder sb = new StringBuilder("{");
            int index = bits.nextSetBit(0);
            if (index != -1) {
                sb.append(taxonList.getTaxon(index).getId());
                index = bits.nextSetBit(index + 1);
            }
            while (index != -1) {
                sb.append(",");
                sb.append(taxonList.getTaxon(index).getId());
                 index = bits.nextSetBit(index + 1);
            }
            sb.append("}");
            return sb.toString();
        }

        class Clade {
            public Clade(BitSet bits) {
                this.bits = bits;
                count = 0;
                credibility = 0.0;
            }

            public int getCount() {
                return count;
            }

            public void setCount(int count) {
                this.count = count;
            }

            public double getCredibility() {
                return credibility;
            }

            public void setCredibility(double credibility) {
                this.credibility = credibility;
            }

            public boolean equals(Object o) {
                if (this == o) return true;
                if (o == null || getClass() != o.getClass()) return false;

                final Clade clade = (Clade) o;

                return !(bits != null ? !bits.equals(clade.bits) : clade.bits != null);

            }

            public int hashCode() {
                return (bits != null ? bits.hashCode() : 0);
            }

            public String toString() {
                return "clade " + bits.toString();
            }

            int count;
            double credibility;
            BitSet bits;
            SimpleNode node;
            CladeSystem conditionalCladeSystem;
        }

        //
        // Private stuff
        //
        TaxonList taxonList = null;
        Map<BitSet, Clade> cladeMap = new HashMap<BitSet, Clade>();

        Tree targetTree;
    }

    int totalTrees = 0;
    int totalTreesUsed = 0;
    double posteriorLimit = 0.0;

    Set<String> attributeNames = new HashSet<String>();

    public static void printTitle() {
        progressStream.println();
        centreLine("TreeSummary " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("MCMC tree set summarizer", 60);
        centreLine("by", 60);
        centreLine("Andrew Rambaut", 60);
        progressStream.println();
        centreLine("Institute of Evolutionary Biology", 60);
        centreLine("University of Edinburgh", 60);
        centreLine("a.rambaut@ed.ac.uk", 60);
        progressStream.println();
        progressStream.println();
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            progressStream.print(" ");
        }
        progressStream.println(line);
    }


    public static void printUsage(Arguments arguments) {

        arguments.printUsage("treesummary", "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: treesummary test.trees out.txt");
        progressStream.println("  Example: treesummary -burnin 100 -heights mean test.trees out.txt");
        progressStream.println();
    }

    //Main method
    public static void main(String[] args) throws IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        String inputFileName = null;
        String outputFileName = null;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in'"),
                        new Arguments.IntegerOption("burninTrees", "the number of trees to be considered as 'burn-in'"),
                        new Arguments.Option("clademap", "show states of all clades over chain length"),
                        new Arguments.RealOption("limit", "the minimum posterior probability for a subtree to be included"),
                        new Arguments.Option("help", "option to print this message")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            progressStream.println(ae);
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("help")) {
            printUsage(arguments);
            System.exit(0);
        }

        int burninStates = -1;
        int burninTrees = -1;
        if (arguments.hasOption("burnin")) {
            burninStates = arguments.getIntegerOption("burnin");
        }
        if (arguments.hasOption("burninTrees")) {
            burninTrees = arguments.getIntegerOption("burninTrees");
        }

        double posteriorLimit = 0.5;
        if (arguments.hasOption("limit")) {
            posteriorLimit = arguments.getRealOption("limit");
        }

        final String[] args2 = arguments.getLeftoverArguments();

        boolean createCladeMap = false;
        boolean createSummaryTree = true;

        if (arguments.hasOption("clademap")) {
            createCladeMap = true;
            createSummaryTree = false;
        }


        switch (args2.length) {
            case 2:
                outputFileName = args2[1];
                // fall to
            case 1:
                inputFileName = args2[0];
                break;
            default: {
                System.err.println("Unknown option: " + args2[2]);
                System.err.println();
                printUsage(arguments);
                System.exit(1);
            }
        }

        new TreeSummary(burninTrees, burninStates, posteriorLimit, createSummaryTree, createCladeMap, inputFileName, outputFileName);

        System.exit(0);
    }

    /**
     * Is x a subset of y?
     * @param x
     * @param y
     * @return
     */
    static boolean isSubSet(BitSet x, BitSet y) {
        y = (BitSet) y.clone();
        y.and(x);
        return y.equals(x);
    }

}

