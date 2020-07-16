/*
 * TreePruner.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.arg.ARGModel;
import dr.util.Version;

import java.io.*;
import java.util.*;

/**
 * @author Philippe Lemey
 * @author Marc Suchard
 * @author Andrew Rambaut
 */
public class TaxaPicker {

    private final static Version version = new BeastVersion();

    // Messages to stderr, output to stdout
    private static PrintStream progressStream = System.err;

    private static final String[] falseTrue = {"false", "true"};
    private static final String NAMECONTENT = "nameContent";

    private TaxaPicker(String inputFileName,
                       String outputFileName,
                       String[] taxaToPrune,
                       boolean basedOnNameContent,
                       int degree,
                       double cutOff
                       ) throws IOException {


        this.degree = degree;
        this.relatives = new HashMap<>();

        List<Tree> trees = new ArrayList<>();

        readTrees(trees, inputFileName);

        List<Taxon> taxa = getTaxaToPrune(trees.get(0), taxaToPrune, basedOnNameContent);

        processTrees(trees, taxa);

        for (String remove : taxaToPrune) {
            relatives.remove(remove);
        }

        for (String name : relatives.keySet()) {
            double probability = (double) relatives.get(name) / ((double) totalTrees);
            if (probability < cutOff) {
                relatives.remove(name);
            }
        }

        writeOutputFile(trees, outputFileName);
    }

    private final int degree;
    private final Map<String, Integer> relatives;

    private void readTrees(List<Tree> trees, String inputFileName) throws IOException {

        progressStream.println("Reading trees (bar assumes 10,000 trees)...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        long stepSize = 10000 / 60;

        FileReader fileReader = new FileReader(inputFileName);
        TreeImporter importer = new NexusImporter(fileReader, false);

        try {
            totalTrees = 0;
            while (importer.hasTree()) {

                Tree tree = importer.importNextTree();
                if (trees == null) {
                    trees = new ArrayList<>();
                }
                trees.add(tree);

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

        progressStream.println("Total trees read: " + totalTrees);
    }

//    private List<Taxon> getTaxaToPrune(Tree tree, String[] names) {
//
//        List<Taxon> taxa = new ArrayList<>();
//        if (names != null) {
//            for (String name : names) {
//
//                int taxonId = tree.getTaxonIndex(name);
//                if (taxonId == -1) {
//                    throw new RuntimeException("Unable to find taxon '" + name + "'.");
//                }
//
//                taxa.add(tree.getTaxon(taxonId));
//            }
//        }
//        return taxa;
//    }
    private List<Taxon> getTaxaToPrune(Tree tree, String[] names, boolean basedOnContent) {

        List<Taxon> taxa = new ArrayList<>();
        if (names != null) {
            for (String name : names) {

                if(!basedOnContent){
                    int taxonId = tree.getTaxonIndex(name);
                    if (taxonId == -1) {
                        throw new RuntimeException("Unable to find taxon '" + name + "'.");
                    }
                    taxa.add(tree.getTaxon(taxonId));
                } else {
                    int counter = 0;
                    for(int i = 0; i < tree.getTaxonCount(); i++) {
                        Taxon taxon = tree.getTaxon(i);
                        String taxonName = taxon.toString();
                        if (taxonName.contains(name)) {
                            taxa.add(taxon);
                            counter ++;
                        }
                    }
                    if (counter == 0){
                        throw new RuntimeException("Unable to find taxon with a name containing '" + name + "'.");
                    }
                }
            }
        }
        return taxa;
    }

    private void processTrees(List<Tree> trees, List<Taxon> taxa) {
        for (Tree tree : trees) {
            processOneTree(tree, taxa);
        }
    }

    private void processOneTree(Tree tree, List<Taxon> taxa) {
        for (Taxon taxon : taxa) {
            processOneTreeForOneTaxon((FlexibleTree)tree, taxon);
        }
    }

    private void processOneTreeForOneTaxon(FlexibleTree tree, Taxon taxon) {
        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            NodeRef tip = tree.getExternalNode(i);
            if (tree.getNodeTaxon(tip) == taxon) {
                processOneTip(tree, tip);
            }
        }
    }

    private void processOneTip(FlexibleTree tree, NodeRef tip) {
        NodeRef parent = tree.getParent(tip);

        NodeRef ancestor = tip;

        int currentDegree = 0;
        while (ancestor != tree.getRoot() && currentDegree < degree) {
            ancestor = tree.getParent(ancestor);
            ++currentDegree;
        }

        recursivelyAddTipsNames(tree, ancestor, 0);
    }

    private void recursivelyAddTipsNames(Tree tree, NodeRef node, int currentDegree) {

        if (!tree.isExternal(node)) {
            if (currentDegree < degree) {
                recursivelyAddTipsNames(tree, tree.getChild(node, 0), currentDegree + 1);
                recursivelyAddTipsNames(tree, tree.getChild(node, 1), currentDegree + 1);
            }
        } else {
            Taxon taxon = tree.getNodeTaxon(node);
            addName(taxon.getId());
        }
    }

    private void addName(String name) {
        int count = 1;
        if (relatives.containsKey(name)) {
            count += relatives.get(name);
        }
        relatives.put(name, count);
    }

    private NodeRef getSibling(FlexibleTree tree, NodeRef parent, NodeRef node) {
        NodeRef sibling = tree.getChild(parent, 0);
        if (sibling == node) {
            sibling = tree.getChild(parent, 1);
        }
        return sibling;
    }

    int totalTrees = 0;

    private void writeOutputFile(List<Tree> trees, String outputFileName) {

        PrintStream ps = null;

        if (outputFileName == null) {
            ps = progressStream;
        } else {
            try {
                ps = new PrintStream(new File(outputFileName));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }

        for (String name : relatives.keySet()) {
            ps.println(name);
        }

        if (ps != null) {
            ps.close();
        }


    }

    private String[] getTreeNames(List<Tree> trees) {
        List<String> names = new ArrayList<>();
        for (Tree tree : trees) {
            names.add(tree.getId());
        }
        return names.toArray(new String[0]);
    }

    public static void printTitle() {
        progressStream.println();
        centreLine("TreePruner " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("Tree pruning tool", 60);
        centreLine("by", 60);
        centreLine("Philippe Lemey, Andrew Rambaut and Marc Suchard", 60);
//        progressStream.println();
//        centreLine("Institute of Evolutionary Biology", 60);
//        centreLine("University of Edinburgh", 60);
//        centreLine("a.rambaut@ed.ac.uk", 60);
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

        arguments.printUsage("treepruner", "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: treepruner -taxaToPrune taxon1,taxon2 input.trees output.trees");
        progressStream.println();
    }

    private static String[] parseTextList(String fileName) {

        List<String> names = new ArrayList<>();

        try {
            BufferedReader br = new BufferedReader(new FileReader(fileName));

            String st;
            while ((st = br.readLine()) != null) {
                names.add(st.trim());
            }

        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(-1);
        }

        return names.toArray(new String[0]);
    }

    //Main method
    public static void main(String[] args) throws IOException {

        String inputFileName = null;
        String outputFileName = null;
        String[] taxaToPrune = null;
        boolean basedOnNameContent = false;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.StringOption("taxaToFind", "list","a list of taxon names to find"),
                        new Arguments.IntegerOption("degree", "degree relatives"),
                        new Arguments.RealOption("cutOff", "cut-off probability"),
                        new Arguments.Option("asFile", "Boolean if taxaToFind is a file"),
                        new Arguments.StringOption(NAMECONTENT, falseTrue, false,
                                "add true noise [default = true])"),
                        new Arguments.Option("help", "option to print this message"),
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

        if (arguments.hasOption("taxaToFind")) {
            if (arguments.hasOption("asFile")) {
                taxaToPrune = parseTextList(arguments.getStringOption("taxaToFind"));
            } else {
                taxaToPrune = Branch2dRateToGrid.parseVariableLengthStringArray(arguments.getStringOption("taxaToFind"));
            }
        }

        double cutOff = 0.8;
        if (arguments.hasOption("cutOff")) {
            cutOff = arguments.getRealOption("cutOff");
        }

        String nameContentString = arguments.getStringOption(NAMECONTENT);
        if (nameContentString != null && nameContentString.compareToIgnoreCase("true") == 0)
            basedOnNameContent = true;

        int degree = arguments.getIntegerOption("degree");

        final String[] args2 = arguments.getLeftoverArguments();

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

        new TaxaPicker(inputFileName, outputFileName, taxaToPrune, basedOnNameContent, degree, cutOff);

        System.exit(0);
    }
}