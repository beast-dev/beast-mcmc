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
import dr.evolution.io.*;
import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.util.Version;

import java.io.*;
import java.util.*;

/**
 * @author Philippe Lemey
 * @author Marc Suchard
 * @author Andrew Rambaut
 */
public class TreePruner {

    private final static Version version = new BeastVersion();

    // Messages to stderr, output to stdout
    private static PrintStream progressStream = System.err;

    private static final String[] falseTrue = {"false", "true"};
    private static final String NAMECONTENT = "nameContent";
    private static final String PRUNEDNODES = "prunedNodes";

    private TreePruner(String inputFileName,
                       String outputFileName,
                       String[] taxaToPrune,
                       boolean basedOnNameContent
                       ) throws IOException {

        List<Tree> trees = new ArrayList<>();

        readTrees(trees, inputFileName);

        List<Taxon> taxa = getTaxaToPrune(trees.get(0), taxaToPrune, basedOnNameContent);

        processTrees(trees, taxa);

        writeOutputFile(trees, outputFileName);
    }

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

    private List<Taxon> getTaxaToPrune(Tree tree, String[] names, boolean basedOnContent) {
        System.out.println("pruning taxa:");
        List<Taxon> taxa = new ArrayList<>();
        if (names != null) {
            for (String name : names) {

                if(!basedOnContent){
                    int taxonId = tree.getTaxonIndex(name);
                    if (taxonId == -1) {
                        throw new RuntimeException("Unable to find taxon '" + name + "'.");
                    }
                    System.out.println(name);
                    taxa.add(tree.getTaxon(taxonId));
                } else {
                    int counter = 0;
                    for(int i = 0; i < tree.getTaxonCount(); i++) {
                        Taxon taxon = tree.getTaxon(i);
                        String taxonName = taxon.toString();
                        if (taxonName.contains(name)) {
                            taxa.add(taxon);
                            System.out.println(taxonName);
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
            setPrunedNodeAnnotation(tree);
            processOneTree(tree, taxa);
        }
    }

    private void setPrunedNodeAnnotation(Tree tree) {
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            FlexibleNode node = (FlexibleNode) tree.getNode(i);
            node.setAttribute(PRUNEDNODES, "0");
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
//                tree.removeTaxon(taxon);
            }
        }
    }

    private void processOneTip(FlexibleTree tree, NodeRef tip) {
        NodeRef parent = tree.getParent(tip);

        if (parent == tree.getRoot()) {
            tree.beginTreeEdit();
            FlexibleNode sibling = (FlexibleNode) getSibling(tree, parent, tip);
            sibling.setParent(null);

            tree.setRoot(sibling);
            // Annotate pruned nodes
            sibling.setAttribute(PRUNEDNODES, Integer.toString(((Integer.valueOf((String) sibling.getAttribute(PRUNEDNODES)) + 1))));
            tree.endTreeEdit();
        } else {

            NodeRef grandParent = tree.getParent(parent);
            NodeRef sibling = getSibling(tree, parent, tip);

            FlexibleNode grandParentNode = (FlexibleNode)grandParent;
            FlexibleNode parentNode = (FlexibleNode) parent;
            FlexibleNode siblingNode = (FlexibleNode) sibling;
//            FlexibleNode tipNode = (FlexibleNode) tip;

            // Remove from topology
            siblingNode.setParent(grandParentNode);
            grandParentNode.removeChild(parentNode);
            grandParentNode.addChild(siblingNode);

            // Adjust branch lengths
            siblingNode.setLength(parentNode.getLength() + siblingNode.getLength());

            // Annotate pruned nodes
            siblingNode.setAttribute(PRUNEDNODES, Integer.toString(((Integer.valueOf((String) siblingNode.getAttribute(PRUNEDNODES)) + 1))));

            // Combine traits
            // TODO
        }
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

        NexusExporter exporter = new NexusExporter(ps);

        if (trees.size() > 0) {
            exporter.exportTrees(trees.toArray(new Tree[0]), true, getTreeNames(trees));
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

    //Main method
    public static void main(String[] args) throws IOException {

        String inputFileName = null;
        String outputFileName = null;
        String[] taxaToPrune = null;
        boolean basedOnNameContent = false;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.StringOption("taxaToPrune", "list","a list of taxon names to prune"),
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

        if (arguments.hasOption("taxaToPrune")) {
            taxaToPrune = Branch2dRateToGrid.parseVariableLengthStringArray(arguments.getStringOption("taxaToPrune"));
        }

        String nameContentString = arguments.getStringOption(NAMECONTENT);
        if (nameContentString != null && nameContentString.compareToIgnoreCase("true") == 0)
            basedOnNameContent = true;

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

        new TreePruner(inputFileName, outputFileName, taxaToPrune, basedOnNameContent);

        System.exit(0);
    }
}