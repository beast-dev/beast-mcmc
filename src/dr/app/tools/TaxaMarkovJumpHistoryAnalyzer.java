/*
 * TaxonMarkovJumpHistory.java
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
public class TaxaMarkovJumpHistoryAnalyzer {

    private final static Version version = new BeastVersion();

    // Messages to stderr, output to stdout
    private static PrintStream progressStream = System.err;

    public static final String HISTORY = "history";
    public static final String BURNIN = "burnin";

    private TaxaMarkovJumpHistoryAnalyzer(String inputFileName,
                                          String outputFileName,
                                          String[] taxaToProcess,
                                          String endState, // if no end state is provided, we will need to go to the root.
                                          String stateAnnotationName,
                                          double mrsd,
                                          int burnin
    ) throws IOException {

        List<Tree> trees = new ArrayList<>();

        readTrees(trees, inputFileName, burnin);

        List<Taxon> taxa = getTaxaToProcess(trees.get(0), taxaToProcess);

        historyStrings = new String[totalUsedTrees*taxa.size()];

        this.mrsd = mrsd;

        this.stateAnnotationName = stateAnnotationName;

        processTrees(trees, taxa, endState, burnin);

        writeOutputFile(outputFileName);
    }

    private void readTrees(List<Tree> trees, String inputFileName, int burnin) throws IOException {

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
                if (totalTrees > burnin){
                    totalUsedTrees++;
                }
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
        if (totalUsedTrees < 1) {
            System.err.println("No trees past burnin (="+burnin+")");
            return;
        }

        progressStream.println("Total trees read: " + totalTrees);
        progressStream.println("Total trees used: " + totalUsedTrees);
    }

    private List<Taxon> getTaxaToProcess(Tree tree, String[] taxaToProcess) {
        List<Taxon> taxa = new ArrayList<>();
        if (taxaToProcess != null) {
            for (String name : taxaToProcess) {
                int taxonId = tree.getTaxonIndex(name);
                if (taxonId == -1) {
                    throw new RuntimeException("Unable to find taxon '" + name + "'.");
                }
                taxa.add(tree.getTaxon(taxonId));
            }
        }
        return taxa;
    }

    private void processTrees(List<Tree> trees, List<Taxon> taxa, String endState, int burnin) {
        if (burnin < 0){
            burnin = 0;
        }
        for (int i = burnin; i <trees.size(); ++i) {
            Tree tree = trees.get(i);
            processOneTree(tree, taxa, endState, i-burnin);
        }
    }

    private void processOneTree(Tree tree, List<Taxon> taxa, String endState, int treeNumber) {
        int counter = 0;
        for (Taxon taxon : taxa) {
//            System.out.println(counter);
//            System.out.println(taxon.getId());
            processOneTreeForOneTaxon(tree, taxon, endState, treeNumber+(totalUsedTrees*counter));
            counter++;
        }
    }

    private void processOneTreeForOneTaxon(Tree tree, Taxon taxon, String endState, int treeNumber) {
        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            NodeRef tip = tree.getExternalNode(i);
            if (tree.getNodeTaxon(tip) == taxon) {
                sb.append(taxon.getId()+",");
                processOneTip(tree, tip, endState, treeNumber);
            }
        }
    }

    private void processOneTip(Tree tree, NodeRef tip, String endState, int treeNumber) {

        String currentState = (String)tree.getNodeAttribute(tip,stateAnnotationName);
        if (currentState==null){
            System.err.println("no "+stateAnnotationName+" annotation for tip");
            System.exit(-1);
        }
        sb.append(currentState+",");
        double nodeTime = tree.getNodeHeight(tip);
        if (mrsd < Double.MAX_VALUE){
            nodeTime = mrsd - nodeTime;
        }
        sb.append(nodeTime);
        while (tip != tree.getRoot()){
            Object[] jumps = readCJH(tip, tree);
            if (jumps != null){
                for (int i = jumps.length-1; i>=0; i--) {
                    Object[] jump = (Object[])jumps[i];
                    if (!currentState.equalsIgnoreCase((String)jump[2])){
                        System.err.println(jump[1]+"\t"+jump[2]);
                        System.err.println("mismatch between states in Markov Jump History");
                        System.exit(-1);
                    }
                    sb.append(",");
                    currentState = ((String)jump[1]);
                    sb.append(currentState+",");
                    double jumpTime = (Double)jump[0];
                    if (mrsd < Double.MAX_VALUE){
                        jumpTime = mrsd - jumpTime;
                    }
                    sb.append(jumpTime);
                }
            }
            tip = tree.getParent(tip);
            if ((endState!=null) && (currentState.equalsIgnoreCase(endState))){
                break;
            }
        }
//        System.out.println(sb.toString());
        historyStrings[treeNumber] = sb.toString();
        sb.setLength(0);
    }

    private static Object[] readCJH(NodeRef node, Tree treeTime){
        if(treeTime.getNodeAttribute(node, HISTORY)!=null){
            Object[] cjh = (Object[])treeTime.getNodeAttribute(node, HISTORY);
            return cjh;
        } else {
            return null;
        }
    }

    private void writeOutputFile(String outputFileName) {
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
        ps.println("#each trajectory starts with the tip name, its state, its height, and then the subsequent states adopted and their transition times up to either the root state or a pre-defined state");
        for (String historyString : historyStrings) {
            ps.println(historyString);
        }
        if (ps != null) {
            ps.close();
        }
    }


    int totalTrees = 0;
    int totalUsedTrees = 0;
    StringBuffer sb = new StringBuffer();
    String[] historyStrings;
    double mrsd;
    String stateAnnotationName;

    public static void printTitle() {
        progressStream.println();
        centreLine("TaxonMarkovJumpHistory " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("tool to get a Markov Jump History for a Taxon", 60);
        centreLine("by", 60);
        centreLine("Philippe Lemey and Marc Suchard", 60);
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

        arguments.printUsage("TaxonMarkovJumpHistory", "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: taxonMarkovJumpHistory -taxaTorocess taxon1,taxon2 input.trees output.trees");
        progressStream.println();
    }

    //Main method
    public static void main(String[] args) throws IOException {

        String inputFileName = null;
        String outputFileName = null;
        String[] taxaToProcess = null;
        String endState = null;
        String stateAnnotationName = "location";
        double mrsd = Double.MAX_VALUE;
        int burnin = -1;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURNIN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.StringOption("taxaToProcess", "list","a list of taxon names to process MJHs"),
                        new Arguments.StringOption("endState", "end_state","a state at which the MJH processing stops"),
                        new Arguments.StringOption("stateAnnotation", "state_annotation_name","The annotation name for the discrete state string"),
                        new Arguments.RealOption("mrsd", "The most recent sampling time to convert heights to times [default=MAX_VALUE]"),
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

        if (arguments.hasOption("taxaToProcess")) {
            taxaToProcess = Branch2dRateToGrid.parseVariableLengthStringArray(arguments.getStringOption("taxaToProcess"));
        }

        if (arguments.hasOption("endState")) {
            endState = arguments.getStringOption("endState");
        }

        if (arguments.hasOption("stateAnnotation")) {
            stateAnnotationName = arguments.getStringOption("stateAnnotation");
        }

        if (arguments.hasOption("mrsd")) {
            mrsd = arguments.getRealOption("mrsd");
        }

        if (arguments.hasOption(BURNIN)) {
            burnin = arguments.getIntegerOption(BURNIN);
            System.err.println("Ignoring a burnin of " + burnin + " trees.");
        }

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

        new TaxaMarkovJumpHistoryAnalyzer(inputFileName, outputFileName, taxaToProcess, endState, stateAnnotationName, mrsd, burnin);

        System.exit(0);
    }
}