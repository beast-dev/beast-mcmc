/*
 * TaxaMarkovJumpHistoryAnalyzer.java
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

package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.util.Version;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Philippe Lemey
 * @author Marc Suchard
 */
public class TaxaMarkovJumpHistoryAnalyzer extends BaseTreeTool {

    private final static Version version = new BeastVersion();

    private static final String HISTORY = "history";
    private static final String BURN_IN = "burnIn";

//    private static final boolean NEW_OUTPUT = true;

    private TaxaMarkovJumpHistoryAnalyzer(String inputFileName,
                                          String outputFileName,
                                          String[] taxaToProcess,
                                          String endState, // if no end state is provided, we will need to go to the root.
                                          double endTime,
                                          String endSister,
                                          String stateAnnotationName,
                                          double mrsd,
                                          int burnIn
    ) throws IOException {

        List<Tree> trees = new ArrayList<>();

        readTrees(trees, inputFileName, burnIn);

        List<Taxon> taxa = getTaxaToProcess(trees.get(0), taxaToProcess);

        this.mrsd = mrsd;
        this.stateAnnotationName = stateAnnotationName;

        this.ps = openOutputFile(outputFileName);
        processTrees(trees, taxa, endState, endTime, endSister, burnIn);
        closeOutputFile(ps);
    }

    private List<Taxon> getTaxaToProcess(Tree tree, String[] taxaToProcess) {
        List<Taxon> taxa = new ArrayList<>();
        if (taxaToProcess != null) {
            for (String name : taxaToProcess) {
                addTaxonByName(tree, taxa, name);
            }
        }
        return taxa;
    }

    private void processTrees(List<Tree> trees, List<Taxon> taxa, String endState, double endTime, String endSister, int burnIn) {
        if (burnIn < 0) {
            burnIn = 0;
        }
        for (int i = burnIn; i < trees.size(); ++i) {
            Tree tree = trees.get(i);
            processOneTree(tree, taxa, endState, endTime, endSister);
        }
    }

    private void processOneTree(Tree tree, List<Taxon> taxa, String endState, double endTime, String endSister) {

        String treeId = tree.getId();
        if (treeId.startsWith("STATE_")) {
            treeId = treeId.replaceFirst("STATE_", "");
        }

        for (Taxon taxon : taxa) {
            processOneTreeForOneTaxon(tree, taxon, endState, endTime, endSister, treeId);
        }
    }

    private void processOneTreeForOneTaxon(Tree tree, Taxon taxon, String endState, double endTime, String endSister, String treeId) {

        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            NodeRef tip = tree.getExternalNode(i);
            if (tree.getNodeTaxon(tip) == taxon) {
                processOneTip(tree, tip, endState, endTime, endSister, treeId, taxon.getId());
            }
        }
    }

    private class Row {
        String taxonId;
        String treeId;
        String location;
        double startTime;
        double endTime;

        private static final String DELIMITER = ",";

        private Row(String taxonId, String treeId,
                    String location, double startTime, double endTime) {
            this.taxonId = taxonId; this.treeId = treeId;
            this.location = location; this.startTime = startTime; this.endTime = endTime;
        }

        public String toString() {
            return taxonId + DELIMITER + treeId + DELIMITER + location + DELIMITER
                    + startTime + DELIMITER + endTime;
        }
    }

    private boolean pathDone(Tree tree, NodeRef node, String currentState, String endState, double endTime, String endSister, double startTime) {

        NodeRef endNode = tree.getRoot();
        if (endSister!=null){
            int taxonId = tree.getTaxonIndex(endSister);
            if (taxonId == -1) {
                throw new RuntimeException("Unable to find taxon '" + endSister + "'.");
            } else {
                for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
                    NodeRef sister = tree.getExternalNode(i);
                    if (tree.getNodeTaxon(sister) == tree.getTaxon(taxonId)) {
                        endNode = TreeUtils.getCommonAncestor(tree, sister, node);
                    }
                }
            }
        }

        if (endTime < Double.MAX_VALUE){
            if (mrsd < Double.MAX_VALUE) {
                return node == tree.getRoot() || currentState.equalsIgnoreCase(endState) || startTime<endTime || node.equals(endNode);
            } else {
                System.err.println("end time "+endTime+" (assumed in absolute time) is ignored because no most recent sampling date (mrsd) provided");
                return node == tree.getRoot() || currentState.equalsIgnoreCase(endState) || node.equals(endNode);
            }
        } else {
            return node == tree.getRoot() || currentState.equalsIgnoreCase(endState) || node.equals(endNode);
        }
    }

    private double adjust(double time) {
        if (mrsd < Double.MAX_VALUE) {
            time = mrsd - time;
        }
        return time;
    }

    private void processOneTip(Tree tree, NodeRef tip, String endState, double endTime, String endSister, String treeId, String taxonId) {

        String currentState = (String) tree.getNodeAttribute(tip, stateAnnotationName);
        if (currentState == null) {
            System.err.println("no " + stateAnnotationName + " annotation for tip");
            System.exit(-1);
        }

        double startTime = adjust(tree.getNodeHeight(tip));

        while (!pathDone(tree, tip, currentState, endState, endTime, endSister, startTime)) {

            Object[] jumps = readCJH(tip, tree);
            if (jumps != null) {

                for (int i = jumps.length - 1; i >= 0; i--) {
                    Object[] jump = (Object[]) jumps[i];

                    if (!currentState.equalsIgnoreCase((String) jump[2])) {
                        System.err.println(jump[1] + "\t" + jump[2]);
                        System.err.println("mismatch between states in Markov Jump History");
                        System.exit(-1);
                    }

                    double jumpTime = adjust((Double) jump[0]);

                    Row row = new Row(taxonId, treeId, currentState, startTime, jumpTime);
                    ps.println(row);

                    startTime = jumpTime;
                    currentState = ((String) jump[1]);
                }
            }
            tip = tree.getParent(tip);
        }

        double tipTime = adjust(tree.getNodeHeight(tip));
        Row row = new Row(taxonId, treeId, currentState, startTime, tipTime);

        ps.println(row);
    }

    private static Object[] readCJH(NodeRef node, Tree treeTime) {
        if (treeTime.getNodeAttribute(node, HISTORY) != null) {
            return (Object[]) treeTime.getNodeAttribute(node, HISTORY);
        } else {
            return null;
        }
    }

    protected PrintStream openOutputFile(String outputFileName) {

        PrintStream ps = super.openOutputFile(outputFileName);
        ps.println("taxonId,treeId,location,startTime,endTime");
        return ps;
    }

    private void closeOutputFile(PrintStream ps) {
        if (ps != null) {
            ps.close();
        }
    }

    private double mrsd;
    private String stateAnnotationName;

    private PrintStream ps;

    public static void printTitle() {
        progressStream.println();
        centreLine("TaxonMarkovJumpHistory " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("tool to get a Markov Jump History for a Taxon", 60);
        centreLine("by", 60);
        centreLine("Philippe Lemey and Marc Suchard", 60);
        progressStream.println();
        progressStream.println();
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("TaxonMarkovJumpHistory", "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: taxonMarkovJumpHistory -taxaToProcess taxon1,taxon2 input.trees output.trees");
        progressStream.println();
    }

    //Main method
    public static void main(String[] args) throws IOException {

        String[] taxaToProcess = null;
        String endState = null;
        double endTime = Double.MAX_VALUE;
        String endSister = null;
        String stateAnnotationName = "location";
        double mrsd = Double.MAX_VALUE;
        int burnIn = -1;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURN_IN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.StringOption("taxaToProcess", "list", "a list of taxon names to process MJHs"),
                        new Arguments.StringOption("endState", "end_state", "a state at which the MJH processing stops"),
                        new Arguments.RealOption("endTime", "a time at which the MJH processing stops"),
                        new Arguments.StringOption("endSister", "end_sister", "to stop at node that is the common ancestor of the relevant taxon and the specified sister taxon"),
                        new Arguments.StringOption("stateAnnotation", "state_annotation_name", "The annotation name for the discrete state string"),
                        new Arguments.RealOption("mrsd", "The most recent sampling time to convert heights to times [default=MAX_VALUE]"),
                        new Arguments.Option("help", "option to print this message"),
                });


        handleHelp(arguments, args, TaxaMarkovJumpHistoryAnalyzer::printUsage);

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

        if (arguments.hasOption("endTime")) {
            endTime = arguments.getRealOption("endTime");
        }

        if (arguments.hasOption("endSister")) {
            endSister = arguments.getStringOption("endSister");
        }

        if (arguments.hasOption(BURN_IN)) {
            burnIn = arguments.getIntegerOption(BURN_IN);
            System.err.println("Ignoring a burn-in of " + burnIn + " trees.");
        }

        String[] fileNames = getInputOutputFileNames(arguments, TaxaMarkovJumpHistoryAnalyzer::printUsage);

        new TaxaMarkovJumpHistoryAnalyzer(fileNames[0], fileNames[1], taxaToProcess,
                endState, endTime, endSister,
                stateAnnotationName,
                mrsd, burnIn);

        System.exit(0);
    }
}