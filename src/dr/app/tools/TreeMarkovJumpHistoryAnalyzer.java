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
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.util.Version;

import java.io.*;
import java.util.*;

/**
 * @author Philippe Lemey
 * @author Marc Suchard
 */
public class TreeMarkovJumpHistoryAnalyzer extends BaseTreeTool {

    private final static Version version = new BeastVersion();

    private static final String HISTORY = "history";
    private static final String BURN_IN = "burnIn";
    private static final String[] falseTrue = {"false", "true"};

//    private static final boolean NEW_OUTPUT = true;

    private TreeMarkovJumpHistoryAnalyzer(String inputFileName,
                                          String outputFileName,
                                          int burnIn,
                                          String[] taxaToIgnore,
                                          String[] statesToIgnore,
                                          boolean basedOnNameContent,
                                          String[] mrcaTaxa,
                                          double mrsd,
                                          String nodeStateAnnotation,
                                          String checkMJforState
     ) throws IOException {

        List<Tree> trees = new ArrayList<>();

        readTrees(trees, inputFileName, burnIn);

        Set ignoredTaxa = getTaxaSet(trees.get(0), taxaToIgnore, basedOnNameContent);
        Set commonAncestorTaxa = getTaxaSet(trees.get(0), mrcaTaxa, false);
        Set ignoredStates = stringArrayToStateSet(statesToIgnore);

        this.mrsd = mrsd;

        this.ps = openOutputFile(outputFileName);
        processTrees(trees, burnIn, ignoredTaxa, ignoredStates, commonAncestorTaxa, nodeStateAnnotation, checkMJforState);
        closeOutputFile(ps);
    }

    private Set getTaxaSet(Tree tree, String[] taxaToProcess, boolean basedOnContent) {
        Set taxa = new HashSet();
        getTaxaToFromName(tree, taxa, taxaToProcess, basedOnContent);
        return taxa;
    }

    private void processTrees(List<Tree> trees, int burnIn, Set ignoredTaxa, Set ignoredStates, Set commonAncestorTaxa, String nodeStateAnnotation, String checkMJforState) {
        if (burnIn < 0) {
            burnIn = 0;
        }
        for (int i = burnIn; i < trees.size(); ++i) {
            Tree tree = trees.get(i);
//            System.out.println(i);
            processOneTree(tree, ignoredTaxa, ignoredStates, commonAncestorTaxa, nodeStateAnnotation, checkMJforState);
        }
    }

    private void processOneTree(Tree tree, Set ignoredTaxa, Set ignoredStates, Set commonAncestorTaxa, String nodeStateAnnotation, String checkMJforState) {

        String treeId = tree.getId();
        if (treeId.startsWith("STATE_")) {
            treeId = treeId.replaceFirst("STATE_", "");
        }

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            NodeRef node = tree.getNode(i);
            if (inClade(tree, node, commonAncestorTaxa, ignoredTaxa)){
                if (nodeToConsider(tree, node, ignoredTaxa)){
                    if (nodeStateAnnotation==null) {
                        Object[] jumps = readCJH(node, tree);

                        if(checkMJforState!=null){
                            if(!tree.isRoot(node)){
                                String nodeState = (String) tree.getNodeAttribute(node, checkMJforState);
                                String parentNodeString = (String) tree.getNodeAttribute(tree.getParent(node), checkMJforState);
                                if (!nodeState.equalsIgnoreCase(parentNodeString)){
                                    if (jumps == null) {
//                                        System.out.println(nodeState+"\t"+parentNodeString);
                                        System.err.println("parent ("+parentNodeString+") and node ("+nodeState+") state are different, but no MJ????");
//                                    System.exit(-1);
                                    } else {
                                        //System.out.println("all fine");
                                    }
                                }
                            }
                        }

                        if (jumps != null) {
                            for (int j = jumps.length - 1; j >= 0; j--) {
                                Object[] jump = (Object[]) jumps[j];
                                if (!ignoredStates.contains(jump[1]) && !ignoredStates.contains(jump[2])){
                                    Row row = new Row(treeId, (String) jump[1], (String) jump[2], adjust((Double) jump[0]));
                                    ps.println(row);
                                }
                            }
                        }
                    } else {
                        Object[] jump = getMJApproximation(node, tree, nodeStateAnnotation);
                        if (jump != null) {
                            if (!ignoredStates.contains(jump[1]) && !ignoredStates.contains(jump[2])) {
                                Row row = new Row(treeId, (String) jump[1], (String) jump[2], adjust((Double) jump[0]));
                                ps.println(row);
                            }
                        }
                    }
                }
            }
        }
    }

    private boolean nodeToConsider(Tree tree, NodeRef node, Set taxa){
        boolean consider = true;
        if (taxa.isEmpty()){
            return consider;
        } else {
            Set descendants = TreeUtils.getDescendantLeaves(tree,node);
            Set taxaDescendants = stringToTaxaSet(tree, descendants);
            taxaDescendants.removeAll(taxa);
            if (taxaDescendants.isEmpty()){
                consider = false;
            }
            return consider;
        }
    }

    private boolean inClade(Tree tree, NodeRef node, Set commonAncestorTaxa, Set ignoredTaxa){
        boolean inClade = false;
        if (commonAncestorTaxa.isEmpty()){
            inClade = true;
        } else {
            Set descendants = TreeUtils.getDescendantLeaves(tree,node);
            Set taxaDescendants = stringToTaxaSet(tree, descendants);
            taxaDescendants.removeAll(ignoredTaxa);
            taxaDescendants.removeAll(commonAncestorTaxa);
            if (taxaDescendants.isEmpty()){
                inClade = true;
//                System.out.println("in clade!");
            }
        }
        return inClade;
    }

    private Set stringToTaxaSet(Tree tree, Set taxa){
        Set taxaStrings = new HashSet();
        Iterator<String> itr = taxa.iterator();
        while(itr.hasNext()){
            String descendant = itr.next();
            Taxon taxon = tree.getTaxon(tree.getTaxonIndex(descendant));
            taxaStrings.add(taxon);
        }
        return taxaStrings;
    }

    private Set stringArrayToStateSet(String[] states){
        Set stateStringSet = new HashSet();
        if (states==null){
            return stateStringSet;
        } else {
            for (int i=0; i< states.length; i++){
                stateStringSet.add(states[i]);
            }
        }
        return stateStringSet;
    }

    private class Row {
        String treeId;
        String startLocation;
        String endLocation;
        double time;

        private static final String DELIMITER = ",";

        private Row(String treeId,
                    String startLocation, String endLocation,
                            double time) {
            this.treeId = treeId;
            this.startLocation = startLocation; this.endLocation = endLocation;
            this.time = time;
        }

        public String toString() {
            return treeId + DELIMITER
                    + startLocation + DELIMITER + endLocation + DELIMITER
                    + time;
        }
    }


    private double adjust(double time) {
        if (mrsd < Double.MAX_VALUE) {
            time = mrsd - time;
        }
        return time;
    }


    private static Object[] readCJH(NodeRef node, Tree treeTime) {
        if (treeTime.getNodeAttribute(node, HISTORY) != null) {
            return (Object[]) treeTime.getNodeAttribute(node, HISTORY);
        } else {
            return null;
        }
    }

    private static Object[] getMJApproximation(NodeRef node, Tree treeTime, String nodeStateAnnotation) {
        String nodeState = (String) treeTime.getNodeAttribute(node, nodeStateAnnotation);
        String parentNodeString = (String) treeTime.getNodeAttribute(treeTime.getParent(node), nodeStateAnnotation);
        if (nodeState.equalsIgnoreCase(parentNodeString)) {
            return null;
        } else {
            double nodeHeight = treeTime.getNodeHeight(node);
            double parentNodeHeight = treeTime.getNodeHeight(treeTime.getParent(node));
            double x = (Math.random() * (parentNodeHeight - nodeHeight)) + nodeHeight;
            return new Object[] {x,parentNodeString,nodeState};
        }
    }

    protected PrintStream openOutputFile(String outputFileName) {

        PrintStream ps = super.openOutputFile(outputFileName);
        ps.println("treeId,startLocation,endLocation,time");
        return ps;
    }

    private void closeOutputFile(PrintStream ps) {
        if (ps != null) {
            ps.close();
        }
    }

    private double mrsd;

    private PrintStream ps;

    public static void printTitle() {
        progressStream.println();
        centreLine("TreeMarkovJumpHistory " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("tool to get a Markov Jump History for a tree", 60);
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

        String[] taxaToIgnore = null;
        String[] statesToIgnore = null;
        double mrsd = Double.MAX_VALUE;
        int burnIn = -1;
        String[] mrcaTaxa = null;
        String nodeStateAnnotation = null;
        String checkMJforState = null;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURN_IN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.StringOption("taxaToIgnore", "list", "a list of taxon names that defines parts of trees to ignore in MJH processing"),
                        new Arguments.StringOption("statesToIgnore", "list", "a list of location states to ignore in MJH processing"),
                        new Arguments.RealOption("mrsd", "The most recent sampling time to convert heights to times [default=MAX_VALUE]"),
                        new Arguments.StringOption(NAME_CONTENT, falseTrue, false,
                                "taxa names contain string [default = false])"),
                        new Arguments.StringOption("mrcaTaxa", "list", "a list of taxon names that defines a clade to focus on for MJH processing"),
                        new Arguments.StringOption("nodeStateAnnotation", "String", "use node state annotations as poor proxy to MJs based on a annotation string for the discrete trait"),
                        new Arguments.StringOption("checkMJforState", "String", "check if there is a MJ on branches with different parent and node state for this annotation"),
                        new Arguments.Option("help", "option to print this message"),
                });


        handleHelp(arguments, args, TaxaMarkovJumpHistoryAnalyzer::printUsage);

        if (arguments.hasOption("taxaToIgnore")) {
            taxaToIgnore = Branch2dRateToGrid.parseVariableLengthStringArray(arguments.getStringOption("taxaToIgnore"));
        }

        if (arguments.hasOption("statesToIgnore")) {
            statesToIgnore = Branch2dRateToGrid.parseVariableLengthStringArray(arguments.getStringOption("statesToIgnore"));
        }

        if (arguments.hasOption("mrsd")) {
            mrsd = arguments.getRealOption("mrsd");
        }

        if (arguments.hasOption(BURN_IN)) {
            burnIn = arguments.getIntegerOption(BURN_IN);
            System.err.println("Ignoring a burn-in of " + burnIn + " trees.");
        }

        if (arguments.hasOption("nodeStateAnnotation")) {
            nodeStateAnnotation = arguments.getStringOption("nodeStateAnnotation");
        }

        if (arguments.hasOption("checkMJforState")) {
            checkMJforState = arguments.getStringOption("checkMJforState");
        }
        if (arguments.hasOption("mrcaTaxa")) {
            mrcaTaxa = Branch2dRateToGrid.parseVariableLengthStringArray(arguments.getStringOption("mrcaTaxa"));
        }


        boolean basedOnNameContent = parseBasedOnNameContent(arguments);

        String[] fileNames = getInputOutputFileNames(arguments, TaxaMarkovJumpHistoryAnalyzer::printUsage);

        new TreeMarkovJumpHistoryAnalyzer(fileNames[0], fileNames[1], burnIn,
                taxaToIgnore,
                statesToIgnore,
                basedOnNameContent,
                mrcaTaxa,
                mrsd,
                nodeStateAnnotation,
                checkMJforState
        );

        System.exit(0);
    }
}