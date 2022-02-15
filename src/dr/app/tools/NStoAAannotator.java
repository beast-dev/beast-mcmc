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
import dr.evolution.tree.*;
import dr.util.Version;
import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.GeneticCode;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import dr.evolution.datatype.*;

/**
 * @author Philippe Lemey
 * @author Marc Suchard
 */
public class NStoAAannotator extends BaseTreeTool {

    private final static Version version = new BeastVersion();

    private static final String NS = "all_N";
    private static final String BURN_IN = "burnIn";
    private static final String AA_ANNOT = "aa_";
    private static final String JOINTSTATE = "jointState";
    private static final String AA_JOINT = "aa_joint";
    private static final String[] falseTrue = {"false", "true"};

    private NStoAAannotator(String inputFileName,
                            String outputFileName,
                            String[] sitesToAnnotate,
                            int burnIn,
                            boolean jointState
    ) throws IOException {

        List<Tree> trees = new ArrayList<>();
        readTrees(trees, inputFileName, burnIn);
        outputTrees = new Tree[trees.size()-burnIn];

        processTrees(trees, sitesToAnnotate, burnIn);
        if (jointState){
            annotateJointState(sitesToAnnotate);
        }

        NexusExporter exporter = new NexusExporter(new PrintStream(outputFileName));
        exporter.exportTrees(outputTrees);

    }

    private void processTrees(List<Tree> trees, String[] sitesToAnnotate, int burnIn) {
        if (burnIn < 0) {
            burnIn = 0;
        }
        for (int i = burnIn; i < trees.size(); ++i) {
            Tree tree = trees.get(i);
            FlexibleTree treeToAnnotate = annotateOneTree(new FlexibleTree(tree, true), sitesToAnnotate);
            outputTrees[i - burnIn] = treeToAnnotate;
        }
    }

    private FlexibleTree annotateOneTree(FlexibleTree tree, String[] sitesToAnnotate) {
        for (int i = 0; i < sitesToAnnotate.length; i++) {
            String site = sitesToAnnotate[i];
            String rootState = findRootCodonState(tree, site);
            if (rootState == null){
                System.err.println("no nonsynonymous changes find for site "+site+"!!");
            }
            tree.setNodeAttribute(tree.getRoot(),AA_ANNOT+site, codonToAA(rootState));
            for (int j = 0; j < tree.getExternalNodeCount(); ++j) {
                FlexibleNode flexNode = (FlexibleNode)tree.getExternalNode(j);
                boolean annotated = hasAttribute(flexNode,AA_ANNOT+site);
                if(hasNSatSiteInPath(tree, flexNode, site)){
                    // as we have set the root state, the root should always be annotated
                    //first get the recentState from the most recent NS and propagate if no NS or NS at relevant site. modify recentState according to NS encountered
                    String recentState = findYoungestNSrecentState(tree,flexNode,site);
                    while(!annotated) {
                        Object[] NS = readNS(flexNode, tree);
                        if (NS != null) {
                            boolean NSatSite =  false;
                            for (int k = 0; k < NS.length; k++) {
                                Object[] singleNS = (Object[]) NS[k];
                                if (site.equalsIgnoreCase(singleNS[0].toString())) {
                                    if(!recentState.equalsIgnoreCase(singleNS[3].toString())) {
                                        System.err.println("current 'to State' does not match 'to State' in NS");
                                    }
                                    tree.setNodeAttribute(flexNode,AA_ANNOT+site, codonToAA(singleNS[3].toString()));
                                    //we re-set the 'recentState' to the 'from state' of the NS on this brance, which will be subsequently used if there are no NS in the ancestral brancehs
                                    recentState = singleNS[2].toString();
                                    NSatSite = true;
                                }
                            }
                            //if none of the NS are at the correct site
                            if(!NSatSite){
                                tree.setNodeAttribute(flexNode,AA_ANNOT+site, codonToAA(recentState));
                            }
                        //if there are no NS
                        } else {
                            tree.setNodeAttribute(flexNode,AA_ANNOT+site, codonToAA(recentState));
                        }
                        flexNode = (FlexibleNode)tree.getParent(flexNode);
                        annotated = hasAttribute(flexNode,AA_ANNOT+site);
                    }
                } else {
                     while(!annotated) {
                        tree.setNodeAttribute(flexNode,AA_ANNOT+site, codonToAA(rootState));
                        flexNode = (FlexibleNode)tree.getParent(flexNode);
                        annotated = hasAttribute(flexNode,AA_ANNOT+site);
                    }
                }
            }
        }
        return tree;
    }

    private void annotateJointState(String[] sitesToAnnotate) {
        for (int i = 0; i < outputTrees.length; ++i) {
            Tree tree = outputTrees[i];
            FlexibleTree treeToAnnotate = annotateJointStateOneTree(new FlexibleTree(tree, true), sitesToAnnotate);
            outputTrees[i] = treeToAnnotate;
        }
    }
    private FlexibleTree annotateJointStateOneTree(FlexibleTree tree, String[] sitesToAnnotate){
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            FlexibleNode flexNode = (FlexibleNode)tree.getExternalNode(i);
            String jointState = "";
            for (int j = 0; j < sitesToAnnotate.length; ++j) {
                char siteState = (Character) tree.getNodeAttribute(flexNode, AA_ANNOT+sitesToAnnotate[j]);
                jointState = jointState + siteState;
            }
            tree.setNodeAttribute(flexNode,AA_JOINT, jointState);
        }
        return tree;
    }

    private boolean hasAttribute(FlexibleNode flexNode, String annotation){
        if (flexNode.getAttribute(annotation) != null) {
            return true;
        } else {
            return false;
        }
    }

    private static Object[] readNS(NodeRef node, Tree treeTime) {
        if(treeTime.isRoot(node)){
            System.err.println("root node!");
        }
        if (treeTime.getNodeAttribute(node, NS) != null) {
            return (Object[]) treeTime.getNodeAttribute(node, NS);
        } else {
            return null;
        }
    }

    private char codonToAA(String codon) {
        int nuc1Int = Nucleotides.INSTANCE.getState(codon.charAt(0));
        int nuc2Int = Nucleotides.INSTANCE.getState(codon.charAt(1));
        int nuc3Int = Nucleotides.INSTANCE.getState(codon.charAt(2));
        int codonInt =  (nuc1Int * 16) + (nuc2Int * 4) + nuc3Int;
        int aaInt = geneticCode.getAminoAcidState(codons.getCanonicalState(codonInt));
        return AminoAcids.INSTANCE.getChar(aaInt);
    }

    private String findRootCodonState(Tree tree, String site){
        String codonState = null;
        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            NodeRef extNode = tree.getExternalNode(i);
            if(hasNSatSiteInPath(tree,extNode,site)){
                codonState = findOldestNSfromState(tree,extNode,site);
            }
        }
        return codonState;
    }

    private boolean hasNSatSiteInPath(Tree tree, NodeRef extNode, String site){
        boolean hasNS = false;
        while(!hasNS && !tree.isRoot(extNode)) {
            Object[] NS = readNS(extNode, tree);
            if (NS != null) {
                for (int i = 0; i < NS.length; i ++) {
                    Object[] singleNS = (Object[]) NS[i];
                    if(site.equalsIgnoreCase(singleNS[0].toString())){
                        hasNS = true;
                    }
                }
            }
            extNode = tree.getParent(extNode);
        }
        return hasNS;
    }

    private String findOldestNSfromState(Tree tree, NodeRef extNode, String site){
        String fromState = null;
        double NStime = 0;
        while(!tree.isRoot(extNode)) {
            Object[] NS = readNS(extNode, tree);
            if (NS != null) {
                for (int i = 0; i < NS.length; i++) {
                    Object[] singleNS = (Object[]) NS[i];
                    if (site.equalsIgnoreCase(singleNS[0].toString())) {
                        double singleNStime = Double.parseDouble(singleNS[1].toString());
                        if(singleNStime > NStime){
                            NStime = singleNStime;
                            fromState = singleNS[2].toString();
                        }
                    }
                }
            }
            extNode = tree.getParent(extNode);
        }
        return fromState;
    }

    private String findYoungestNSrecentState(Tree tree, NodeRef extNode, String site){
        String recentState = null;
        boolean NSfound = false;
        while(!NSfound && !tree.isRoot(extNode)) {
            Object[] NS = readNS(extNode, tree);
            if (NS != null) {
                for (int i = 0; i < NS.length; i++) {
                    Object[] singleNS = (Object[]) NS[i];
                    if (site.equalsIgnoreCase(singleNS[0].toString())) {
                        recentState = singleNS[3].toString();
                        NSfound = true;
                    }
                }
            }
            extNode = tree.getParent(extNode);
        }
        return recentState;
    }

    private Tree[] outputTrees;
    GeneticCode geneticCode = GeneticCode.UNIVERSAL;
    Codons codons = Codons.UNIVERSAL;

    public static void printTitle() {
        progressStream.println();
        centreLine("NStoAAannotator " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("tool to annotate AA node states based on NS substitutions", 60);
        centreLine("by", 60);
        centreLine("Philippe Lemey and Marc Suchard", 60);
        progressStream.println();
        progressStream.println();
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("NStoAAannotator", "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: NStoAAannotator -sitesToAnnotate site1,site2 input.trees output.trees");
        progressStream.println();
    }

    //Main method
    public static void main(String[] args) throws IOException {

        String[] sitesToAnnotate = null;
        int burnIn = -1;
        boolean annotateJointState = false;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURN_IN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.StringOption("sitesToAnnotate", "list", "a list of sites to annotate"),
                        new Arguments.StringOption(JOINTSTATE, falseTrue, false,
                                "annotate the joint State for the sites [default = false])"),
                        new Arguments.Option("help", "option to print this message"),
                });

        handleHelp(arguments, args, NStoAAannotator::printUsage);

        if (arguments.hasOption("sitesToAnnotate")) {
            sitesToAnnotate = Branch2dRateToGrid.parseVariableLengthStringArray(arguments.getStringOption("sitesToAnnotate"));
        }

        if (arguments.hasOption(BURN_IN)) {
            burnIn = arguments.getIntegerOption(BURN_IN);
            System.err.println("Ignoring a burn-in of " + burnIn + " trees.");
        }

        String jointStateString = arguments.getStringOption(JOINTSTATE);
        if (jointStateString != null && jointStateString.compareToIgnoreCase("true") == 0)
            annotateJointState = true;

        String[] fileNames = getInputOutputFileNames(arguments, NStoAAannotator::printUsage);

        new NStoAAannotator(fileNames[0], fileNames[1], sitesToAnnotate, burnIn, annotateJointState);

        System.exit(0);
    }
}