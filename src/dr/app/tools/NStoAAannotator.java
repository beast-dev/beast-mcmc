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

    private NStoAAannotator(String inputFileName,
                            String outputFileName,
                            String[] sitesToAnnotate,
                            int burnIn
    ) throws IOException {

        List<Tree> trees = new ArrayList<>();
        readTrees(trees, inputFileName, burnIn);
        outputTrees = new Tree[trees.size()-burnIn];

        processTrees(trees, sitesToAnnotate, burnIn);

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

        //TODO: annotate all nodes and not justt the ones for a branch with a NS substitution
        for (int i = 0; i < tree.getNodeCount(); i++) {
            FlexibleNode flexNode = (FlexibleNode)tree.getNode(i);
            Object[] NS = readNS(flexNode, tree);
            if (NS != null) {
                for (int j = 0; j < NS.length; j ++) {
                    Object[] singleNS = (Object[]) NS[j];
                    for (int k = 0; k < sitesToAnnotate.length; k ++) {
                        if (sitesToAnnotate[k].equalsIgnoreCase(singleNS[0].toString())){
//                           tree.setNodeAttribute(flexNode,"codon"+sitesToAnnotate[k], singleNS[3]);
//                            tree.setNodeAttribute(flexNode.getParent(),
//                                    "codon"+sitesToAnnotate[k], singleNS[2]);
                            tree.setNodeAttribute(flexNode,"aa_"+sitesToAnnotate[k], codonToAA((String)singleNS[3]));
                            tree.setNodeAttribute(flexNode.getParent(),
                                    "aa_"+sitesToAnnotate[k], codonToAA((String)singleNS[3]));
                        }
                    }
                }
            }
        }
        return tree;
    }

    private static Object[] readNS(NodeRef node, Tree treeTime) {
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

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURN_IN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.StringOption("sitesToAnnotate", "list", "a list of sites to annotate"),
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

        String[] fileNames = getInputOutputFileNames(arguments, NStoAAannotator::printUsage);

        new NStoAAannotator(fileNames[0], fileNames[1], sitesToAnnotate, burnIn);

        System.exit(0);
    }
}