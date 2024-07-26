/*
 * PersistenceSummarizer.java
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
import dr.util.Pair;
import dr.util.Version;

import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * @author Philippe Lemey
 * @author Marc Suchard
 */
public class TreeMRCASummarizer extends BaseTreeTool {

    private final static Version version = new BeastVersion();

    private static final String BURN_IN = "burnIn";

    private TreeMRCASummarizer(String inputFileName,
                               String outputFileName,
                               int pairCount
    ) throws IOException {

        SequentialTreeReader treeReader = new SequentialTreeReader(inputFileName, 0);

        this.ps = openOutputFile(outputFileName);
        processTrees(treeReader, pairCount);
        closeOutputFile(ps);
    }

    private void processTrees(SequentialTreeReader treeReader, int pairCount) throws IOException {
        int burnIn = 0;

        Tree tree = treeReader.getTree(0);

        Set<Pair<Taxon, Taxon>> pairs = new LinkedHashSet<>();
        while (pairs.size() < pairCount) {
            int r = (int)(Math.random() * tree.getTaxonCount());
            Taxon taxon1 = tree.getTaxon(r);
            Taxon taxon2 = taxon1;
            while (taxon2 == taxon1) {
                r = (int)(Math.random() * tree.getTaxonCount());
                taxon2 = tree.getTaxon(r);
            }
            pairs.add(new Pair<>(taxon1, taxon2));
        }

        ps.println("state");
        for (Pair<Taxon, Taxon> pair : pairs) {
            ps.println("\t" + pair.first + "+" + pair.second);
        }
        ps.println();


        int index = 1;
        while (treeReader.getTree(index) != null) {
            tree = treeReader.getTree(index);
            processOneTree(tree, pairs);
            index++;
        }
    }

    private void processOneTree(Tree tree, Collection<Pair<Taxon, Taxon>> pairs) {

        String treeId = tree.getId();
        if (treeId.startsWith("STATE_")) {
            treeId = treeId.replaceFirst("STATE_", "");
        }
        long state = Long.parseLong(treeId);

        ps.println(state);
        for (Pair<Taxon, Taxon> pair : pairs) {
            NodeRef tip1 = tree.getExternalNode(tree.getTaxonIndex(pair.first));
            NodeRef tip2 = tree.getExternalNode(tree.getTaxonIndex(pair.second));
            NodeRef mrca = TreeUtils.getCommonAncestor(tree, tip1, tip2);
//            double pathLength = TreeUtils.getPathLength(tree, tip1, tip2);
            ps.println("\t" + tree.getNodeHeight(mrca));
        }
        ps.println();

    }

    private void closeOutputFile(PrintStream ps) {
        if (ps != null) {
            ps.close();
        }
    }

    private PrintStream ps;

    public static void printTitle() {
        progressStream.println();
        centreLine("TreeMRCASummarizer " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("logs tMRCA for random pairs of taxa to a file", 60);
        centreLine("by", 60);
        centreLine("Andrew Rambaut", 60);
        progressStream.println();
        progressStream.println();
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("treemrcasummarizer", "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: treemrcasummarizer -pairs 100 input.trees output.log");
        progressStream.println();
    }

    private static final String PAIRS = "pairs";

    //Main method
    public static void main(String[] args) throws IOException, Arguments.ArgumentException {

        int pairs = 10;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(PAIRS, "the number of pairs of taxa [default = 10]"),
                        new Arguments.Option("help", "option to print this message"),
                });


        handleHelp(arguments, args, TreeMRCASummarizer::printUsage);

        if (arguments.hasOption(PAIRS)) {
            pairs = arguments.getIntegerOption(PAIRS);
        }
        System.err.println("Processing " + pairs + " pairs of taxa.");

        String[] fileNames = getInputOutputFileNames(arguments, TreeMRCASummarizer::printUsage);
        if (fileNames == null) {
            System.err.println("Missing input and output file names");
            printUsage(arguments);
            System.exit(1);
        }
        if (fileNames.length == 1) {
            System.err.println("Missing output file name");
            printUsage(arguments);
            System.exit(1);
        }

        new TreeMRCASummarizer(fileNames[0], fileNames[1], pairs);
        System.exit(0);
    }
}