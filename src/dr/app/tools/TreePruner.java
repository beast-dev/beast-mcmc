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
public class TreePruner extends BaseTreeTool {

    private final static Version version = new BeastVersion();

    private static final String[] falseTrue = {"false", "true"};
    private static final String NAME_CONTENT = "nameContent";

    private TreePruner(String inputFileName,
                       String outputFileName,
                       String[] taxaToPrune,
                       boolean basedOnNameContent
                       ) throws IOException {

        List<Tree> trees = new ArrayList<>();

        readTrees(trees, inputFileName, 0);

        List<Taxon> taxa = getTaxaToPrune(trees.get(0), taxaToPrune, basedOnNameContent);

        processTrees(trees, taxa);

        writeOutputFile(trees, outputFileName);
    }

    private List<Taxon> getTaxaToPrune(Tree tree, String[] names, boolean basedOnContent) {

        List<Taxon> taxa = new ArrayList<>();
        getTaxaToFromName(tree, taxa, names, basedOnContent);
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

        if (parent == tree.getRoot()) {

            //sibling must become new root
            NodeRef sibling = getSibling(tree, parent, tip);

            FlexibleNode siblingNode = (FlexibleNode) sibling;
            siblingNode.setParent(null);

            tree.beginTreeEdit();
            tree.setRoot(sibling);
            tree.endTreeEdit();

            //throw new RuntimeException("Still need to handle this situation");

        } else {

            NodeRef grandParent = tree.getParent(parent);
            NodeRef sibling = getSibling(tree, parent, tip);

            FlexibleNode grandParentNode = (FlexibleNode)grandParent;
            FlexibleNode parentNode = (FlexibleNode) parent;
            FlexibleNode siblingNode = (FlexibleNode) sibling;
            //FlexibleNode tipNode = (FlexibleNode) tip;

            // Remove from topology
            siblingNode.setParent(grandParentNode);
            grandParentNode.removeChild(parentNode);
            grandParentNode.addChild(siblingNode);

            // Adjust branch lengths
            siblingNode.setLength(parentNode.getLength() + siblingNode.getLength());

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

    private void writeOutputFile(List<Tree> trees, String outputFileName) {


        PrintStream ps = openOutputFile(outputFileName);

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
        progressStream.println();
        progressStream.println();
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("treepruner", "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: treepruner -taxaToPrune taxon1,taxon2 input.trees output.trees");
        progressStream.println();
    }

    //Main method
    public static void main(String[] args) throws IOException {

        String[] taxaToPrune = null;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.StringOption("taxaToPrune", "list","a list of taxon names to prune"),
                        new Arguments.StringOption(NAME_CONTENT, falseTrue, false,
                                "add true noise [default = true])"),
                        new Arguments.Option("help", "option to print this message"),
                });

        handleHelp(arguments, args, TreePruner::printUsage);

        if (arguments.hasOption("taxaToPrune")) {
            taxaToPrune = Branch2dRateToGrid.parseVariableLengthStringArray(arguments.getStringOption("taxaToPrune"));
        }

        boolean basedOnNameContent = parseBasedOnNameContent(arguments);

        String[] fileNames = getInputOutputFileNames(arguments, TreePruner::printUsage);

        new TreePruner(fileNames[0], fileNames[1], taxaToPrune, basedOnNameContent);

        System.exit(0);
    }
}