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

    private final static boolean USE_R = false;

    private static boolean forceIntegerToDiscrete = false;

    private double maxState = 1;


    // Messages to stderr, output to stdout
    private static PrintStream progressStream = System.err;


    /**
     * Burnin can be specified as the number of trees or the number of states
     * (one or other should be zero).
     *
     * @param burninTrees
     * @param burninStates
     * @param inputFileName
     * @param outputFileName
     * @throws IOException
     */
    public TreePruner(final int burninTrees,
                      final long burninStates,
                      String inputFileName,
                      String outputFileName,
                      String[] taxaToPrune

    ) throws IOException {

        int burnin = -1;

        totalTrees = 10000;
        totalTreesUsed = 0;

        progressStream.println("Reading trees (bar assumes 10,000 trees)...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        long stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        FileReader fileReader = new FileReader(inputFileName);
        TreeImporter importer = new NexusImporter(fileReader, false);

        try {
            totalTrees = 0;
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                addTree(tree);

                long state = Long.MAX_VALUE;

                if (burninStates > 0) {
                    // if burnin has been specified in states, try to parse it out...
                    String name = tree.getId().trim();

                    if (name != null && name.length() > 0 && name.startsWith("STATE_")) {
                        state = Long.parseLong(name.split("_")[1]);
                        maxState = state;
                    }
                }

                if (totalTrees >= burninTrees && state >= burninStates) {
                    // if either of the two burnin thresholds have been reached...

                    if (burnin < 0) {
                        // if this is the first time this point has been reached,
                        // record the number of trees this represents for future use...
                        burnin = totalTrees;
                    }
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

        progressStream.println("Total trees read: " + totalTrees);
        if (burninTrees > 0) {
            progressStream.println("Ignoring first " + burninTrees + " trees" +
                    (burninStates > 0 ? " (" + burninStates + " states)." : "."));
        } else if (burninStates > 0) {
            progressStream.println("Ignoring first " + burninStates + " states (" + burnin + " trees).");
        }
        progressStream.println();

        processTrees(taxaToPrune);

        writeOutputFile(outputFileName);

    }

    private void processTrees(String[] taxaToPrune) {
        // TODO
    }

    private void addTree(Tree tree) {
        if (trees == null) {
            trees = new ArrayList<>();
        }
        trees.add(tree);
    }

    int totalTrees = 0;
    int totalTreesUsed = 0;

    private void writeOutputFile(String outputFileName) {

        PrintStream ps = null;
        try {
            ps = new PrintStream(new File(outputFileName));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        NexusExporter exporter = new NexusExporter(ps);

        if (trees.size() > 0) {
            exporter.exportTrees(trees.toArray(new Tree[trees.size()]), true, getTreeNames(trees));
        }

        ps.close();
    }

    private List<Tree> trees;

    private String[] getTreeNames(List<Tree> trees) {
        List<String> names = new ArrayList<>();
        for (Tree tree : trees) {
            names.add(tree.getId());
        }
        return names.toArray(new String[names.size()]);
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
        progressStream.println("  Example: treepruner 1000000 test.trees taxa.txt output.trees");
        progressStream.println();
    }

    private static String[] parseVariableLengthStringArray(String inString) throws Arguments.ArgumentException{

        List<String> returnList = new ArrayList<String>();
        StringTokenizer st = new StringTokenizer(inString, ",");
        while (st.hasMoreTokens()) {
            try {
            returnList.add(st.nextToken());
            } catch (NumberFormatException e) {
                throw new Arguments.ArgumentException();
            }
        }

        if (returnList.size() > 0) {
            String[] stringArray = new String[returnList.size()];
            stringArray = returnList.toArray(stringArray);
            return stringArray;
        }
        return null;
    }


    //Main method
    public static void main(String[] args) throws IOException {

        String inputFileName = null;
        String outputFileName = null;
        String[] taxaToPrune = null;
        long burninStates = -1;
        int burninTrees = -1;

        if (args.length == 0) {

            final String versionString = version.getVersionString();
            String nameString = "TreePruner " + versionString;
            String aboutString = "<html><center><p>" + versionString + ", " + version.getDateString() + "</p>" +
                    "<p>by<br>" +
                    "Philippe Lemey, Andrew Rambaut and Marc Suchard</p>" +
//                    "<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
//                    "<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
//                    "<p>Department of Computer Science, University of Auckland<br>" +
//                    "<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +
                    "<p>Part of the BEAST package:<br>" +
                    "<a href=\"http://beast.community\">http://beast.community</a></p>" +
                    "</center></html>";

            progressStream = System.out;

            printTitle();

            try {
                new TreePruner(
                        burninTrees,
                        burninStates,
                        inputFileName,
                        outputFileName,
                        taxaToPrune
                        );

            } catch (Exception ex) {
                System.err.println("Exception: " + ex.getMessage());
            }

        }

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.LongOption("burnin", "the number of states to be considered as 'burn-in'"),
                        new Arguments.IntegerOption("burninTrees", "the number of trees to be considered as 'burn-in'"),
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

        if (arguments.hasOption("burnin")) {
            burninStates = arguments.getLongOption("burnin");
        }
        if (arguments.hasOption("burninTrees")) {
            burninTrees = arguments.getIntegerOption("burninTrees");
        }

        if (arguments.hasOption("taxaToPrune")) {
            try {
                taxaToPrune = parseVariableLengthStringArray(arguments.getStringOption("taxaToPrune"));
            } catch (Arguments.ArgumentException e) {
                System.err.println("Error reading " + arguments.getStringOption("taxaToPrune"));
            }
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

        new TreePruner(burninTrees, burninStates, inputFileName, outputFileName, taxaToPrune);

        System.exit(0);
    }
}


