/*
 * TreeLogAnalyser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.app.util.Utils;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeTraceAnalysis;
import dr.util.Version;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class TreeLogAnalyser {

    private final static Version version = new BeastVersion();

    static boolean combine = true;

    public TreeLogAnalyser(int burnin, String inputFileName, String outputFileName, String trueTreeFileName,
                           String exportFileName, double minSupport, double credibleSetProbability, int maxExport, boolean verbose) throws IOException {

        List<File> files = new ArrayList<File>();
        File inputFile = new File(inputFileName);

        if (inputFile.isDirectory()) {
            System.out.println("Analysing all tree files below directory: " + inputFileName);

            collectFiles(inputFile, files);
        } else if (inputFile.isFile()) {
            System.out.println("Analysing tree file: " + inputFileName);

            files.add(inputFile);
        } else {
            System.err.println(inputFileName + " does not exist!");
            System.exit(0);
        }

        if( files.size() == 0 ) {
           System.err.println("No valid files");
           System.exit(0);
        }

        if (outputFileName != null) {
            FileOutputStream outputStream = new FileOutputStream(outputFileName);
            System.setOut(new PrintStream(outputStream));
        }

        Tree trueTree = null;
        if (trueTreeFileName != null) {
            NexusImporter importer = new NexusImporter(new FileReader(trueTreeFileName));
            try {
                trueTree = importer.importNextTree();
            } catch (Importer.ImportException e) {
                throw new IOException(e.getMessage());
            }
        }

        analyze(files, burnin, trueTree, verbose, exportFileName, minSupport, credibleSetProbability, maxExport, new boolean[]{true});
    }

    private static void collectFiles(File file, List<File> files) {

        if (file.isFile()) {
            if (file.getName().endsWith(".tre") || file.getName().endsWith(".trees") || file.getName().endsWith(".t")) {
                files.add(file);
            }
        } else {
            File[] listFiles = file.listFiles();
            for (File listFile : listFiles) {
                collectFiles(listFile, files);
            }
        }
    }

    private static void analyze(List<File> files, int burnin, Tree tree, boolean verbose, String exportFileName,
                                double minSupport, double credibleSetProbability, int maxExport, boolean[] drawHeader) {

        if (combine) {
            try {
                Reader[] readers = new Reader[files.size()];
                for (int i = 0; i < readers.length; i++) {
                    readers[i] = new FileReader(files.get(i));
                }
                TreeTraceAnalysis analysis = TreeTraceAnalysis.analyzeLogFile(readers, burnin, verbose);
                if (exportFileName != null) {
                    PrintStream exportStream = new PrintStream(exportFileName);
                    //System.err.println("Exporting trees ...");
                    analysis.export(exportStream, minSupport, maxExport, verbose);
                } else {
                    if (verbose) {
                        analysis.report(0.05, credibleSetProbability, (int)(minSupport+.5));
                    } else {
                        final String name = files.size() > 1 ? "combined" : files.get(0).toString();
                        analysis.shortReport(name, tree, drawHeader[0]);
                        drawHeader[0] = false;
                    }
                }

            } catch (IOException ioe) {
                //
            }
        } else {
            for (File file : files) {
                try {
                    final Reader[] readers = {new FileReader(file)};
                    TreeTraceAnalysis analysis = TreeTraceAnalysis.analyzeLogFile(readers, burnin, verbose);
                    if (verbose) {
                        analysis.report((int)(minSupport+.5));
                    } else {
                        analysis.shortReport(file.toString(), tree, drawHeader[0]);
                        drawHeader[0] = false;
                    }
                } catch (IOException ioe) {
                    //
                }
            }
        }
    }

    public static void printTitle() {
        System.out.println();
        centreLine("TreeLogAnalyser " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("MCMC Output analysis", 60);
        centreLine("by", 60);
        centreLine("Alexei Drummond and Andrew Rambaut", 60);
        System.out.println();
        centreLine("Department of Computer Science", 60);
        centreLine("University of Auckland", 60);
        centreLine("alexei@cs.auckland.ac.nz", 60);
        System.out.println();
        centreLine("Institute of Evolutionary Biology", 60);
        centreLine("University of Edinburgh", 60);
        centreLine("a.rambaut@ed.ac.uk", 60);
        System.out.println();
        System.out.println();
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            System.out.print(" ");
        }
        System.out.println(line);
    }


    public static void printUsage(Arguments arguments) {

        arguments.printUsage("treeloganalyser", "<input-file-name> [<true-tree-file-name> [<output-file-name>]]");
        System.out.println();
        System.out.println("  Example: treeloganalyser test.trees trueTree.tree out.txt");
        System.out.println();
    }

    //Main method
    public static void main(String[] args) throws java.io.IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in' [default = none]"),
                        new Arguments.StringOption("export", "file-name", "name of file to export"),
                        new Arguments.RealOption("limit", "don't export trees with support lower than limit [default = 0.0]"),
                        new Arguments.RealOption("probability", "credible set probability limit [default = 0.95]"),
                        new Arguments.IntegerOption("max", "export no more than max trees [default = all]"),
                        new Arguments.Option("short", "use this option to produce a short report"),
                        new Arguments.Option("help", "option to print this message")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println(ae);
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("help")) {
            printUsage(arguments);
            System.exit(0);
        }

        int burnin = -1;
        if (arguments.hasOption("burnin")) {
            burnin = arguments.getIntegerOption("burnin");
        }

        boolean shortReport = arguments.hasOption("short");

        String exportFileName = null;
        if (arguments.hasOption("export")) {
            exportFileName = arguments.getStringOption("export");
        }

        double minSupport = 0.0;
        if (arguments.hasOption("limit")) {
            minSupport = arguments.getRealOption("limit");
        }

        double credibleSetProbability = 0.95;
        if (arguments.hasOption("probability")) {
            credibleSetProbability = arguments.getRealOption("probability");
        }

        int maxExport = -1;
        if (arguments.hasOption("max")) {
            maxExport = arguments.getIntegerOption("max");
        }

        String inputFileName = null;
        String trueTreeFileName = null;
        String outputFileName = null;

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 3) {
            System.err.println("Unknown option: " + args2[2]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (args2.length > 0) {
            inputFileName = args2[0];
        }

        if (args2.length > 1) {
            trueTreeFileName = args2[1];
        }

        if (args2.length > 2) {
            outputFileName = args2[2];
        }

        if (inputFileName == null) {
            // No input file name was given so throw up a dialog box...
            inputFileName = Utils.getLoadFileName("TreeLogAnalyser " + version.getVersionString() + " - Select log file to analyse");
        }

        new TreeLogAnalyser(burnin, inputFileName, outputFileName, trueTreeFileName, exportFileName,
                minSupport, credibleSetProbability, maxExport, !shortReport);

        System.exit(0);
    }
}

