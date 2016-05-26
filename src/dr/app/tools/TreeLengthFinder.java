/*
 * TreeLengthFinder.java
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
import dr.app.treestat.statistics.TreeLength;
import dr.inference.trace.TraceException;
import dr.util.Version;
import dr.evolution.io.TreeImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.Importer;
import dr.evolution.tree.Tree;

import java.io.*;

/**
 * @author Wai Lok Sibon Li
 */

public class TreeLengthFinder {

    private final static Version version = new BeastVersion();


    public TreeLengthFinder(int burnin, String inputFileName, String outputFileName/*, boolean verbose,
                       boolean hpds, boolean ess, boolean stdErr,
                       String marginalLikelihood*/) throws IOException, TraceException {

        File parentFile = new File(inputFileName);

        if (parentFile.isDirectory()) {
            System.out.println("Analysing all trees files below directory: " + inputFileName);
        } else if (parentFile.isFile()) {
            System.out.println("Analysing tree file: " + inputFileName);
        } else {
            System.err.println(inputFileName + " does not exist!");
            System.exit(0);
        }

        if (outputFileName != null) {
            FileOutputStream outputStream = new FileOutputStream(outputFileName);
            System.setOut(new PrintStream(outputStream));
        }

        analyze(parentFile, burnin/*, verbose, new boolean[]{true}, hpds, ess, stdErr, marginalLikelihood*/);
    }

    /**
     * Recursively analyzes log files.
     *
     * @param file       the file to analyze (if this is a directory then the files within it are analyzed)
     * @param burnin     the burnin to use
     * @throws dr.inference.trace.TraceException
     *          if the trace file is in the wrong format or corrupted
     */
    private void analyze(File file, int burnin) throws TraceException {

        if (file.isFile()) {
            try {

                String name = file.getCanonicalPath();

                report(name, burnin);
                //TraceAnalysis.report(name, burnin, marginalLikelihood);

            } catch (IOException e) {
                //e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        } else {
            File[] files = file.listFiles();
            for (File f : files) {
                if (f.isDirectory()) {
                    analyze(f, burnin);
                //} else if (f.getName().endsWith(".tree") || f.getName().endsWith(".p")) {
                } else if (f.getName().endsWith(".trees")) {
                    analyze(f, burnin);
                }
            }
        }
    }


    /**
     * Recursively analyzes trees files.
     *
     * @param name       the file to analyze (if this is a directory then the files within it are analyzed)
     * @param burnin     the burnin to use
     */
    private void report(String name, int burnin) {
        double treeLength = 0.0;
        int count = 0;

        try {
            FileReader fileReader = new FileReader(new File(name));
            TreeImporter importer = new NexusImporter(fileReader);
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();
                if(count>=burnin) {
                    treeLength += TreeLength.FACTORY.createStatistic().getSummaryStatistic(tree)[0];

                }
                count++;
            }
            treeLength /= (count-burnin);
            System.out.println(name + "\t" + burnin + "\t" + treeLength);
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
        } catch (IOException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
        }
    }

    public static void printTitle() {
        System.out.println();
        centreLine("TreeLengthFinder " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("MCMC Output analysis", 60);
        centreLine("by", 60);
        centreLine("Andrew Rambaut and Alexei J. Drummond", 60);
        System.out.println();
        centreLine("Institute of Evolutionary Biology", 60);
        centreLine("University of Edinburgh", 60);
        centreLine("a.rambaut@ed.ac.uk", 60);
        System.out.println();
        centreLine("Department of Computer Science", 60);
        centreLine("University of Auckland", 60);
        centreLine("alexei@cs.auckland.ac.nz", 60);
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

        //arguments.printUsage("loganalyser", "[-burnin <burnin>] [-short][-hpd] [-std] [<input-file-name> [<output-file-name>]]");
        arguments.printUsage("treelengthfinder", "[-burnin <burnin>][<input-file-name> [<output-file-name>]]");
        System.out.println();
        System.out.println("  Example: treelengthfinder test.tree");
        System.out.println("  Example: treelengthfinder -burnin 10000 trees.log out.txt");
        System.out.println();

    }

    //Main method
    public static void main(String[] args) throws IOException, TraceException {

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in'"),
                        /*new Arguments.Option("short", "use this option to produce a short report"),
                        new Arguments.Option("hpd", "use this option to produce hpds for each trace"),
                        new Arguments.Option("ess", "use this option to produce ESSs for each trace"),
                        new Arguments.Option("stdErr", "use this option to produce standard Error"),
                        new Arguments.StringOption("marginal", "trace_name", "specify the trace to use to calculate the marginal likelihood"),*/
//				new Arguments.Option("html", "format output as html"),
//				new Arguments.Option("svg", "generate svg graphics"),
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




        /*boolean hpds = arguments.hasOption("hpd");
        boolean ess = arguments.hasOption("ess");
        boolean stdErr = arguments.hasOption("stdErr");
        boolean shortReport = arguments.hasOption("short");

        String marginalLikelihood = null;
        if (arguments.hasOption("marginal")) {
            marginalLikelihood = arguments.getStringOption("marginal");
        }*/

        String inputFileName = null;
        String outputFileName = null;

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 2) {
            System.err.println("Unknown option: " + args2[2]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (args2.length > 0) {
            inputFileName = args2[0];
        }
        if (args2.length > 1) {
            outputFileName = args2[1];
        }

        if (inputFileName == null) {
            // No input file name was given so throw up a dialog box...
            inputFileName = Utils.getLoadFileName("TreeLengthFinder " + version.getVersionString() + " - Select tree file to analyse");
        }

        if(burnin==-1) {
            System.out.println("Enter number of trees to burn-in (integer): ");
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            burnin = Integer.parseInt(br.readLine());
        }

        new TreeLengthFinder(burnin, inputFileName, outputFileName/*, !shortReport, hpds, ess, stdErr, marginalLikelihood*/);

        System.exit(0);
    }
}