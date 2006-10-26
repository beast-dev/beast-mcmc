/*
 * LogAnalyser.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
import dr.inference.trace.TraceAnalysis;
import dr.util.Version;

import java.io.*;

public class LogAnalyser {

	private final static Version version = new BeastVersion();

	public LogAnalyser(int burnin, String inputFileName, String outputFileName, boolean verbose) throws java.io.IOException {

		File inputFile = new File(inputFileName);

        if (inputFile.isDirectory()) {
            System.out.println("Analysing all log files below directory: " + inputFileName);
        } else if (inputFile.isFile()) {
            System.out.println("Analysing log file: " + inputFileName);
        } else {
            System.err.println(inputFileName + " does not exist!");
            System.exit(0);
        }

		if (outputFileName != null) {
			FileOutputStream outputStream = new FileOutputStream(outputFileName);
			System.setOut(new PrintStream(outputStream));
		}

        analyze(inputFile, burnin, verbose, new boolean[] {true});
	}

    /**
     * Recursively analyzes log files.
     * @param file the file to analyze (if this is a directory then the files within it are analyzed)
     * @param burnin the burnin to use
     * @param verbose if true then a full report is done on each log file, otherwise only a single line report is made
     * @param drawHeader if boolean value in the zeroth position of this array is true then a head is drawn for the short reports.
     */
    private static void analyze(File file, int burnin, boolean verbose, boolean[] drawHeader) {

        if (file.isFile()) {
            try {

                if (verbose) {
                    TraceAnalysis.report(new FileReader(file), burnin);
                } else {
                    TraceAnalysis.shortReport(file.getName(), new FileReader(file), burnin, drawHeader[0]);
                    drawHeader[0] = false;
                }
            } catch (IOException e) {
                //e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        } else {
            File[] files = file.listFiles();
            for (int i = 0; i < files.length; i++) {
                if (files[i].isDirectory()) {
                    analyze(files[i], burnin, verbose, drawHeader);
                } else if (files[i].getName().endsWith(".log") || files[i].getName().endsWith(".p")) {
                    analyze(files[i], burnin, verbose, drawHeader);
                } else {
                    if (verbose) System.out.println("Ignoring file: " + files[i]);
                }
            }
        }
    }

    public static void printTitle() {
        System.out.println();
        centreLine("LogAnalyser v1.4, 2002-2006", 60);
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
        for (int i = 0; i < n1; i++) { System.out.print(" "); }
        System.out.println(line);
    }

 
	public static void printUsage(Arguments arguments) {

		arguments.printUsage("loganalyser", "[-burnin <burnin>] [-short] [<input-file-name> [<output-file-name>]]");
		System.out.println();
		System.out.println("  Example: loganalyser test.log");
		System.out.println("  Example: loganalyser -burnin 10000 trees.log out.txt");
		System.out.println();

	}

	//Main method
	public static void main(String[] args) throws java.io.IOException {

		printTitle();

		Arguments arguments = new Arguments(
			new Arguments.Option[] {
				new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in'"),
				new Arguments.Option("short", "use this option to produce a short report"),
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

		boolean shortReport = arguments.hasOption("short");

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
			inputFileName = Utils.getLoadFileName("LogAnalyser v1.4 - Select log file to analyse");
		}

		new LogAnalyser(burnin, inputFileName, outputFileName, !shortReport);

		System.exit(0);
	}
}

