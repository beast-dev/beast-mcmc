/*
 * LogAnalyser.java
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
import dr.app.util.Utils;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceAnalysis;
import dr.inference.trace.TraceException;
import dr.util.Version;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Locale;

import static dr.inference.trace.TraceAnalysis.analyzeLogFile;

public class OperatorAnalyser {

    private final static Version version = new BeastVersion();


    public OperatorAnalyser(long burnin, String opsFileName, String logFileName, String outputFileName) throws IOException, TraceException {

        if (outputFileName != null) {
            FileOutputStream outputStream = new FileOutputStream(outputFileName);
            System.setOut(new PrintStream(outputStream));
        }

        analyze(opsFileName, logFileName, burnin);
    }


    /**
     * Analyzes ops/log files.
     */
    private void analyze(String opsFileName, String logFileName, long burnin) throws TraceException, IOException {
        LogFileTraces traces = analyzeLogFile(logFileName, burnin);
    }

    public static void printTitle() {
        System.out.println();
        centreLine("OpsAnalyser " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("BEAST operator performance analysis", 60);
        centreLine("by", 60);
        centreLine("Andrew Rambaut", 60);
        System.out.println();
        centreLine("Institute of Ecology and Evolution", 60);
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

        arguments.printUsage("opanalyser", "[-burnin <burnin>] <operator-file-name> <log-file-name>");
        System.out.println();
        System.out.println("  Example: opanalyser test.ops test.log");
        System.out.println("  Example: opanalyser -burnin 10000 test.ops test.log out.txt");
        System.out.println();

    }

    //Main method
    public static void main(String[] args) throws IOException, TraceException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in'"),
//                        new Arguments.Option("short", "use this option to produce a short report"),
//                        new Arguments.Option("hpd", "use this option to produce hpds for each trace"),
//                        new Arguments.Option("ess", "use this option to produce ESSs for each trace"),
//                        new Arguments.Option("stdErr", "use this option to produce standard Error"),
//                        new Arguments.StringOption("marginal", "trace_name", "specify the trace to use to calculate the marginal likelihood"),
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

//        boolean hpds = arguments.hasOption("hpd");
//        boolean ess = arguments.hasOption("ess");
//        boolean stdErr = arguments.hasOption("stdErr");
//        boolean shortReport = arguments.hasOption("short");

        String opsFileName = null;
        String logFileName = null;
        String outputFileName = null;

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 3) {
            System.err.println("Unknown option: " + args2[3]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (args2.length > 0) {
            opsFileName = args2[0];
        }
        if (args2.length > 1) {
            logFileName = args2[1];
        }
        if (args2.length > 2) {
            outputFileName = args2[2];
        }

            new OperatorAnalyser(burnin, opsFileName, logFileName, outputFileName);

        System.exit(0);
    }
}

