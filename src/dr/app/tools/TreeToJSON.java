/*
 * RootToTip.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.tools;

import dr.app.beauti.options.DateGuesser;
import dr.app.pathogen.TemporalRooting;
import dr.app.util.Arguments;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.stats.Regression;
import dr.stats.Variate;
import dr.util.Version;

import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/*
 * Essentially a command line version of PathOgen. Written to
 * perform the analysis on sets of trees.
 * 
 * @author Andrew Rambaut
 */

public class TreeToJSON {

    private final static Version version = new Version() {
        @Override
        public String getVersion() {
            return "1.0";
        }

        @Override
        public String getVersionString() {
            return "v1.0";
        }

        @Override
        public String getBuildString() {
            return "";
        }

        @Override
        public String getDateString() {
            return "2015";
        }

        @Override
        public String[] getCredits() {
            return new String[0];
        }

        @Override
        public String getHTMLCredits() {
            return "";
        }
    };

    public TreeToJSON(String dateOrder, String inputFileName, String outputFileName) throws IOException {

        System.out.println("Reading tree(s)...");

        boolean firstTree = true;
        FileReader fileReader = new FileReader(inputFileName);
        TreeImporter importer = new NexusImporter(fileReader);

        Tree tree = null;

        DateGuesser dg = new DateGuesser();
        dg.fromLast = false;
        if (dateOrder.equals("FIRST")) {
            dg.order = 0;
        } else if (dateOrder.equals("LAST")) {
            dg.order = 0;
            dg.fromLast = true;
        } else {
            dg.order = Integer.parseInt(dateOrder) - 1;
            if (dg.order < 0 || dg.order > 100) {
                System.err.println("Error Parsing order of date field: " + dateOrder);
            }
        }

        TaxonList taxa = null;

        try {
            while (importer.hasTree()) {
                Tree tree1 = importer.importNextTree();

                if (firstTree) {
                    taxa = tree1;
                    tree = tree1;

                    dg.guessDates(taxa);

                    firstTree = false;
                }

            }
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }
        fileReader.close();

        System.out.println("Writing tree(s)...");

        PrintStream printStream;

        if (outputFileName != null) {
            printStream = new PrintStream(outputFileName);
        } else {
            printStream = new PrintStream(System.out);
        }


        AuspiceJSONExporter exporter = new AuspiceJSONExporter(printStream);
        exporter.exportTree(tree);

        printStream.close();

        System.out.println("Done.");

    }

    int totalTrees = 0;
    int totalTreesUsed = 0;

    public static void printTitle() {
        System.out.println();
        centreLine("TreeToJSON " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("Tree file to JSON format converter", 60);
        centreLine("by", 60);
        centreLine("Andrew Rambaut", 60);
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

        arguments.printUsage("treetojson", "<input-file-name> [<output-file-name>]");
        System.out.println();
        System.out.println("  Example: treetojson test.tree test.json");
        System.out.println();
    }

    //Main method
    public static void main(String[] args) throws IOException {

        String inputFileName = null;
        String outputFileName = null;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.StringOption("dateorder", "date_order", "order of date field in taxon name: first, last, 1, 2 etc. [default = last]"),
                        new Arguments.Option("help", "option to print this message"),
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

        String dateOrder = "LAST";
        if (arguments.hasOption("dateorder")) {
            dateOrder = arguments.getStringOption("dateorder").toUpperCase();
        }

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 2) {
            System.err.println("Unknown option: " + args2[2]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (args2.length == 0) {
            System.err.println("Missing input file name");
            printUsage(arguments);
            System.exit(1);
        }


        inputFileName = args2[0];
        if (args2.length == 2) {
            outputFileName = args2[1];
        }

        new TreeToJSON(dateOrder,
                inputFileName,
                outputFileName
        );

        System.exit(0);
    }

}