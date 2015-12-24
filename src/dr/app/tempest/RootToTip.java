/*
 * RootToTip.java
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

package dr.app.tempest;

import dr.app.beauti.options.DateGuesser;
import dr.app.util.Arguments;
import dr.app.tools.NexusExporter;
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
 * Essentially a command line version of TempEst. Written to
 * perform the analysis on sets of trees.
 * 
 * @author Andrew Rambaut
 */

public class RootToTip {

    private final static Version version = new Version() {
        @Override
        public String getVersion() {
            return "1.5";
        }

        @Override
        public String getVersionString() {
            return "v1.5";
        }

        @Override
        public String getBuildString() {
            return "";
        }

        @Override
        public String getDateString() {
            return "2003-2015";
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

    public RootToTip(int burnin, String dateOrder, final boolean keepRoot, String outgroup,
                     boolean writeTree, String inputFileName, String outputFileName) throws IOException {

        System.out.println("Reading tree(s)...");

        boolean firstTree = true;
        FileReader fileReader = new FileReader(inputFileName);
        TreeImporter importer = new NexusImporter(fileReader);

        List<Regression> regressions = new ArrayList<Regression>();
        List<Tree> trees = new ArrayList<Tree>();

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
        TemporalRooting temporalRooting = null;

        try {
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (firstTree) {
                    taxa = tree;

                    dg.guessDates(taxa);

                    temporalRooting = new TemporalRooting(taxa);

                    firstTree = false;
                }

                if (totalTrees >= burnin) {
                    Tree rootedTree = tree;

                    if (!keepRoot) {
                        rootedTree = temporalRooting.findRoot(tree, TemporalRooting.RootingFunction.CORRELATION);
                    }

                    regressions.add(temporalRooting.getRootToTipRegression(rootedTree));

                    if (writeTree) {
                        trees.add(rootedTree);
                    }
                    totalTreesUsed += 1;
                }
                totalTrees += 1;

            }
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }
        fileReader.close();

        PrintWriter printWriter;

        if (!writeTree && outputFileName != null) {
            printWriter = new PrintWriter(outputFileName);
        } else {
            printWriter = new PrintWriter(System.out);
        }

        if (regressions.size() == 1) {
            Regression r = regressions.get(0);

            Variate dates = r.getXData();
            Variate distances = r.getYData();

            printWriter.println("date\tdistance");
            for (int i = 0; i < dates.getCount(); i++) {
                printWriter.println(dates.get(i) + "\t" + distances.get(i));
            }
            printWriter.println();
            printWriter.println("Regression slope = " + r.getGradient());
            printWriter.println("X-Intercept = " + r.getXIntercept());
            printWriter.println("Y-Intercept = " + r.getYIntercept());
            printWriter.println("Residual mean squared = " + r.getResidualMeanSquared());
            printWriter.println("R^2 = " + r.getRSquared());
            printWriter.println("Correlation coefficient = " + r.getCorrelationCoefficient());

        } else {
            printWriter.println("tree\tslope\tx-intercept\ty-intercept\tcorrelation");
            int i = 1;
            for (Regression r : regressions) {
                printWriter.print("\t" + i);
                printWriter.print("\t" + r.getGradient());
                printWriter.print("\t" + r.getXIntercept());
                printWriter.println("\t" + r.getYIntercept());
                printWriter.println("\t" + r.getCorrelationCoefficient());
            }

        }

        printWriter.close();

        if (writeTree) {
            PrintStream printStream;

            if (outputFileName != null) {
                printStream = new PrintStream(outputFileName);
            } else {
                printStream = new PrintStream(System.out);
            }

            NexusExporter exporter = new NexusExporter(printStream);
            Tree[] treeArray = new Tree[trees.size()];
            trees.toArray(treeArray);

            exporter.exportTrees(treeArray);

            printStream.close();
        }
    }

    int totalTrees = 0;
    int totalTreesUsed = 0;

    public static void printTitle() {
        System.out.println();
        centreLine("RootToTip " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("Root to tip distance vs. time of sampling", 60);
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

        arguments.printUsage("roottotip", "<input-file-name> [<output-file-name>]");
        System.out.println();
        System.out.println("  Example: roottotip -burnin 100 test.trees rootToTip.txt");
        System.out.println();
    }

    //Main method
    public static void main(String[] args) throws IOException {

        String inputFileName = null;
        String outputFileName = null;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("burnin", "the number of trees to be ignored as 'burn-in' [default = 0]"),
                        new Arguments.StringOption("dateorder", "date_order", "order of date field in taxon name: first, last, 1, 2 etc. [default = last]"),
//                        new Arguments.StringOption("outgroup", "{taxon list}", "one or more taxa that will be used to root the tree(s) [default = find root]"),
                        new Arguments.Option("keeproot", "keep the existing root of the input trees [default = estimate root]"),
                        new Arguments.Option("writetree", "Write the optimally rooted tree to the output file"),
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

        int burnin = 0;
        if (arguments.hasOption("burnin")) {
            burnin = arguments.getIntegerOption("burnin");
        }

        String dateOrder = "LAST";
        if (arguments.hasOption("dateorder")) {
            dateOrder = arguments.getStringOption("dateorder").toUpperCase();
        }

        String outgroup = null;
        if (arguments.hasOption("outgroup")) {
            outgroup = arguments.getStringOption("dateorder").toUpperCase();
        }

        boolean keepRoot = arguments.hasOption("keeproot");

        boolean writeTree = arguments.hasOption("writetree");

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

        new RootToTip(burnin,
                dateOrder,
                keepRoot,
                outgroup,
                writeTree,
                inputFileName,
                outputFileName
        );

        System.exit(0);
    }

}