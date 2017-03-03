/*
 * TopologyTracer.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.*;
import dr.util.Version;
import jebl.evolution.treemetrics.BilleraMetric;
import jebl.evolution.treemetrics.CladeHeightMetric;
import jebl.evolution.treemetrics.RobinsonsFouldMetric;

import java.io.*;
import java.util.ArrayList;
import java.util.Locale;

/**
 * @author Guy Baele
 */
public class TopologyTracer {

    private final static Version version = new BeastVersion();

    private static final String STATE = "state";
    private static final String RFDISTANCE = "RFdistance (pseudo)";
    private static final String BILLERA_METRIC = "BilleraMetric";
    private static final String CLADE_HEIGHT = "cladeHeight";
    private static final String BRANCH_SCORE_METRIC = "branchScoreMetric";
    private static final String PATH_DIFFERENCE = "pathDifference (pseudo)";
    private static final String KC_METRIC = "KCmetric";

    // output to stdout
    private static PrintStream progressStream = System.out;

    public TopologyTracer(int burnin, String treeFile, String outputFile, ArrayList<Double> lambdaValues) {

        try {

            progressStream.println("Reading & processing trees ...");

            BufferedReader reader = new BufferedReader(new FileReader(treeFile));

            String line = reader.readLine();

            TreeImporter importer;
            if (line.toUpperCase().startsWith("#NEXUS")) {
                importer = new NexusImporter(reader);
            } else {
                importer = new NewickImporter(reader);
            }

            //pick first tree as focal tree
            Tree focalTree = importer.importNextTree();

            ArrayList<Long> treeStates = new ArrayList<Long>();
            ArrayList<String> treeIds = new ArrayList<String>();

            ArrayList<Double> jeblRFDistances = new ArrayList<Double>();
            ArrayList<Double> billeraMetric = new ArrayList<Double>();
            ArrayList<Double> cladeHeightMetric = new ArrayList<Double>();
            ArrayList<Double> branchScoreMetric = new ArrayList<Double>();
            ArrayList<Double> pathDifferenceMetric = new ArrayList<Double>();
            ArrayList<ArrayList<Double>> kcMetrics = new ArrayList();
            for (int i = 0; i < lambdaValues.size(); i++) {
                kcMetrics.add(new ArrayList<Double>());
            }

            //take into account first distance of focal tree to itself
            treeStates.add((long)0);

            jeblRFDistances.add(new RobinsonsFouldMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(focalTree)));
            billeraMetric.add(new BilleraMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(focalTree)));
            cladeHeightMetric.add(new CladeHeightMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(focalTree)));
            branchScoreMetric.add(new BranchScoreMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(focalTree)));
            pathDifferenceMetric.add(new SPPathDifferenceMetric().getMetric(focalTree, focalTree));
            ArrayList<Double> allKCMetrics = (new KCPathDifferenceMetric().getMetric(focalTree, focalTree, lambdaValues));
            for (int i = 0; i < allKCMetrics.size(); i++) {
                kcMetrics.get(i).add(allKCMetrics.get(i));
            }

            int numberOfTrees = 1;

            while (importer.hasTree()) {

                //no need to keep trees in memory
                Tree tree = importer.importNextTree();
                treeIds.add(tree.getId());
                treeStates.add(Long.parseLong(tree.getId().split("_")[1]));

                //TODO Does the BEAST/JEBL code report half the RF distance?
                jeblRFDistances.add(new RobinsonsFouldMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(tree)));
                billeraMetric.add(new BilleraMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(tree)));
                cladeHeightMetric.add(new CladeHeightMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(tree)));
                branchScoreMetric.add(new BranchScoreMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(tree)));
                pathDifferenceMetric.add(new SPPathDifferenceMetric().getMetric(focalTree, tree));
                allKCMetrics = (new KCPathDifferenceMetric().getMetric(focalTree, tree, lambdaValues));
                for (int i = 0; i < allKCMetrics.size(); i++) {
                    kcMetrics.get(i).add(allKCMetrics.get(i));
                }

                //TODO Last tree is not being processed?
                //System.out.println(tree.getId());

                numberOfTrees++;

                if (numberOfTrees % 25 == 0) {
                    progressStream.println(numberOfTrees + " trees parsed ...");
                    progressStream.flush();
                }

            }

            progressStream.println("\nWriting log file ...");

            BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
            BeastVersion version = new BeastVersion();
            writer.write("# BEAST " + version.getVersionString() + "\n");
            writer.write(STATE + "\t" + RFDISTANCE + "\t" + BILLERA_METRIC + "\t" +
                    BRANCH_SCORE_METRIC + "\t" + CLADE_HEIGHT + "\t");
            for (Double l : lambdaValues) {
                writer.write(KC_METRIC + "-" + l + "\t");
            }
            writer.write(PATH_DIFFERENCE + "\n");

            for (int i = 0; i < treeStates.size(); i++) {
                writer.write(treeStates.get(i) + "\t");
                writer.write(jeblRFDistances.get(i) + "\t");
                writer.write(billeraMetric.get(i) + "\t");
                writer.write(branchScoreMetric.get(i) + "\t");
                writer.write(cladeHeightMetric.get(i) + "\t");
                for (int j = 0; j < lambdaValues.size(); j++) {
                    writer.write(kcMetrics.get(j).get(i) + "\t");
                }
                writer.write(pathDifferenceMetric.get(i) + "\n");
            }

            progressStream.println("Done.");
            progressStream.flush();
            progressStream.close();

            writer.flush();
            writer.close();

        } catch (FileNotFoundException fnf) {
            System.err.println(fnf);
        } catch (IOException ioe) {
            System.err.println(ioe);
        } catch (Importer.ImportException ime) {
            System.err.println(ime);
        }

    }

    public static void printUsage(Arguments arguments) {
        arguments.printUsage("TopologyTracer", "<input-file-name> <output-file-name>");
        System.out.println();
        System.out.println("  Example: treeloganalyser test.trees ess-values.log");
        System.out.println();
    }

    public static void main(String[] args) {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in' [default = none]"),
                        new Arguments.RealOption("lambda", "the lambda value to be used for the 'Kendall-Colijn metric' [default = {0,0.5,1}]"),
                        new Arguments.Option("help", "option to print this message")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println(ae);
            printUsage(arguments);
            System.exit(1);
        }

        int burnin = 0;

        if (arguments.hasOption("burnin")) {
            burnin = arguments.getIntegerOption("burnin");
        }

        ArrayList<Double> lambdaValues = new ArrayList<Double>();
        lambdaValues.add(0.0);
        lambdaValues.add(0.5);
        lambdaValues.add(1.0);

        if (arguments.hasOption("lambda")) {
            lambdaValues.add(arguments.getRealOption("lambda"));
        }

        String inputFileName = null;
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
            outputFileName = args2[1];
        }

        if (inputFileName == null) {
            // No input file name was given so throw up a dialog box...
            inputFileName = Utils.getLoadFileName("TopologyTracer " + version.getVersionString() + " - Select log file to analyse");
        }

        new TopologyTracer(burnin, inputFileName, outputFileName, lambdaValues);

        System.exit(0);


    }

}
