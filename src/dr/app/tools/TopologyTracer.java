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
import dr.evolution.tree.treemetrics.*;
import dr.util.Version;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * @author Guy Baele
 * @author Andrew Rambaut
 */
public class TopologyTracer {

    private final static Version version = new BeastVersion();

    private final static boolean PROFILE = true;

    private static final String STATE = "state";
    private static final String RFDISTANCE = "RFdistance";
    private static final String CLADE_HEIGHT = "cladeHeight";
    private static final String BRANCH_SCORE_METRIC = "rootedBranchScoreMetric";
    private static final String PATH_DIFFERENCE = "SteelPenny";
    private static final String KC_METRIC = "KCmetric";

    // output to stdout
    private static PrintStream progressStream = System.out;

    public TopologyTracer(int burnin, String treeFile, String userProvidedTreeFile, String outputFile, ArrayList<Double> lambdaValues) {

        try {

            long startTime = System.currentTimeMillis();

            progressStream.println("Reading & processing trees ...");

            BufferedReader reader = new BufferedReader(new FileReader(treeFile));

            String line = reader.readLine();

            TreeImporter importer;
            if (line.toUpperCase().startsWith("#NEXUS")) {
                importer = new NexusImporter(reader);
            } else {
                importer = new NewickImporter(reader);
            }

            Tree focalTree = null;
            boolean userProvidedTree = false;
            if (!userProvidedTreeFile.equals("")) {
                userProvidedTree = true;
                progressStream.println("User-provided focal tree.");
                //get tree from user provided tree file
                BufferedReader focalReader = new BufferedReader(new FileReader(userProvidedTreeFile));
                TreeImporter userImporter;
                String userLine = focalReader.readLine();
                if (userLine.toUpperCase().startsWith("#NEXUS")) {
                    userImporter = new NexusImporter(focalReader);
                } else {
                    userImporter = new NewickImporter(focalReader);
                }
                focalTree = userImporter.importNextTree();
            } else {
                //pick first tree as focal tree
                focalTree = importer.importNextTree();
            }

            List<TreeMetric> treeMetrics = new ArrayList<TreeMetric>();
            treeMetrics.add(new RobinsonFouldsMetric());
            treeMetrics.add(new CladeHeightMetric());
            treeMetrics.add(new RootedBranchScoreMetric());
            treeMetrics.add(new SteelPennyPathDifferenceMetric());
            for (double lambda:  lambdaValues) {
                treeMetrics.add(new KendallColijnPathDifferenceMetric(lambda));
            }


            List<Long> treeStates = new ArrayList<Long>();
            List<String> treeIds = new ArrayList<String>();
            List<List<Double>> metricValues = new ArrayList<List<Double>>();
            for (TreeMetric treeMetric : treeMetrics) {
                metricValues.add(new ArrayList<Double>());
            }

            if (!userProvidedTree) {
                //take into account first distance of focal tree to itself
                treeIds.add(focalTree.getId());
                treeStates.add((long) 0);

                int i = 0;
                for (TreeMetric treeMetric : treeMetrics) {
                    metricValues.get(i).add(treeMetric.getMetric(focalTree, focalTree));
                    i++;
                }

            }

            int numberOfTrees = 1;

            long[] timings = new long[6];
            long beforeTime, afterTime;

            while (importer.hasTree()) {

                //no need to keep trees in memory
                Tree tree = importer.importNextTree();
                treeIds.add(tree.getId());
                treeStates.add(Long.parseLong(tree.getId().split("_")[1]));

                int i = 0;
                for (TreeMetric treeMetric : treeMetrics) {
                    metricValues.get(i).add(treeMetric.getMetric(focalTree, tree));
                    i++;
                }

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
            writer.write(STATE);
            for (TreeMetric treeMetric : treeMetrics) {
                writer.write("\t" + treeMetric);
            }
            writer.write("\n");

            for (int i = 0; i < treeStates.size(); i++) {
                writer.write(Long.toString(treeStates.get(i)));
                for (List<Double> values : metricValues) {
                    writer.write("\t" + values.get(i));
                }
                writer.write("\n");
            }

            progressStream.println("Done.");

            long endTime = System.currentTimeMillis();

            progressStream.println("\nAnalyzed " + treeStates.size() + " trees, took " + (endTime-startTime)/1000.0 + " seconds.\n");

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
                        new Arguments.StringOption("tree", "tree file name", "a focal tree provided by the user [default = first tree in .trees file]"),
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

        String providedFileName = "";
        if (arguments.hasOption("tree")) {
            providedFileName = arguments.getStringOption("tree");
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

        new TopologyTracer(burnin, inputFileName, providedFileName, outputFileName, lambdaValues);

        System.exit(0);


    }

}
