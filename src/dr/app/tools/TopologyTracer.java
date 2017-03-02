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
import dr.evolution.tree.BranchScoreMetric;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeMetrics;
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
    private static final String RFDISTANCE = "RFdistance";
    //private static final String JEBLRFDISTANCE = "jeblRFdistance";
    private static final String BILLERA_METRIC = "BilleraMetric";
    private static final String CLADE_HEIGHT = "cladeHeight";
    private static final String BRANCH_SCORE_METRIC = "branchScoreMetric";

    // Messages to stderr, output to stdout
    private static PrintStream progressStream = System.out;

    //TODO Nothing is being done with the burnin argument as of yet; let Tracer take care of it
    public TopologyTracer(int burnin, String treeFile, String outputFile) {

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

            //ArrayList<Double> rfDistances = new ArrayList<Double>();
            ArrayList<Double> jeblRFDistances = new ArrayList<Double>();
            ArrayList<Double> billeraMetric = new ArrayList<Double>();
            ArrayList<Double> cladeHeightMetric = new ArrayList<Double>();
            ArrayList<Double> branchScoreMetric = new ArrayList<Double>();

            //TODO Still need mean path difference metric?

            //take into account first distance of focal tree to itself
            treeStates.add((long)0);

            //rfDistances.add(TreeMetrics.getRobinsonFouldsDistance(focalTree, focalTree));
            jeblRFDistances.add(new RobinsonsFouldMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(focalTree)));
            billeraMetric.add(new BilleraMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(focalTree)));
            cladeHeightMetric.add(new CladeHeightMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(focalTree)));
            branchScoreMetric.add(new BranchScoreMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(focalTree)));

            int numberOfTrees = 1;

            while (importer.hasTree()) {

                //no need to keep trees in memory
                Tree tree = importer.importNextTree();
                treeIds.add(tree.getId());
                treeStates.add(Long.parseLong(tree.getId().split("_")[1]));
                //TODO Does the BEAST/JEBL code report half the RF distance?
                //rfDistances.add(TreeMetrics.getRobinsonFouldsDistance(focalTree, tree));

                jeblRFDistances.add(new RobinsonsFouldMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(tree)));

                billeraMetric.add(new BilleraMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(tree)));

                cladeHeightMetric.add(new CladeHeightMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(tree)));

                branchScoreMetric.add(new BranchScoreMetric().getMetric(Tree.Utils.asJeblTree(focalTree), Tree.Utils.asJeblTree(tree)));

                //TODO Last tree is not being processed?
                //System.out.println(tree.getId());

                numberOfTrees++;

                if (numberOfTrees % 25 == 0) {
                    progressStream.println(numberOfTrees + " trees parsed ...");
                    progressStream.flush();
                }

            }

            /*for (int i = 0; i < treeIds.size(); i++) {
                System.out.println(treeIds.get(i) + " " + distances.get(i));
            }*/

            progressStream.println("\nWriting log file ...");

            BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
            BeastVersion version = new BeastVersion();
            writer.write("# BEAST " + version.getVersionString() + "\n");
            writer.write(STATE + "\t" + RFDISTANCE + "\t" + BILLERA_METRIC + "\t" +
                    BRANCH_SCORE_METRIC + "\t" + CLADE_HEIGHT + "\n");

            for (int i = 0; i < treeStates.size(); i++) {

                writer.write(treeStates.get(i) + "\t");
                //writer.write(rfDistances.get(i) + "\t");
                writer.write(jeblRFDistances.get(i) + "\t");
                writer.write(billeraMetric.get(i) + "\t");
                writer.write(branchScoreMetric.get(i) + "\t");
                writer.write(cladeHeightMetric.get(i) + "\n");
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

        new TopologyTracer(burnin, inputFileName, outputFileName);

        System.exit(0);


    }

}
