/*
 * CoalGenApp.java
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

package dr.app.coalgen;

import dr.evolution.io.Importer;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import jam.framework.SingleDocApplication;
import jebl.evolution.coalescent.EmpiricalDemographicFunction;
import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NewickExporter;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.io.TreeImporter;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.treesimulation.CoalescentIntervalGenerator;
import jebl.evolution.treesimulation.IntervalGenerator;
import jebl.evolution.treesimulation.TreeSimulator;

import javax.swing.*;
import java.io.*;
import java.util.Arrays;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class CoalGenApp {
    public CoalGenApp() {

        try {

            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

            java.net.URL url = CoalGenApp.class.getResource("/images/CoalGen.png");
            Icon icon = null;

            if (url != null) {
                icon = new ImageIcon(url);
            }

            String nameString = "CoalGen";
            String aboutString = "Coalescent Tree Simulator\nVersion 2.0\n \nCopyright 2004-2009 Andrew Rambaut and Alexei Drummond\nAll Rights Reserved.";


            SingleDocApplication app = new SingleDocApplication(new CoalGenMenuFactory(), nameString, aboutString, icon);


            CoalGenFrame frame = new CoalGenFrame(nameString);
            app.setDocumentFrame(frame);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void simulate(String inputFileName, String treesFileName, String outputFileName) throws IOException, TraceException, Importer.ImportException {
        File logFile = new File(inputFileName);

        System.out.println("Loading trace file: " + inputFileName);
        LogFileTraces traces = new LogFileTraces(inputFileName, logFile);
        traces.loadTraces();
        traces.setBurnIn(0);
        System.out.println(traces.getStateCount() + " states loaded");

        System.out.println();
        System.out.println("Opening trees file: " + treesFileName);
//        BufferedReader reader = new BufferedReader(new FileReader(treesFileName));

        System.out.println("Simulating...");
        System.out.println("0              25             50             75            100");
        System.out.println("|--------------|--------------|--------------|--------------|");

        int stepSize = traces.getStateCount() / 60;
        if (stepSize < 1) stepSize = 1;

        PrintWriter writer = new PrintWriter(new FileWriter(outputFileName));

        FileReader fileReader = new FileReader(treesFileName);
        TreeImporter importer = new NexusImporter(fileReader);

            EmpiricalDemographicFunction demo = null;
            IntervalGenerator intervals = null;
            TreeSimulator sim = null;

            CoalGenData data = new CoalGenData();
            data.traces = traces;
            data.setDemographicModel(7); // const stepwise bsp
            data.setupSkyline();

            double[] popSizes = new double[data.popSizeCount];
            double[] groupSizes = new double[data.groupSizeCount];


            NewickExporter exporter = new NewickExporter(writer);
            int count = 0;
            try {
                while (importer.hasTree()) {
                    RootedTree inTree = (RootedTree)importer.importNextTree();

                    if (sim == null) {
                        setSamplingTimes(inTree);
                        sim = new TreeSimulator(inTree.getTaxa(), "date");
                    }

                    data.getNextSkyline(popSizes, groupSizes);

                    double[] times = getTimes(inTree, groupSizes);
                    demo = new EmpiricalDemographicFunction(popSizes, times, true);
                    intervals = new CoalescentIntervalGenerator(demo);

                    RootedTree outTree = sim.simulate(intervals);

                    exporter.exportTree(outTree);
                    writer.println();
                    writer.flush();

                    if (count > 0 && count % stepSize == 0) {
                        System.out.print("*");
                        System.out.flush();
                    }
                    count++;
                }
            } catch (ImportException e) {
                e.printStackTrace();
            }

        fileReader.close();


        writer.close();

    }

    private static void setSamplingTimes(final RootedTree tree) {
        setSamplingTimes(tree, tree.getRootNode(), 0.0);

        // the dates are distances from root so make them heights
        double maxLength = 0.0;
        for (Taxon taxon : tree.getTaxa()) {
            double length = (Double)taxon.getAttribute("length");
            if (length > maxLength) {
                maxLength = length;
            }
        }
        for (Taxon taxon : tree.getTaxa()) {
            double length = (Double)taxon.getAttribute("length");
            taxon.setAttribute("date", maxLength - length);
        }
    }

    private static void setSamplingTimes(final RootedTree tree, final Node node, double length) {
        if (tree.isExternal(node)) {
            tree.getTaxon(node).setAttribute("length", length);
        } else {
            for (Node child : tree.getChildren(node)) {
                double l = tree.getLength(child);
                setSamplingTimes(tree, child, length + l);
            }
        }
    }

    private static double[] getTimes(final RootedTree tree, final double[] groupSizes) {
        double[] heights = new double[tree.getInternalNodes().size()];

        int i = 0;
        for (Node node : tree.getInternalNodes()) {
            heights[i] = tree.getHeight(node);
            i++;
        }
        Arrays.sort(heights);

        if (groupSizes != null) {
            double[] allHeights = heights;
            heights = new double[groupSizes.length];
            int k = 0;
            for (int j = 0; j < groupSizes.length; j++) {
                k += groupSizes[j];
                heights[j] = allHeights[k - 1];
            }
        }

        return heights;
    }


    // Main entry point
    static public void main(String[] args) {
        if (args.length > 1) {

            if (args.length != 3) {
                System.err.println("Usage: coalgen <input_file> <tree_file> <output_file>");
                return;
            }

            String inputFileName = args[0];
            String treeFileName = args[1];
            String outputFileName = args[2];

            try {
                CoalGenApp.simulate(inputFileName, treeFileName, outputFileName);
            } catch (IOException e) {
                e.printStackTrace();
            } catch (TraceException e) {
                e.printStackTrace();
            } catch (Importer.ImportException e) {
                e.printStackTrace();
            }

        } else {
            System.setProperty("com.apple.macos.useScreenMenuBar", "true");
            System.setProperty("apple.laf.useScreenMenuBar", "true");

            new CoalGenApp();
        }
    }
}