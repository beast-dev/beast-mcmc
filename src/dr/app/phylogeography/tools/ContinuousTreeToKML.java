/*
 * ContinuousTreeToKML.java
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

package dr.app.phylogeography.tools;

import dr.app.util.Arguments;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.io.TreeImporter;
import jebl.evolution.io.ImportException;
import jebl.evolution.trees.RootedTree;

import java.io.*;

/**
 * @author Philippe Lemey
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 */
public class ContinuousTreeToKML {
    public static final String HELP = "help";
    public static final String ANNOTATION = "annotation";
    public static final String ALTITUDE = "altitude";
    public static final String MRSD = "mrsd";
    public static final String SLICES = "slices";
    public static final String SLICEBW = "slicebw";
    public static final String SLICEMIDPOINT = "slicemidpoint";

    private static final String commandName = "continuous_tree_to_kml";
    private static final PrintStream progressStream = System.out;

    public static final String[] falseTrue = new String[] {"false","true"};

    public static void printUsage(Arguments arguments) {

        arguments.printUsage(commandName, "[<inputTree-file-name>] [<output-file-name>]"); //TODO: set this right
        progressStream.println();
        progressStream.println("  Example: " + commandName + " input.tre output.kml");
        progressStream.println();
    }

    public static void main(String[] args) {

        double altitude = 0;            // altitutude of the root of the 3D trees
        double mostRecentDate = 2010;  // required to convert heights to calendar dates
        String coordinateLabel = "loc";

        boolean makeTreeSlices = false;
        double[] sliceTimes = null;
        double treeSliceBranchWidth = 3;
        boolean showBranchAtMidPoint = false; // shows complete branch for slice if time is more recent than the branch's midpoint

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.StringOption(ANNOTATION, "location annotation label", "specifies the label used for location coordinates annotation [default=location]"),
                        new Arguments.RealOption(ALTITUDE,"specifies the altitude of the root of the 3D tree [default=no 3D tree]"),
                        new Arguments.RealOption(MRSD,"specifies the most recent sampling data in fractional years to rescale time [default=2010]"),

                        new Arguments.StringOption(SLICES,"time","specifies a slice time-list [default=none]"),
                        new Arguments.StringOption(SLICEMIDPOINT, falseTrue, false,
                                "shows complete branch for sliced tree if time is more recent than the branch's midpoint [default=false"),

                        new Arguments.Option(HELP, "option to print this message")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            progressStream.println(ae);
            printUsage(arguments);
            System.exit(1);
        }

        if (args.length == 0 || arguments.hasOption(HELP)) {
            printUsage(arguments);
            System.exit(0);
        }


        if (arguments.hasOption(MRSD)) {
            mostRecentDate = arguments.getRealOption(MRSD);
        }

        if (arguments.hasOption(ALTITUDE)) {
            altitude = arguments.getRealOption(ALTITUDE);
        }

        String annotationLabel = arguments.getStringOption(ANNOTATION);
        if (annotationLabel != null){
            coordinateLabel = annotationLabel;
        }

        String sliceString = arguments.getStringOption(SLICES);
        if (sliceString != null) {
            makeTreeSlices = true;
            try{
                sliceTimes = DiscreteTreeToKML.parseVariableLengthDoubleArray(sliceString);
            } catch (Arguments.ArgumentException ae){
                System.err.println("error reading slice heights");
                ae.printStackTrace();
                return;
            }
            makeTreeSlices = true;
        }

        if (arguments.hasOption(SLICEBW)) {
            treeSliceBranchWidth = arguments.getRealOption(SLICEBW);
        }

        String midpointString = arguments.getStringOption(SLICEMIDPOINT);
        if (midpointString != null && midpointString.compareToIgnoreCase("true") == 0) {
            showBranchAtMidPoint = true;
        }

        final String[] args2 = arguments.getLeftoverArguments();

        String inputFileName = null;
        String outputFileName = null;

        switch (args2.length) {
            case 0:
                printUsage(arguments);
                System.exit(1);
            case 1:
                inputFileName = args2[0];
                outputFileName = inputFileName + ".kml";
                break;
            case 2:
                inputFileName = args2[0];
                outputFileName = args2[1];
                break;
            default: {
                System.err.println("Unknown option: " + args2[2]);
                System.err.println();
                printUsage(arguments);
                System.exit(1);
            }
        }

        RootedTree tree = null;

        try {
            TreeImporter importer = new NexusImporter(new FileReader(inputFileName));
            tree = (RootedTree)importer.importNextTree();
        } catch (ImportException e) {
            e.printStackTrace();
            return;
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

        ContinuousKML exporter = new ContinuousKML(tree, inputFileName, altitude, mostRecentDate, coordinateLabel);

        try {
            BufferedWriter out1 = new BufferedWriter(new FileWriter(outputFileName));
            StringBuffer buffer = new StringBuffer();

            //we write the general tree stuff, but when making slices we do not include everything in the buffer compilation
            exporter.writeTreeToKML();

            if (makeTreeSlices) {
                for (int i = 0; i < sliceTimes.length; i++) {
//                    System.out.println(sliceTimes[i]);
                    exporter.writeTreeToKML(sliceTimes[i], treeSliceBranchWidth, showBranchAtMidPoint);
                }
            }

            exporter.compileBuffer(buffer,makeTreeSlices);
            out1.write(buffer.toString());
            out1.close();
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

    }
}
