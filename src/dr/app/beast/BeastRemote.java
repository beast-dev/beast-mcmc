/*
 * BeastRemote.java
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

package dr.app.beast;

import dr.app.util.Arguments;
import dr.inference.parallel.MPIServices;
import dr.math.MathUtils;
import mpi.MPI;

import java.io.File;

/**
 * @author Marc Suchard
 */
public class BeastRemote extends BeastMain {

    public BeastRemote(File inputFile, BeastConsoleApp consoleApp, boolean verbose, boolean parserWarning) {
        super(inputFile, consoleApp, 0, verbose, parserWarning, true, null, false, null, 0);
    }

    public static void terminateSlaves() {
        int size = MPI.COMM_WORLD.Size();
        for (int i = 1; i < size; i++) {
            System.err.println("Sending kill to process " + i);
            MPIServices.requestTermination(i);
        }

    }

    public static void main(String[] oldArgs) throws java.io.IOException {

        // First populate args from MPI.WORLD

        // String[] args = null;
        MPI.Init(oldArgs);
        System.err.println(oldArgs[0]);
        System.err.println(oldArgs[1]);
        System.err.println(oldArgs[2]);

        int rank = MPI.COMM_WORLD.Rank();
        int argLength = oldArgs.length - 3;
        String[] args = new String[argLength];
        for (int i = 0; i < argLength; i++) {

            args[i] = oldArgs[i + 3];
            System.err.println(i + " : " + args[i]);
            if (args[i].contains(".xml")) { // append rank
                args[i] = args[i].replace(".xml", Integer.toString(rank) + ".xml");
                System.err.println("Attempting to load: " + args[i]);
            }
        }

        /*  for (String str : args)
                      System.err.println(str);*/
        //System.exit(1);


        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.Option("verbose", "verbose XML parsing messages"),
                        new Arguments.Option("warnings", "Show warning messages about BEAST XML file"),
                        new Arguments.Option("window", "provide a console window"),
                        new Arguments.Option("working", "change working directory to input file's directory"),
                        new Arguments.LongOption("seed", "specify a random number generator seed"),
                        new Arguments.Option("help", "option to print this message")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println();
            System.out.println(ae.getMessage());
            System.out.println();
            printTitle();
            printUsage(arguments);
            MPI.Finalize();
            System.exit(1);
        }

        if (arguments.hasOption("help")) {
            printTitle();
            printUsage(arguments);
            MPI.Finalize();
            System.exit(0);
        }

        boolean verbose = arguments.hasOption("verbose");
        boolean parserWarning = arguments.hasOption("pwarning"); // if dev, then auto turn on, otherwise default to turn off
//        boolean window = arguments.hasOption("window");
        boolean working = arguments.hasOption("working");

        long seed = MathUtils.getSeed();
        if (arguments.hasOption("seed")) {
            seed = arguments.getLongOption("seed");
            MathUtils.setSeed(seed);
        }

//		if (System.getProperty("dr.app.beast.main.window", "false").toLowerCase().equals("true")) {
//			window = true;
//		}

        BeastConsoleApp consoleApp = null;

/*		if (window) {
			System.setProperty("com.apple.macos.useScreenMenuBar", "true");
			System.setProperty("apple.laf.useScreenMenuBar", "true");
			System.setProperty("apple.awt.showGrowBox", "true");

			javax.swing.Icon icon = IconUtils.getIcon(BeastMain.class, "images/beast.png");

			String nameString = "BEAST " + version.getVersionString();
			String aboutString = "<html><center><p>Bayesian Evolutionary Analysis Sampling Trees<br>" +
					"Version " + version.getVersionString() + ", 2002-2006</p>" +
					"<p>by<br>" +
					"Alexei J. Drummond and Andrew Rambaut</p>" +
					"<p>Department of Computer Science, University of Auckland<br>" +
					"<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +
					"<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
					"<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
					"<p><a href=\"http://evolve.zoo.ox.ac.uk/beast/\">http://evolve.zoo.ox.ac.uk/beast/</a></p>" +
					"<p>Source code distributed under the GNU LGPL:<br>" +
					"<a href=\"http://github.com/beast-dev/beast-mcmc/\">http://github.com/beast-dev/beast-mcmc/</a></p>" +
					"<p>Additional programming by:<br>" +
					"Roald Forsberg, Gerton Lunter, Sidney Markowitz, Oliver Pybus</p>" +
					"<p>Thanks to Korbinian Strimmer for use of his code</p>" +
					"</center></html>";

			consoleApp = new BeastConsoleApp(nameString, aboutString, icon);
		}*/      // Remote can never be interactive

        String inputFileName;

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 1) {
            System.err.println("Unknown option: " + args2[1]);
            System.err.println();
            printTitle();
            printUsage(arguments);
            MPI.Finalize();
            System.exit(1);
        }

        File inputFile = null;

        if (args2.length > 0) {
            inputFileName = args2[0];
            inputFile = new File(inputFileName);
        }

        /*	if (inputFileName == null) {
                        // No input file name was given so throw up a dialog box...
                        inputFile = Utils.getLoadFile("BEAST " + version.getVersionString() + " - Select XML input file");
                    }*/

        if (inputFile != null && working) {
            System.setProperty("user.dir", inputFile.getParent());
        }

        printTitle();

        System.out.println();
        System.out.println("Random number seed: " + seed);
        System.out.println();

        new BeastRemote(inputFile, consoleApp, verbose, parserWarning);
        //System.err.println("Did I get here");
        if (rank == 0)
            terminateSlaves();
        MPI.Finalize();
    }


}
