/*
 * DensityPlotter.java
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

package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.Tree;
import dr.util.Version;

import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

/*
 * @author Marc Suchard and Andrew Rambaut
 */

public class DensityPlotter {

	private final static Version version = new BeastVersion();

	public DensityPlotter(int burnin,
	                      String inputFileName,
	                      String outputFileName,
	                      String trait1AttributeName,
	                      String trait2AttributeName,
	                      int timeBinCount,
	                      int value1BinCount,
	                      int value2BinCount,
	                      double timeUpper,
	                      double timeLower,
	                      double value1Upper,
	                      double value1Lower,
	                      double value2Upper,
	                      double value2Lower,
	                      boolean printHeaders,
	                      boolean outputTIFF,
	                      boolean logScale

	) throws IOException {

		System.out.println("Reading trees...");

		DensityMap[] densityMaps;

		if (trait2AttributeName != null) {
			densityMaps = new DensityMap[timeBinCount];
			for (int i = 0; i < densityMaps.length; i++) {
				densityMaps[i] = new DensityMap(i, value1BinCount, value2BinCount,
						value1Upper, value1Lower, value2Upper, value2Lower, logScale);
			}
		} else {
			densityMaps = new DensityMap[]{
					new DensityMap(0, timeBinCount, value1BinCount,
							timeUpper, timeLower, value1Upper, value1Lower, logScale)

			};
		}

		boolean firstTree = true;
		double maxTreeHeight = 0.0;
		FileReader fileReader = new FileReader(inputFileName);
		TreeImporter importer = new NexusImporter(fileReader);
		try {
			while (importer.hasTree()) {
				Tree tree = importer.importNextTree();

				if (firstTree) {
					firstTree = false;
				}

				if (totalTrees >= burnin) {
					if (trait2AttributeName != null) {
						for (int i = 0; i < densityMaps.length; i++) {
							densityMaps[i].calibrate(tree, trait1AttributeName, trait2AttributeName);
						}
					} else {
						densityMaps[0].calibrate(tree, trait1AttributeName);
					}
					if (tree.getNodeHeight(tree.getRoot()) > maxTreeHeight) {
						maxTreeHeight = tree.getNodeHeight(tree.getRoot());
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

		double startTime = 0.0;
		double endTime = maxTreeHeight;
		if (timeUpper != Double.POSITIVE_INFINITY) {
			endTime = timeUpper;
		}
		if (timeLower != Double.NEGATIVE_INFINITY) {
			startTime = timeLower;
		}
		double deltaTime = (endTime - startTime) / (double) (timeBinCount);

		// If we want a density plot then we have to read the trees
		// again - the first time was to get the range of values,
		// this read actually creates the map.
		fileReader = new FileReader(inputFileName);
		importer = new NexusImporter(fileReader);
		try {
			totalTrees = 0;
			while (importer.hasTree()) {
				Tree tree = importer.importNextTree();

				if (totalTrees >= burnin) {
					if (trait2AttributeName != null) {
						double sampleTime = startTime;
						for (int i = 0; i < densityMaps.length; i++) {
							densityMaps[i].addTree(tree, sampleTime, trait1AttributeName, trait2AttributeName);
							sampleTime += deltaTime;
						}
					} else {
						densityMaps[0].addTree(tree, trait1AttributeName);
					}
				}
				totalTrees += 1;

			}
		} catch (Importer.ImportException e) {
			System.err.println("Error Parsing Input Tree: " + e.getMessage());
			return;
		}

//		PrintWriter printWriter = null;
		if (trait2AttributeName != null) {
			for (int i = 0; i < densityMaps.length; i++) {
				if (outputTIFF) {
					densityMaps[i].writeAsTIFF(outputFileName + "." + String.format("%03d", i) + ".tif");

				} else {
					PrintWriter printWriter = new PrintWriter(outputFileName + "." + String.format("%03d", i));
					printWriter.println(densityMaps[i].toString(printHeaders));
					printWriter.close();
				}

			}
		} else {
			PrintWriter printWriter = new PrintWriter(outputFileName);
			printWriter.println(densityMaps[0].toString(printHeaders));
			printWriter.close();
		}
//		printWriter.close();


	}

	int totalTrees = 0;
	int totalTreesUsed = 0;

	public static void printTitle() {
		System.out.println();
		centreLine("DensityPlotter " + version.getVersionString() + ", " + version.getDateString(), 60);
		centreLine("BEAST time vs. parameter density analysis", 60);
		centreLine("by", 60);
		centreLine("Andrew Rambaut, Marc A. Suchard and Alexei J. Drummond", 60);
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

		arguments.printUsage("densityplotter", "<input-file-name> [<output-file-name>]");
		System.out.println();
		System.out.println("  Example: densityplotter -burnin 100 -trait rate test.trees density.plot");
		System.out.println();
	}

	//Main method
	public static void main(String[] args) throws IOException {

		String inputFileName = null;
		String outputFileName = null;

		printTitle();

		Arguments arguments = new Arguments(
				new Arguments.Option[]{
						new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in' [default = 0]"),
						new Arguments.StringOption("trait", "trait_name", "specifies an attribute to use to create a density map [default = rate]"),
						new Arguments.StringOption("trait2", "trait_name", "specifies a second attribute to use to create a density map [default = rate]"),
						new Arguments.IntegerOption("time_bins", "the number of bins for the time axis of the density map [default = 100]"),
						new Arguments.IntegerOption("value_bins", "the number of bins for the value axis of the density map [default = 20]"),
						new Arguments.RealOption("time_upper", "the upper time bound for the density map [default = max tree age]"),
						new Arguments.RealOption("time_lower", "the lower time bound for the density map [default = 0]"),
						new Arguments.RealOption("value_upper", "the upper value bound for the density map [default = max value]"),
						new Arguments.RealOption("value_lower", "the lower value bound for the density map [default = min value]"),
						new Arguments.RealOption("value2_upper", "the upper second value bound for the density map [default = max value]"),
						new Arguments.RealOption("value2_lower", "the lower second value bound for the density map [default = min value]"),
						new Arguments.StringOption("headers", "with_headers", "prints row/column labels in output [default = true"),
						new Arguments.Option("tiff", "output in TIFF format"),
						new Arguments.Option("help", "option to print this message"),
						new Arguments.Option("logScale", "transform trait to log scale")
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

		String trait1AttributeName = "rate";
		String trait2AttributeName = null;
		int timeBinCount = 100;
		int valueBinCount = 25;
		int value2BinCount = 25;
		double timeUpper = Double.POSITIVE_INFINITY;
		double timeLower = Double.NEGATIVE_INFINITY;
		double valueUpper = Double.POSITIVE_INFINITY;
		double valueLower = Double.NEGATIVE_INFINITY;
		double value2Upper = Double.POSITIVE_INFINITY;
		double value2Lower = Double.NEGATIVE_INFINITY;
		boolean logScale = false;

		if (arguments.hasOption("logScale"))
			logScale = true;

		if (arguments.hasOption("trait")) {
			trait1AttributeName = arguments.getStringOption("trait");
		}

		if (arguments.hasOption("trait2")) {
			trait2AttributeName = arguments.getStringOption("trait2");
		}

		if (arguments.hasOption("time_bins")) {
			timeBinCount = arguments.getIntegerOption("time_bins");
		}

		if (arguments.hasOption("value_bins")) {
			valueBinCount = arguments.getIntegerOption("value_bins");
		}
		value2BinCount = valueBinCount;

		if (arguments.hasOption("value2_bins")) {
			value2BinCount = arguments.getIntegerOption("value2_bins");
		}

		if (arguments.hasOption("time_upper")) {
			timeUpper = arguments.getRealOption("time_upper");
		}

		if (arguments.hasOption("time_lower")) {
			timeLower = arguments.getRealOption("time_lower");
		}

		if (arguments.hasOption("value_upper")) {
			valueUpper = arguments.getRealOption("value_upper");
		}

		if (arguments.hasOption("value_lower")) {
			valueLower = arguments.getRealOption("value_lower");
		}

		if (arguments.hasOption("value2_upper")) {
			value2Upper = arguments.getRealOption("value2_upper");
		}

		if (arguments.hasOption("value2_lower")) {
			value2Lower = arguments.getRealOption("value2_lower");
		}

		boolean printHeaders = true;
		if (arguments.hasOption("headers")) {
			String text = arguments.getStringOption("headers");
			if (text.toUpperCase().compareTo("FALSE") == 0)
				printHeaders = false;
		}

		boolean outputTIFF = false;
		if (arguments.hasOption("tiff"))
			outputTIFF = true;

		String[] args2 = arguments.getLeftoverArguments();

		if (args2.length > 2) {
			System.err.println("Unknown option: " + args2[2]);
			System.err.println();
			printUsage(arguments);
			System.exit(1);
		}

		if (args2.length == 2) {
			inputFileName = args2[0];
			outputFileName = args2[1];
		} else {
			System.err.println("Missing input or output file name");
			printUsage(arguments);
			System.exit(1);
		}

		new DensityPlotter(burnin,
				inputFileName,
				outputFileName,
				trait1AttributeName,
				trait2AttributeName,
				timeBinCount,
				valueBinCount,
				value2BinCount,
				timeUpper,
				timeLower,
				valueUpper,
				valueLower,
				value2Upper,
				value2Lower,
				printHeaders, outputTIFF, logScale
		);

		System.exit(0);
	}

}