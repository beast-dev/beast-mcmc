/*
 * NexusConverter.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.app.beast.BeastVersion;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.util.Version;

import java.io.EOFException;
import java.io.FileReader;
import java.io.FileWriter;

public class NexusConverter {

	private final static Version version = new BeastVersion();

	public NexusConverter(String inputFileName, String outputFileName, String fileNameStem)
							throws java.io.IOException {

		System.out.println("Converting NEXUS file, '" + inputFileName + "', to BEAST XML file, '" + outputFileName + "'");

		BeastOptions beastOptions = new BeastOptions();

		beastOptions.fileNameStem = fileNameStem;

		try {
			FileReader reader = new FileReader(inputFileName);

			NexusApplicationImporter importer = new NexusApplicationImporter(reader);

			boolean done = false;

			beastOptions.alignment = null;
			beastOptions.tree = null;
			beastOptions.taxonList = null;

			while (!done) {
				try {

					NexusImporter.NexusBlock block = importer.findNextBlock();

					if (block == NexusImporter.TAXA_BLOCK) {

						if (beastOptions.taxonList != null) {
							throw new NexusImporter.MissingBlockException("TAXA block already defined");
						}

						beastOptions.taxonList = importer.parseTaxaBlock();

					} else if (block == NexusImporter.CALIBRATION_BLOCK) {
						if (beastOptions.taxonList == null) {
							throw new NexusImporter.MissingBlockException("TAXA or DATA block must be defined before a CALIBRATION block");
						}

						importer.parseCalibrationBlock(beastOptions.taxonList);

					} else if (block == NexusImporter.CHARACTERS_BLOCK) {

						if (beastOptions.taxonList == null) {
							throw new NexusImporter.MissingBlockException("TAXA block must be defined before a CHARACTERS block");
						}

						if (beastOptions.alignment != null) {
							throw new NexusImporter.MissingBlockException("CHARACTERS or DATA block already defined");
						}

						beastOptions.alignment = (SimpleAlignment)importer.parseCharactersBlock(beastOptions.taxonList);

					} else if (block == NexusImporter.DATA_BLOCK) {

						if (beastOptions.alignment != null) {
							throw new NexusImporter.MissingBlockException("CHARACTERS or DATA block already defined");
						}

						// A data block doesn't need a taxon block before it
						// but if one exists then it will use it.
						beastOptions.alignment = (SimpleAlignment)importer.parseDataBlock(beastOptions.taxonList);
						if (beastOptions.taxonList == null) {
							beastOptions.taxonList = beastOptions.alignment;
						}

					} else if (block == NexusImporter.TREES_BLOCK) {

						if (beastOptions.taxonList == null) {
							throw new NexusImporter.MissingBlockException("TAXA or DATA block must be defined before a TREES block");
						}

						if (beastOptions.tree != null) {
							throw new NexusImporter.MissingBlockException("TREES block already defined");
						}

						Tree[] trees = importer.parseTreesBlock(beastOptions.taxonList);
						if (trees.length > 0) {
							beastOptions.tree = trees[0];
						}

					} else if (block == NexusApplicationImporter.PAUP_BLOCK) {

						importer.parsePAUPBlock(beastOptions);

					} else if (block == NexusApplicationImporter.MRBAYES_BLOCK) {

						importer.parseMrBayesBlock(beastOptions);

					} else if (block == NexusApplicationImporter.RHINO_BLOCK) {

						importer.parseRhinoBlock(beastOptions);

					} else {
						// Ignore the block..
					}

				} catch (EOFException ex) {
					done = true;
				}
			}

			if (beastOptions.alignment == null) {
				throw new NexusImporter.MissingBlockException("DATA or CHARACTERS block is missing");
			}
		} catch (Importer.ImportException ime) {
			System.out.println("Error parsing imported NEXUS file:");
			System.out.println(ime);
		} catch (java.io.IOException ioe) {
			System.out.println("Error reading or writing a file:");
			System.out.println(ioe);
		} catch (Exception ex) {
			System.out.println("Fatal exception (email the authors)");
			ex.printStackTrace(System.out);
		}

		FileWriter fw = new FileWriter(outputFileName);
		beastOptions.generateXML(fw);
		fw.close();
	}

	public static void printTitle() {

		System.out.println("+-----------------------------------------------\\");
		System.out.println("|           NexusConverter v1.0 2003            |\\");

		String versionString = "BEAST Library: " + version.getVersionString();
		System.out.print("|");
		int n = 47 - versionString.length();
		int n1 = n / 2;
		int n2 = n1 + (n % 2);
		for (int i = 0; i < n1; i++) { System.out.print(" "); }
		System.out.print(versionString);
		for (int i = 0; i < n2; i++) { System.out.print(" "); }
		System.out.println("||");

		System.out.println("|       Alexei Drummond and Andrew Rambaut      ||");
		System.out.println("|              University of Oxford             ||");
		System.out.println("|      http://evolve.zoo.ox.ac.uk/beast/        ||");
		System.out.println("\\-----------------------------------------------\\|");
		System.out.println(" \\-----------------------------------------------\\");
		System.out.println();
	}

	public static void printUsage(Arguments arguments) {

		arguments.printUsage("nexusconverter", "<input-file-name> [<output-file-name>]");
		System.out.println();
		System.out.println("  Example: nexusconverter test.nex");
		System.out.println("  Example: nexusconverter test.nex test.xml");
		System.out.println();

	}

	//Main method
	public static void main(String[] args) throws java.io.IOException {

		boolean help = false;

		printTitle();

		Arguments arguments = new Arguments(
			new Arguments.Option[] {
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
			System.exit(1);
		}

		String inputFileName = null;
		String outputFileName = null;

		String[] args2 = arguments.getLeftoverArguments();

		if (args2.length > 2) {
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
			inputFileName = Utils.getLoadFileName("NexusConverter v1.0 - Select NEXUS input file to convert");
		}

		String fileNameStem = Utils.trimExtensions(inputFileName,
									new String[] {"nex", "NEX", "tre", "TRE", "nexus", "NEXUS"});

		if (outputFileName == null) {
			// No ouput file name was given so construct one from the input name...

			outputFileName = fileNameStem + ".xml";
		}

		new NexusConverter(inputFileName, outputFileName, fileNameStem);

		System.out.println("done.");
		System.exit(0);
	}
}

