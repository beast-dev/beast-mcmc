/*
 * LogCombiner.java
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

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.util.Version;
import org.virion.jam.console.ConsoleApplication;

import javax.swing.*;
import java.io.*;
import java.text.DecimalFormat;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class LogCombiner {

	private final static Version version = new BeastVersion();

	public LogCombiner(int[] burnins, int resample, String[] inputFileNames, String outputFileName, boolean treeFiles, boolean convertToDecimal, boolean useScale, double scale) throws IOException {

		PrintStream outputStream = new PrintStream(new FileOutputStream(outputFileName));

		boolean firstFile = true;
		int stateCount = 0;
		int stateStep = -1;

		String[] titles = null;

		System.out.println();
		for (int i = 0; i < inputFileNames.length; i++) {
			File inputFile = new File(inputFileNames[i]);

			if (!inputFile.exists()) {
				System.err.println(inputFileNames[i] + " does not exist!");
				return;
			}  else if (inputFile.isDirectory()) {
				System.err.println(inputFileNames[i] + " is a directory.");
				return;
			}

			int burnin = burnins[0];
			if (burnins.length > i) {
				burnin = burnins[i];
			}

			if (burnin > 0) {
				System.out.print("Combining file: '" + inputFileNames[i] + "' removing burnin: " + burnin);
			} else {
				System.out.print("Combining file: '" + inputFileNames[i] + "' without removing burnin");
			}

			if (resample > 0) {
				System.out.print(", resampling with frequency: " + resample);
			}

            if (useScale) {
                System.out.println(", rescaling by: " + scale);
            } else {
                System.out.println();
            }

            BufferedReader reader = new BufferedReader(new FileReader(inputFile));
			if (treeFiles) {
				int lineCount = 1;

				String line = reader.readLine();
				while (line != null && !line.trim().toUpperCase().startsWith("TREE")) {
					if (firstFile) {
						outputStream.println(line);
					}
					line = reader.readLine();
					lineCount++;
				}

				while (line != null && !line.trim().toUpperCase().startsWith("END;")) {
					String[] parts = line.split("\\s");
					if (parts.length < 5) {
						System.err.println("Parsing error at line " + lineCount + ", file " + inputFileNames[i] + ": tree missing?\r"+line);
						return;
					}
					String stateString = parts[1];
					String tree = parts[4];

					parts = stateString.split("_");

					int state = Integer.parseInt(parts[1]);

					if (stateStep < 0 && state > 0) {
						stateStep = state;
					}

					if (state >= burnin) {
						if (stateStep > 0) {
							stateCount += stateStep;
						}

						if (resample < 0 || stateCount % resample == 0) {
							if (convertToDecimal || useScale) {
								tree = reformatNumbers(tree, convertToDecimal, useScale, scale);
							}

							outputStream.println("tree STATE_" + stateCount +" = [&R] " + tree);
						}
					}

					line = reader.readLine();
					lineCount++;
				}
			} else {
				int lineCount = 1;
				String line = reader.readLine();
				if (firstFile) {
					titles = line.split("\t");
					outputStream.println(line);
				} else {
					String[] newTitles = line.split("\t");
					if (newTitles.length != titles.length) {
						System.err.println("ERROR: The number of columns in file, " + inputFileNames[i]+ ", does not match that of the first file");
						return;
					}
					for (int k = 0; k < newTitles.length; k++) {
						if (!newTitles[k].equals(titles[k])) {
							System.err.println("WARNING: The column heading, " + newTitles[k] + " in file, " + inputFileNames[i]+ ", does not match the first file's heading, " + titles[k]);
						}
					}
				}
				line = reader.readLine();
				lineCount++;

				while (line != null) {
					String[] parts = line.split("\\s");
					int state = Integer.parseInt(parts[0]);

					if (stateStep < 0 && state > 0) {
						stateStep = state;
					}

					if (state >= burnin) {
						if (stateStep > 0) {
							stateCount += stateStep;
						}

						if (resample < 0 || stateCount % resample == 0) {
							outputStream.print(stateCount);
							for (int j = 1; j < parts.length; j++) {
								String value = parts[j];

								if (useScale) {
									if (titles[j].equals("clock.rate") || titles[j].startsWith("skyline.popSize")) {
										value = reformatNumbers(value, convertToDecimal, true, 1.0 / scale);
									} else if (titles[j].equals("treeModel.rootHeight")) {
										value = reformatNumbers(value, convertToDecimal, true, scale);
									}
								} else {
									value = reformatNumbers(value, convertToDecimal, false, 1.0);
								}
								outputStream.print("\t" + value);
							}
							outputStream.println();
						}
					}

					line = reader.readLine();
					lineCount++;
				}
			}

			firstFile = false;
		}

		if (treeFiles) {
			outputStream.println("end;");
		}
        System.out.println("Finished.");
	}

	private static final DecimalFormat decimalFormatter = new DecimalFormat("#.############");
    private static final DecimalFormat scientificFormatter = new DecimalFormat("#.############E0");

	private String reformatNumbers(String line, boolean convertDecimal, boolean useScale, double scale) {
		StringBuffer outLine = new StringBuffer();

		Pattern pattern = Pattern.compile("\\d+\\.\\d+(E[\\-\\d\\.]+)?");
		Matcher matcher = pattern.matcher(line);

		int lastEnd = 0;
		while(matcher.find()) {
			int start = matcher.start();
			String token = matcher.group();
			double value = Double.parseDouble(token);
			if (useScale) {
				value *= scale;
			}
			String outToken = (convertDecimal ? decimalFormatter.format(value) : scientificFormatter.format(value));

			outLine.append(line.substring(lastEnd, start));
			outLine.append(outToken);

			lastEnd = matcher.end();
		}
		outLine.append(line.substring(lastEnd));
		return outLine.toString();
	}

	public static void printTitle() {

		System.out.println("+-----------------------------------------------\\");
		System.out.println("|              LogCombiner v1.4 2006            |\\");
		System.out.println("|              MCMC Output Combiner             ||");

		String versionString = "BEAST Library: " + version.getVersionString();
		System.out.print("|");
		int n = 47 - versionString.length();
		int n1 = n / 2;
		int n2 = n1 + (n % 2);
		for (int i = 0; i < n1; i++) { System.out.print(" "); }
		System.out.print(versionString);
		for (int i = 0; i < n2; i++) { System.out.print(" "); }
		System.out.println("||");

		System.out.println("|       Andrew Rambaut and Alexei Drummond      ||");
		System.out.println("|              University of Oxford             ||");
		System.out.println("|      http://evolve.zoo.ox.ac.uk/beast/        ||");
		System.out.println("\\-----------------------------------------------\\|");
		System.out.println(" \\-----------------------------------------------\\");
		System.out.println();
	}

	public static void printUsage(Arguments arguments) {

		arguments.printUsage("logcombiner", "[-trees] [-decimal] [-burnin <burnin>] [-resample <frequency>] [-scale <scale>] <input-file-name1> [<input-file-name2> ...] <output-file-name>");
		System.out.println();
		System.out.println("  Example: loganalyser test1.log test2.log combined.log");
		System.out.println("  Example: loganalyser -burnin 10000 test1.log test2.log combined.log");
		System.out.println();

	}

	//Main method
	public static void main(String[] args) throws IOException {

		boolean treeFiles = false;
		boolean convertToDecimal = false;

		int burnin = -1;
		int resample = -1;
		double scale = 1.0;
		boolean useScale = false;

		if (args.length == 0) {
			System.setProperty("com.apple.macos.useScreenMenuBar","true");
			System.setProperty("apple.laf.useScreenMenuBar","true");
			System.setProperty("apple.awt.showGrowBox","true");

			java.net.URL url = LogCombiner.class.getResource("/images/utility.png");
			javax.swing.Icon icon = null;

			if (url != null) {
				icon = new javax.swing.ImageIcon(url);
			}

			String nameString = "LogCombiner v1.4";
			String aboutString = "LogCombiner v1.4\n" +
					"©2006 Andrew Rambaut & Alexei Drummond\n" +
					"University of Oxford";

			ConsoleApplication consoleApp = new ConsoleApplication(nameString, aboutString, icon, true);

			printTitle();

			LogCombinerDialog dialog = new LogCombinerDialog(new JFrame());

			if (!dialog.showDialog("LogCombiner v1.4")) {
				return;
			}

			treeFiles = dialog.isTreeFiles();
			convertToDecimal = dialog.convertToDecimal();
			if (dialog.isResampling()) {
				resample = dialog.getResampleFrequency();
			}

			String[] inputFiles = dialog.getFileNames();
			int[] burnins = dialog.getBurnins();

			String outputFileName = dialog.getOutputFileName();

			if (outputFileName == null) {
				System.err.println("No output file specified");
			}

			try {
				new LogCombiner(burnins, resample, inputFiles, outputFileName, treeFiles, convertToDecimal, useScale, scale);

			} catch (Exception ex) {
				System.err.println("Exception: " + ex.getMessage());
			}
            System.out.println("Finished - Quit program to exit.");
			while (true) {
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
			}
		} else {
			printTitle();

			Arguments arguments = new Arguments(
					new Arguments.Option[] {
							new Arguments.Option("trees", "use this option to combine tree log files"),
							new Arguments.Option("decimal", "this option converts numbers from scientific to decimal notation"),
							new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in'"),
							new Arguments.IntegerOption("resample", "resample the log files to this frequency (the original sampling frequency must be a factor of this value)"),
							new Arguments.RealOption("scale", "a scaling factor that will multiply any time units by this value"),
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
				System.exit(0);
			}

			treeFiles = arguments.hasOption("trees");

			convertToDecimal = arguments.hasOption("decimal");

			burnin = -1;
			if (arguments.hasOption("burnin")) {
				burnin = arguments.getIntegerOption("burnin");
			}

			resample = -1;
			if (arguments.hasOption("resample")) {
				resample = arguments.getIntegerOption("resample");
			}

			scale = 1.0;
			useScale = false;
			if (arguments.hasOption("scale")) {
				scale = arguments.getRealOption("scale");
				useScale = true;
			}

			String[] inputFileNames = null;
			String outputFileName = null;

			String[] args2 = arguments.getLeftoverArguments();

			if (args2.length < 2) {
				System.err.println("Requires a minimum of 1 input filename and 1 output filename");
				System.err.println();
				printUsage(arguments);
				System.exit(1);
			}

			inputFileNames = new String[args2.length - 1];
			for (int i = 0; i < inputFileNames.length; i++) {
				inputFileNames[i] = args2[i];
			}
			outputFileName = args2[args2.length - 1];

			new LogCombiner(new int[] { burnin }, resample, inputFileNames, outputFileName, treeFiles, convertToDecimal, useScale, scale);
		}

		System.exit(0);
	}
}

