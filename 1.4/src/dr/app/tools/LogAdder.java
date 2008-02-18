package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.util.Version;

import java.io.*;
import java.text.DecimalFormat;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id:$
 */
public class LogAdder {

	private final static Version version = new BeastVersion();

	public LogAdder(String operation, String logFileName, String column1, String column2, String outputFile, String outColumn) throws IOException {

		System.out.println("Processing file: '" + logFileName + " (" + column1 + ", " + column2 +")");
		System.out.println();

		PrintWriter logWriter = new PrintWriter(new FileOutputStream(outputFile));

		String[] titles = null;

		System.out.println();
		File logFile = new File(logFileName);

		if (!logFile.exists()) {
			System.err.println(logFileName + " does not exist!");
			return;
		}

		BufferedReader reader = new BufferedReader(new FileReader(logFile));

		String line = reader.readLine();
		titles = line.split("\t");

		int columnIndex1 = -1;
		int columnIndex2 = -1;

		logWriter.print(titles[0]);

		for (int i = 1; i < titles.length; i++) {
			if (column1.equalsIgnoreCase(titles[i])) {
				columnIndex1 = i;
			}
			if (column2.equalsIgnoreCase(titles[i])) {
				columnIndex2 = i;
			}
			logWriter.print("\t" + titles[i]);
		}
		logWriter.println("\t" + outColumn);

		if (columnIndex1 < 0) {
			System.err.println("Column " + column1 + " is not found");
			return;
		}

		if (columnIndex2 < 0) {
			System.err.println("Column " + column2 + " is not found");
			return;
		}

		boolean done = false;

		double[] values = null;

		while (!done) {

			line = reader.readLine();

			if (line == null) {
				done = true;
			} else {

				String[] cols = line.split("\t");
				if (cols.length != titles.length) {
					System.err.println("Row does not have the right number of columns");
					return;
				}

				double value1 = Double.parseDouble(cols[columnIndex1]);
				double value2 = Double.parseDouble(cols[columnIndex2]);

				logWriter.print(cols[0]);

				for (int i = 1; i < titles.length; i++) {
					logWriter.print("\t" + cols[i]);
				}

				double value3 = 0;

				if (operation.equalsIgnoreCase("-sum")) {
					value3 = value1 + value2;
				} else if (operation.equalsIgnoreCase("-product")) {
					value3 = value1 * value2;
				} else if (operation.equalsIgnoreCase("-difference")) {
					value3 = value1 - value2;
				} else if (operation.equalsIgnoreCase("-ratio")) {
					value3 = value1 / value2;
				} else {
					System.err.println("Unknown operation, " + operation);
					return;
				}

				logWriter.println("\t" + value3);
			}
		}


		logWriter.close();
	}

	private static final DecimalFormat decimalFormatter = new DecimalFormat("#.############");
	private static final DecimalFormat scientificFormatter = new DecimalFormat("#.############E0");

	public static void printTitle() {
		System.out.println();
		centreLine("LogAdder v1.4.2, 2006-2007", 60);
		centreLine("MCMC Output Cleaner", 60);
		centreLine("by", 60);
		centreLine("Andrew Rambaut and Alexei J. Drummond", 60);
		System.out.println();
		centreLine("Institute of Evolutionary Biology", 60);
		centreLine("University of Edinburgh", 60);
		centreLine("a.rambaut@ed.ac.uk", 60);
		System.out.println();
		centreLine("Department of Computer Science", 60);
		centreLine("University of Auckland", 60);
		centreLine("alexei@cs.auckland.ac.nz", 60);
		System.out.println();
		System.out.println();
	}

	public static void centreLine(String line, int pageWidth) {
		int n = pageWidth - line.length();
		int n1 = n / 2;
		for (int i = 0; i < n1; i++) { System.out.print(" "); }
		System.out.println(line);
	}


	public static void printUsage() {

		System.out.println("Usage: logadder -sum|-difference|-product|-ratio <log-file> <column1> <column2> <output-file> <output-column>");
		System.out.println();
		System.out.println();

	}

	//Main method
	public static void main(String[] args) throws IOException {

		printTitle();

		if (args.length != 6) {
			printUsage();
			System.exit(1);
		}


		new LogAdder(args[0], args[1], args[2], args[3], args[4], args[5]);

		System.out.println("Finished.");

		System.exit(0);
	}
}