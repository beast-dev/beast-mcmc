package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.evolution.io.*;
import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.util.Version;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id:$
 */
public class LogCleaner {

	private final static Version version = new BeastVersion();

	public LogCleaner(String logFileName, String treeFileName, String fileStem) throws IOException {

		System.out.println("Cleaning files: '" + logFileName + " & " + treeFileName);
		System.out.println();

		PrintWriter logWriter = new PrintWriter(new FileOutputStream(fileStem + ".log"));
		PrintWriter treeWriter = new PrintWriter(new FileOutputStream(fileStem + ".trees"));

		boolean firstFile = true;
		boolean firstTree = true;
		int stateCount = 0;
		int stateStep = -1;

		String[] titles = null;

		System.out.println();
		File logFile = new File(logFileName);

		if (!logFile.exists()) {
			System.err.println(logFileName + " does not exist!");
			return;
		}

		File treeFile = new File(treeFileName);
		if (!treeFile.exists()) {
			System.err.println(treeFileName + " does not exist!");
			return;
		}


		TreeImporter importer = new NexusImporter(new FileReader(treeFile));
		BufferedReader reader = new BufferedReader(new FileReader(logFile));

		String line = reader.readLine();
		titles = line.split("\t");
		logWriter.println(line);

		boolean done = false;
		boolean holdLog = false;
		boolean holdTree = false;

		double[] values = null;
		
		while (!done) {

			boolean skip = false;
			int logState = 0;

			if (!holdLog) {
				line = reader.readLine();

				values = new double[titles.length - 1];

				if (line == null) {
					done = true;
				} else {

					String[] cols = line.split("\t");
					if (cols.length != titles.length) {
						skip = true;
					}

					try {
						logState = Integer.parseInt(cols[0]);
						if (stateStep < 0 && logState > 0) {
							stateStep = logState;
						}

						if (stateStep > 0) {
							stateCount += stateStep;
						}

						for (int i = 1; i < cols.length; i++) {
							values[i - 1] = Double.parseDouble(cols[i]);
						}

					} catch(NumberFormatException nfe) {
						skip = true;
					}
				}
			}

			Tree tree = null;

			if (holdTree) {
				try {
					if (!importer.hasTree()) {
						done = true;
					}
					tree = importer.importNextTree();
				} catch (Importer.ImportException e) {
					skip = true;
				}
			}

			if (!skip) {
				if (firstTree) {
					startLog(tree, treeWriter);
					firstTree = false;
				}

				String name = tree.getId();
				// split on underscore in STATE_xxxx
				String[] bits = name.split("_");
				int treeState = Integer.parseInt(bits[1]);
				if (logState != treeState) {
					System.out.println("Log & Trees file are out of alignment: " + logState + " & " + treeState);
				}

				writeLog(stateCount, values, logWriter);
				writeTree(stateCount, tree, treeWriter);
			}
		}

		stopLog(treeWriter);
		treeWriter.close();
		logWriter.close();
	}

	private Map taxonMap = new HashMap();

	private void startLog(Tree tree, PrintWriter writer) {

		int taxonCount = tree.getTaxonCount();
		writer.println("#NEXUS");
		writer.println("");
		writer.println("Begin taxa;");
		writer.println("\tDimensions ntax=" + taxonCount + ";");
		writer.println("\tTaxlabels");
		for (int i = 0; i < taxonCount; i++) {
			writer.println("\t\t" + tree.getTaxon(i).getId());
		}
		writer.println("\t\t;");
		writer.println("End;");
		writer.println("");
		writer.println("Begin trees;");

// This is needed if the trees use numerical taxon labels
		writer.println("\tTranslate");
		for (int i = 0; i < taxonCount; i++) {
			int k = i + 1;
			Taxon taxon = tree.getTaxon(i);
			taxonMap.put(taxon.getId(), new Integer(k));
			if (k < taxonCount) {
				writer.println("\t\t" + k + " " + taxon.getId() + ",");
			} else {
				writer.println("\t\t" + k + " " + taxon.getId());
			}
		}
		writer.println("\t\t;");
	}

	private void writeLog(int state, double[] values, PrintWriter writer) {
		writer.print(state);
		for (int i = 0; i < values.length; i++) {
			writer.print("\t");
			writer.print(values[i]);
		}
		writer.println();
	}

	private void writeTree(int state, Tree tree, PrintWriter writer) {

		StringBuffer buffer = new StringBuffer("tree STATE_");
		buffer.append(state);
		Double lnP = (Double)tree.getAttribute("lnP");
		if (lnP != null) {
			buffer.append(" [&lnP=").append(lnP).append("]");
		}

		buffer.append(" = [&R] ");

		writeTree(tree, tree.getRoot(), taxonMap, buffer);

		buffer.append(";");
		writer.println(buffer.toString());
	}

	private void writeTree(Tree tree, NodeRef node, Map taxonMap, StringBuffer buffer) {

		NodeRef parent = tree.getParent(node);

		if (tree.isExternal(node)) {
			String taxon = tree.getNodeTaxon(node).getId();
			Integer taxonNo = (Integer)taxonMap.get(taxon);
			if (taxonNo == null) {
				throw new IllegalArgumentException("Taxon, " + taxon + ", not recognized from first tree file");
			}
			buffer.append(taxonNo);
		} else {
			buffer.append("(");
			writeTree(tree, tree.getChild(node, 0), taxonMap, buffer);
			for (int i = 1; i < tree.getChildCount(node); i++) {
				buffer.append(",");
				writeTree(tree, tree.getChild(node, i), taxonMap, buffer);
			}
			buffer.append(")");
		}

		boolean hasAttribute = false;
		Iterator iter = tree.getNodeAttributeNames(node);
		while (iter != null && iter.hasNext()) {
			String name = (String)iter.next();
			Object value = tree.getNodeAttribute(node, name);

			if (!hasAttribute) {
				buffer.append(":[&");
				hasAttribute = true;
			} else {
				buffer.append(",");
			}
			buffer.append(name).append("=").append(value);
		}

		if (hasAttribute) {
			buffer.append("]");
		}

		if (parent != null) {
			if (!hasAttribute) {
				buffer.append(":");
			}
			double length = tree.getBranchLength(node);
			buffer.append(scientificFormatter.format(length));
		}
	}

	private void stopLog(PrintWriter writer) {
		writer.println("End;");
	}

	private static final DecimalFormat decimalFormatter = new DecimalFormat("#.############");
	private static final DecimalFormat scientificFormatter = new DecimalFormat("#.############E0");

	public static void printTitle() {
		System.out.println();
		centreLine("LogCleaner v1.4.2, 2006-2007", 60);
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

		System.out.println("Usage: logcleaner <log-file-name> <tree-file-name> <output-file-stem>");
		System.out.println();
		System.out.println("  Example: logcleaner test1.log test1.tree cleaned");
		System.out.println();

	}

	//Main method
	public static void main(String[] args) throws IOException {

		printTitle();

		if (args.length != 3) {
			printUsage();
			System.exit(1);
		}


		new LogCleaner(args[0], args[1], args[2]);

		System.out.println("Finished.");

		System.exit(0);
	}
}
