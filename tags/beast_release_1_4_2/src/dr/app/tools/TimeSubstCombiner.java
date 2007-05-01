package dr.app.tools;

import dr.util.Version;
import dr.app.beast.BeastVersion;
import dr.evolution.io.*;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;

import java.io.*;
import java.util.*;
import java.text.DecimalFormat;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id:$
 */
public class TimeSubstCombiner {

	private final static Version version = new BeastVersion();

	public TimeSubstCombiner(String timeTreeFileName, String substTreeFileName, String outputFileName) throws IOException {

		System.out.println("Combining files: '" + timeTreeFileName + " & " + substTreeFileName);
		System.out.println();

		PrintWriter treeWriter = new PrintWriter(new FileOutputStream(outputFileName));

		boolean firstTree = true;

		System.out.println();
		File timeTreeFile = new File(timeTreeFileName);

		if (!timeTreeFile.exists()) {
			System.err.println(timeTreeFileName + " does not exist!");
			return;
		}

		File substTreeFile = new File(substTreeFileName);
		if (!substTreeFile.exists()) {
			System.err.println(substTreeFileName + " does not exist!");
			return;
		}


		TreeImporter timeTreeImporter = new NexusImporter(new FileReader(timeTreeFile));
		TreeImporter substTreeImporter = new NexusImporter(new FileReader(substTreeFile));

		try {
			while (timeTreeImporter.hasTree()) {

				Tree timeTree = timeTreeImporter.importNextTree();
				Tree substTree = substTreeImporter.importNextTree();

				if (firstTree) {
					startLog(timeTree, treeWriter);
					firstTree = false;
				}

				String name = timeTree.getId();
				// split on underscore in STATE_xxxx
				String[] bits = name.split("_");
				int treeState = Integer.parseInt(bits[1]);

				writeTree(treeState, timeTree, substTree, treeWriter);

			}

			stopLog(treeWriter);
			treeWriter.close();
		} catch (Importer.ImportException e) {
			System.err.println(e.getMessage());
			return;
		}
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

	private void writeTree(int state, Tree timeTree, Tree substTree, PrintWriter writer) {

		StringBuffer buffer = new StringBuffer("tree STATE_");
		buffer.append(state);
		Double lnP = (Double)timeTree.getAttribute("lnP");
		if (lnP != null) {
			buffer.append(" [&lnP=").append(lnP).append("]");
		}

		buffer.append(" = [&R] ");

		writeTree(timeTree, timeTree.getRoot(), substTree, substTree.getRoot(), taxonMap, buffer);

		buffer.append(";");
		writer.println(buffer.toString());
	}

	private void writeTree(Tree timeTree, NodeRef timeNode, Tree substTree, NodeRef substNode, Map taxonMap, StringBuffer buffer) {

		NodeRef parent = timeTree.getParent(timeNode);

		if (timeTree.isExternal(timeNode)) {
			String taxon = timeTree.getNodeTaxon(timeNode).getId();
			Integer taxonNo = (Integer)taxonMap.get(taxon);
			if (taxonNo == null) {
				throw new IllegalArgumentException("Taxon, " + taxon + ", not recognized from first tree file");
			}
			buffer.append(taxonNo);
		} else {
			buffer.append("(");
			writeTree(timeTree, timeTree.getChild(timeNode, 0), substTree, substTree.getChild(substNode, 0), taxonMap, buffer);
			for (int i = 1; i < timeTree.getChildCount(timeNode); i++) {
				buffer.append(",");
				writeTree(timeTree, timeTree.getChild(timeNode, i), substTree, substTree.getChild(substNode, i), taxonMap, buffer);
			}
			buffer.append(")");
		}

		double time = 0;
		double subst = 0;

		if (parent != null) {
			time = timeTree.getBranchLength(timeNode);
			subst = substTree.getBranchLength(substNode);
			buffer.append(":[&rate=");
			buffer.append(scientificFormatter.format(subst/time));

			Iterator iter = timeTree.getNodeAttributeNames(timeNode);
			while (iter != null && iter.hasNext()) {
				String name = (String)iter.next();
				Object value = timeTree.getNodeAttribute(timeNode, name);

				buffer.append(",").append(name).append("=").append(value);
			}

			buffer.append("]");

			double length = timeTree.getBranchLength(timeNode);
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
		centreLine("TimeSubstCombiner v1.4.2, 2006-2007", 60);
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

		System.out.println("Usage: TimeSubstCombiner <time-tree-file-name> <subst-tree-file-name> <output-file-name>");
		System.out.println();

	}

	//Main method
	public static void main(String[] args) throws IOException {

		printTitle();

		if (args.length != 3) {
			printUsage();
			System.exit(1);
		}


		new TimeSubstCombiner(args[0], args[1], args[2]);

		System.out.println("Finished.");

		System.exit(0);
	}
}
