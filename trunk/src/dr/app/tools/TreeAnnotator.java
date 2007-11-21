/*
 * TreeAnnotator.java
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
import dr.evolution.io.*;
import dr.evolution.tree.*;
import dr.evolution.util.TaxonList;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;
import dr.util.Version;
import org.virion.jam.console.ConsoleApplication;

import javax.swing.*;
import java.io.*;
import java.util.*;

public class TreeAnnotator {

	private final static Version version = new BeastVersion();

	public final static int MAX_CLADE_CREDIBILITY = 0;
	public final static int MAX_SUM_CLADE_CREDIBILITY = 1;
	public final static int USER_TARGET_TREE = 2;

	public final static int KEEP_HEIGHTS = 0;
	public final static int MEAN_HEIGHTS = 1;
	public final static int MEDIAN_HEIGHTS = 2;

	public TreeAnnotator(int burnin,
	                     int heightsOption,
	                     double posteriorLimit,
	                     int targetOption,
	                     String targetTreeFileName,
	                     String inputFileName,
	                     String outputFileName,
	                     String densityAttributeName,
	                     int timeBinCount,
	                     int valueBinCount,
	                     double timeLimit,
	                     String densityFileName
	) throws IOException {

		this.posteriorLimit = posteriorLimit;

		attributeNames.add("height");
		attributeNames.add("length");

		System.out.println("Reading trees...");

		DensityMap densityMap = null;
		if (densityAttributeName != null) {
			densityMap = new DensityMap(timeBinCount, valueBinCount, timeLimit,0.01);
		}

		CladeSystem cladeSystem = null;

		if (outputFileName != null) {
			cladeSystem = new CladeSystem();
		}

		boolean firstTree = true;
		FileReader fileReader = new FileReader(inputFileName);
		TreeImporter importer = new NexusImporter(fileReader);
		try {
			while (importer.hasTree()) {
				Tree tree = importer.importNextTree();

				if (firstTree) {
					setupAttributes(tree);
					firstTree = false;
				}

				if (totalTrees >= burnin) {
					if (cladeSystem != null) {
						cladeSystem.add(tree);
					}

					if (densityMap != null) {
						densityMap.calibrate(tree, densityAttributeName);
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

		if (densityMap != null) {
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
						densityMap.addTree(tree, densityAttributeName);
					}
					totalTrees += 1;

				}
			} catch (Importer.ImportException e) {
				System.err.println("Error Parsing Input Tree: " + e.getMessage());
				return;
			}
		}

		if (densityMap != null && densityFileName != null) {
			PrintWriter printWriter = new PrintWriter(densityFileName);
			printWriter.println(densityMap.toString());
			printWriter.close();
		}

		if (outputFileName == null) {
			return;
		}

		cladeSystem.calculateCladeCredibilities(totalTreesUsed);

		System.out.println("Total trees read: " + totalTrees);
		if (burnin > 0) {
			System.out.println("Ignoring first " + burnin + " trees.");
		}

		MutableTree targetTree;

		if (targetOption == USER_TARGET_TREE) {
			if (targetTreeFileName != null) {
				System.out.println("Reading user specified target tree, " + targetTreeFileName);

				importer = new NexusImporter(new FileReader(targetTreeFileName));
				try {
					targetTree = new FlexibleTree(importer.importNextTree());
				} catch (Importer.ImportException e) {
					System.err.println("Error Parsing Target Tree: " + e.getMessage());
					return;
				}
			} else {
				System.err.println("No user target tree specified.");
				return;
			}
		} else if (targetOption == MAX_CLADE_CREDIBILITY) {
			System.out.println("Finding maximum credibility tree...");
			targetTree = new FlexibleTree(summarizeTrees(burnin, cladeSystem, inputFileName, false));
		} else if (targetOption == MAX_SUM_CLADE_CREDIBILITY) {
			System.out.println("Finding maximum sum clade credibility tree...");
			targetTree = new FlexibleTree(summarizeTrees(burnin, cladeSystem, inputFileName, true));
		} else {
			throw new RuntimeException("Unknown target tree option");
		}

		System.out.println("Annotating target tree...");
		cladeSystem.annotateTree(targetTree, targetTree.getRoot(), null, heightsOption);

		System.out.println("Writing annotated tree....");
		if (outputFileName != null) {
			NexusExporter exporter = new NexusExporter(new PrintStream(new FileOutputStream(outputFileName)));
			exporter.exportTree(targetTree);
		} else {
			NexusExporter exporter = new NexusExporter(System.out);
			exporter.exportTree(targetTree);
		}

	}

	private void setupAttributes(Tree tree) {
		for (int i = 0; i < tree.getNodeCount(); i++) {
			NodeRef node = tree.getNode(i);
			Iterator iter = tree.getNodeAttributeNames(node);
			if (iter != null) {
				while (iter.hasNext()) {
					String name = (String)iter.next();
					attributeNames.add(name);
				}
			}
		}
	}

	private Tree summarizeTrees(int burnin, CladeSystem cladeSystem, String inputFileName, boolean useSumCladeCredibility) throws IOException {

		Tree bestTree = null;
		double bestScore = Double.NEGATIVE_INFINITY;

		System.out.println("Analyzing " + totalTreesUsed + " trees...");
		System.out.println("0              25             50             75            100");
		System.out.println("|--------------|--------------|--------------|--------------|");

		int stepSize = totalTrees / 60;
		if (stepSize < 1) stepSize = 1;

		int counter = 0;
		TreeImporter importer = new NexusImporter(new FileReader(inputFileName));
		try {
			while (importer.hasTree()) {
				Tree tree = importer.importNextTree();

				if (counter >= burnin) {
					double score = scoreTree(tree, cladeSystem, useSumCladeCredibility);
//                    System.out.println(score);
					if (score > bestScore) {
						bestTree = tree;
						bestScore = score;
					}
				}
				if (counter > 0 && counter % stepSize == 0) {
					System.out.print("*");
					System.out.flush();
				}
				counter++;
			}
		} catch (Importer.ImportException e) {
			System.err.println("Error Parsing Input Tree: " + e.getMessage());
			return null;
		}
		System.out.println();
		System.out.println();
		if (useSumCladeCredibility) {
			System.out.println("Highest Sum Clade Credibility: " + bestScore);
		} else {
			System.out.println("Highest Log Clade Credibility: " + bestScore);
		}

		return bestTree;
	}

	private double scoreTree(Tree tree, CladeSystem cladeSystem, boolean useSumCladeCredibility) {
		if (useSumCladeCredibility) {
			return cladeSystem.getSumCladeCredibility(tree, tree.getRoot(), null);
		} else {
			return cladeSystem.getLogCladeCredibility(tree, tree.getRoot(), null);
		}
	}

	private class CladeSystem {
		//
		// Public stuff
		//

		/**
		 */
		public CladeSystem() {
		}

		/**
		 * adds all the clades in the tree
		 */
		public void add(Tree tree) {
			if (taxonList == null) {
				taxonList = tree;
			}

			// Recurse over the tree and add all the clades (or increment their
			// frequency if already present). The root clade is added too (for
			// annotation purposes).
			addClades(tree, tree.getRoot(), null);
		}

		public Map getCladeMap() {
			return cladeMap;
		}

		public Clade getClade(NodeRef node) {
			return null;
		}

		private void addClades(Tree tree, NodeRef node, BitSet bits) {

			BitSet bits2 = new BitSet();

			if (tree.isExternal(node)) {

				int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
				bits2.set(index);

			} else {

				for (int i = 0; i < tree.getChildCount(node); i++) {

					NodeRef node1 = tree.getChild(node, i);

					addClades(tree, node1, bits2);
				}
			}

			addClade(bits2, tree, node);

			if (bits != null) {
				bits.or(bits2);
			}
		}

		private void addClade(BitSet bits, Tree tree, NodeRef node) {
			Clade clade = cladeMap.get(bits);
			if (clade == null) {
				clade = new Clade(bits);
				cladeMap.put(bits, clade);
			}
			clade.setCount(clade.getCount() + 1);

			if (attributeNames != null) {
				if (clade.attributeLists == null) {
					clade.attributeLists = new List[attributeNames.size()];
					for (int i = 0; i < attributeNames.size(); i++) {
						clade.attributeLists[i] = new ArrayList();
					}
				}

				for (int i = 0; i < attributeNames.size(); i++) {
					String attributeName = attributeNames.get(i);
					Object value;
					if (attributeName.equals("height")) {
						value = tree.getNodeHeight(node);
					} else if (attributeName.equals("length")) {
						value = tree.getBranchLength(node);
					} else {
						value = tree.getNodeAttribute(node, attributeName);
					}

					//if (value == null) {
					//    System.out.println("attribute " + attributeNames[i] + " is null.");
					//}

					if (value != null) {
						clade.attributeLists[i].add(value);
					}
				}
			}
		}

		public void calculateCladeCredibilities(int totalTreesUsed) {
			for (Clade clade : cladeMap.values()) {
				clade.setCredibility(((double) clade.getCount()) / totalTreesUsed);
			}
		}

		public double getSumCladeCredibility(Tree tree, NodeRef node, BitSet bits) {

			double sum = 0.0;

			if (tree.isExternal(node)) {

				int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
				bits.set(index);
			} else {

				BitSet bits2 = new BitSet();
				for (int i = 0; i < tree.getChildCount(node); i++) {

					NodeRef node1 = tree.getChild(node, i);

					sum += getSumCladeCredibility(tree, node1, bits2);
				}

				sum += getCladeCredibility(bits2);

				if (bits != null) {
					bits.or(bits2);
				}
			}

			return sum;
		}

		public double getLogCladeCredibility(Tree tree, NodeRef node, BitSet bits) {

			double logCladeCredibility = 0.0;

			if (tree.isExternal(node)) {

				int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
				bits.set(index);
			} else {

				BitSet bits2 = new BitSet();
				for (int i = 0; i < tree.getChildCount(node); i++) {

					NodeRef node1 = tree.getChild(node, i);

					logCladeCredibility += getLogCladeCredibility(tree, node1, bits2);
				}

				logCladeCredibility += Math.log(getCladeCredibility(bits2));

				if (bits != null) {
					bits.or(bits2);
				}
			}

			return logCladeCredibility;
		}

		private double getCladeCredibility(BitSet bits) {
			Clade clade = cladeMap.get(bits);
			if (clade == null) {
				return 0.0;
			}
			return clade.getCredibility();
		}

		public void annotateTree(MutableTree tree, NodeRef node, BitSet bits, int heightsOption) {

			BitSet bits2 = new BitSet();

			if (tree.isExternal(node)) {

				int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
				bits2.set(index);

				annotateNode(tree, node, bits2, true, heightsOption);
			} else {

				for (int i = 0; i < tree.getChildCount(node); i++) {

					NodeRef node1 = tree.getChild(node, i);

					annotateTree(tree, node1, bits2, heightsOption);
				}

				annotateNode(tree, node, bits2, false, heightsOption);
			}

			if (bits != null) {
				bits.or(bits2);
			}
		}

		private void annotateNode(MutableTree tree, NodeRef node, BitSet bits, boolean isTip, int heightsOption) {
			Clade clade = cladeMap.get(bits);
			if (clade == null) {
				throw new RuntimeException("Clade missing");
			}

			boolean filter = false;
			if (!isTip) {
				double posterior = clade.getCredibility();
				tree.setNodeAttribute(node, "posterior", posterior);
				if (posterior < posteriorLimit) {
					filter = true;
				}
			}

			for (int i = 0; i < attributeNames.size(); i++) {
				String attributeName = attributeNames.get(i);


				double[] values = new double[clade.attributeLists[i].size()];
				HashMap<String, Integer> hashMap = new HashMap<String, Integer>();

				if (values.length > 0) {
					Object v = clade.attributeLists[i].get(0);

					boolean isHeight = attributeName.equals("height");
					boolean isBoolean = v instanceof Boolean;

					boolean isDiscrete = v instanceof String;

					double minValue = Double.MAX_VALUE;
					double maxValue = -Double.MAX_VALUE;
					for (int j = 0; j < clade.attributeLists[i].size(); j++) {
						if (!isDiscrete) {
							values[j] = ((Number) clade.attributeLists[i].get(j)).doubleValue();
							if (values[j] < minValue) minValue = values[j];
							if (values[j] > maxValue) maxValue = values[j];
						} else {
							String value = (String) clade.attributeLists[i].get(j);
							if (value.startsWith("\"")) {
								value = value.replaceAll("\"", "");
							}
							if (hashMap.containsKey(value)) {
								int count = hashMap.get(value);
								hashMap.put(value, count + 1);
							} else {
								hashMap.put(value, 1);
							}
						}
					}
					if (isHeight) {
						if (heightsOption == MEAN_HEIGHTS) {
							double mean = DiscreteStatistics.mean(values);
							tree.setNodeHeight(node, mean);
						} else if (heightsOption == MEDIAN_HEIGHTS) {
							double median = DiscreteStatistics.median(values);
							tree.setNodeHeight(node, median);
						} else {
							// keep the existing height
						}
					}

					if (!filter) {
						if (!isDiscrete)
							annotateMeanAttribute(tree, node, attributeName, values);
						else
							annotateModeAttribute(tree, node, attributeName, hashMap);
						if (!isBoolean && minValue < maxValue && !isDiscrete) {
							// Basically, if it is a boolean (0, 1) then we don't need the distribution information
							// Likewise if it doesn't vary.
							annotateMedianAttribute(tree, node, attributeName + "_median", values);
							annotateHPDAttribute(tree, node, attributeName + "_95%_HPD", 0.95, values);
							annotateRangeAttribute(tree, node, attributeName + "_range", values);
						}
					}
				}
			}
		}

		private void annotateMeanAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
			double mean = DiscreteStatistics.mean(values);
			tree.setNodeAttribute(node, label, mean);
		}

		private void annotateMedianAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
			double median = DiscreteStatistics.median(values);
			tree.setNodeAttribute(node, label, median);

		}

		private void annotateModeAttribute(MutableTree tree, NodeRef node, String label, HashMap<String, Integer> values) {
			String mode = null;
			int maxCount = 0;
			int totalCount = 0;

			for (String key : values.keySet()) {
				int thisCount = values.get(key);
				if (thisCount == maxCount)
					mode.concat("+" + key);
				else if (thisCount > maxCount) {
					mode = key;
					maxCount = thisCount;
				}
				totalCount += thisCount;
			}
			double freq = (double) maxCount / (double) totalCount;
			tree.setNodeAttribute(node, label, mode);
			tree.setNodeAttribute(node, label + ".prob", freq);
		}

		private void annotateRangeAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
			double min = DiscreteStatistics.min(values);
			double max = DiscreteStatistics.max(values);
			tree.setNodeAttribute(node, label, new Object[]{min, max});
		}

		private void annotateHPDAttribute(MutableTree tree, NodeRef node, String label, double hpd, double[] values) {
			int[] indices = new int[values.length];
			HeapSort.sort(values, indices);

			double minRange = Double.MAX_VALUE;
			int hpdIndex = 0;

			int diff = (int) Math.round(hpd * (double) values.length);
			for (int i = 0; i <= (values.length - diff); i++) {
				double minValue = values[indices[i]];
				double maxValue = values[indices[i + diff - 1]];
				double range = Math.abs(maxValue - minValue);
				if (range < minRange) {
					minRange = range;
					hpdIndex = i;
				}
			}
			double lower = values[indices[hpdIndex]];
			double upper = values[indices[hpdIndex + diff - 1]];
			tree.setNodeAttribute(node, label, new Object[]{lower, upper});
		}

		class Clade {
			public Clade(BitSet bits) {
				this.bits = bits;
				count = 0;
				credibility = 0.0;
			}

			public int getCount() {
				return count;
			}

			public void setCount(int count) {
				this.count = count;
			}

			public double getCredibility() {
				return credibility;
			}

			public void setCredibility(double credibility) {
				this.credibility = credibility;
			}

			public boolean equals(Object o) {
				if (this == o) return true;
				if (o == null || getClass() != o.getClass()) return false;

				final Clade clade = (Clade) o;

				return !(bits != null ? !bits.equals(clade.bits) : clade.bits != null);

			}

			public int hashCode() {
				return (bits != null ? bits.hashCode() : 0);
			}

			int count;
			double credibility;
			BitSet bits;
			List[] attributeLists = null;
		}

		//
		// Private stuff
		//
		TaxonList taxonList = null;
		Map<BitSet, Clade> cladeMap = new HashMap<BitSet, Clade>();
	}

	int totalTrees = 0;
	int totalTreesUsed = 0;
	double posteriorLimit = 0.0;

	List<String> attributeNames = new ArrayList<String>();
	TaxonList taxa = null;

	public static void printTitle() {
		System.out.println();
		centreLine("TreeAnnotator " + version.getVersionString() + ", " + version.getDateString(), 60);
		centreLine("MCMC Output analysis", 60);
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
		for (int i = 0; i < n1; i++) {
			System.out.print(" ");
		}
		System.out.println(line);
	}


	public static void printUsage(Arguments arguments) {

		arguments.printUsage("treeannotator", "<input-file-name> [<output-file-name>]");
		System.out.println();
		System.out.println("  Example: treeannotator test.trees out.txt");
		System.out.println("  Example: treeannotator -burnin 100 -heights mean test.trees out.txt");
		System.out.println("  Example: treeannotator -burnin 100 -target map.tree test.trees out.txt");
		System.out.println();
	}

	//Main method
	public static void main(String[] args) throws IOException {

		String targetTreeFileName = null;
		String inputFileName = null;
		String outputFileName = null;

		if (args.length == 0) {
			System.setProperty("com.apple.macos.useScreenMenuBar", "true");
			System.setProperty("apple.laf.useScreenMenuBar", "true");
			System.setProperty("apple.awt.showGrowBox", "true");

			java.net.URL url = LogCombiner.class.getResource("/images/utility.png");
			javax.swing.Icon icon = null;

			if (url != null) {
				icon = new javax.swing.ImageIcon(url);
			}

			final String versionString = version.getVersionString();
			String nameString = "TreeAnnotator " + versionString;
			String aboutString = "<html><center><p>" + versionString + ", " + version.getDateString() + "</p>" +
					"<p>by<br>" +
					"Andrew Rambaut and Alexei J. Drummond</p>" +
					"<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
					"<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
					"<p>Department of Computer Science, University of Auckland<br>" +
					"<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +
					"<p>Part of the BEAST package:<br>" +
					"<a href=\"http://beast.bio.ed.ac.uk/\">http://beast.bio.ed.ac.uk/</a></p>" +
					"</center></html>";

			/*ConsoleApplication consoleApp = */ new ConsoleApplication(nameString, aboutString, icon, true);

			printTitle();

			TreeAnnotatorDialog dialog = new TreeAnnotatorDialog(new JFrame());

			if (!dialog.showDialog("TreeAnnotator " + versionString)) {
				return;
			}

			int burnin = dialog.getBurnin();
			double posteriorLimit = dialog.getPosteriorLimit();
			int targetOption = dialog.getTargetOption();
			int heightsOption = dialog.getHeightsOption();

			targetTreeFileName = dialog.getTargetFileName();
			if (targetOption == USER_TARGET_TREE && targetTreeFileName == null) {
				System.err.println("No target file specified");
				return;
			}

			inputFileName = dialog.getInputFileName();
			if (inputFileName == null) {
				System.err.println("No input file specified");
				return;
			}

			outputFileName = dialog.getOutputFileName();
			if (outputFileName == null) {
				System.err.println("No output file specified");
				return;
			}

			try {
				new TreeAnnotator(burnin,
						heightsOption,
						posteriorLimit,
						targetOption,
						targetTreeFileName,
						inputFileName,
						outputFileName,
						null,
						0,
						0,
						0.0,
						null);

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
		}

		printTitle();

		Arguments arguments = new Arguments(
				new Arguments.Option[]{
						//new Arguments.StringOption("target", new String[] { "maxclade", "maxtree" }, false, "an option of 'maxclade' or 'maxtree'"),
						new Arguments.StringOption("heights", new String[]{"keep", "median", "mean"}, false, "an option of 'keep', 'median' or 'mean'"),
						new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in'"),
						new Arguments.RealOption("limit", "the minimum posterior probability for a node to be annoated"),
						new Arguments.StringOption("target", "target_file_name", "specifies a user target tree to be annotated"),
						new Arguments.StringOption("density", "density_attribute", "specifies an attribute to use to create a density map"),
						new Arguments.IntegerOption("time_bins", "the number of bins for the time axis of the density map"),
						new Arguments.IntegerOption("value_bins", "the number of bins for the value axis of the density map"),
						new Arguments.RealOption("time_limit", "the upper time bound for the density map"),
						new Arguments.StringOption("density_file", "density_file_name", "specifies a file name for the density map"),
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

		int heights = KEEP_HEIGHTS;
		if (arguments.hasOption("heights")) {
			String value = arguments.getStringOption("heights");
			if (value.equalsIgnoreCase("mean")) {
				heights = MEAN_HEIGHTS;
			} else if (value.equalsIgnoreCase("median")) {
				heights = MEDIAN_HEIGHTS;
			}
		}

		int burnin = -1;
		if (arguments.hasOption("burnin")) {
			burnin = arguments.getIntegerOption("burnin");
		}

		double posteriorLimit = 0.0;
		if (arguments.hasOption("limit")) {
			posteriorLimit = arguments.getRealOption("limit");
		}

		int target = MAX_CLADE_CREDIBILITY;
		if (arguments.hasOption("target")) {
			target = USER_TARGET_TREE;
			targetTreeFileName = arguments.getStringOption("target");
		}

		String densityAttributeName = null;
		int timeBinCount = 100;
		int valueBinCount = 25;
		double timeLimit = -1.0;
		String densityFileName = "density.map";
		if (arguments.hasOption("density")) {
			densityAttributeName = arguments.getStringOption("density");
		}

		if (arguments.hasOption("time_bins")) {
			timeBinCount = arguments.getIntegerOption("time_bins");
		}

		if (arguments.hasOption("value_bins")) {
			valueBinCount = arguments.getIntegerOption("value_bins");
		}

		if (arguments.hasOption("time_limit")) {
			timeLimit = arguments.getRealOption("time_limit");
		}

		if (arguments.hasOption("density_file")) {
			densityFileName = arguments.getStringOption("density_file");
		}

		String[] args2 = arguments.getLeftoverArguments();

		if (args2.length > 2) {
			System.err.println("Unknown option: " + args2[2]);
			System.err.println();
			printUsage(arguments);
			System.exit(1);
		}

		if (args2.length == 2) {
			targetTreeFileName = null;
			inputFileName = args2[0];
			outputFileName = args2[1];
		} else if (args2.length == 1) {
			targetTreeFileName = null;
			inputFileName = args2[0];
			outputFileName = null;
		} else {
			System.err.println("Missing input file name");
			printUsage(arguments);
			System.exit(1);
		}

		new TreeAnnotator(burnin,
				heights,
				posteriorLimit,
				target,
				targetTreeFileName,
				inputFileName,
				outputFileName,
				densityAttributeName,
				timeBinCount,
				valueBinCount,
				timeLimit,
				densityFileName);

		System.exit(0);
	}

}

