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
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;
import dr.util.Version;
import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RVector;
import org.rosuda.JRI.Rengine;
import org.virion.jam.console.ConsoleApplication;

import javax.swing.*;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
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
	                     String outputFileName
	) throws IOException {

		this.posteriorLimit = posteriorLimit;

		attributeNames.add("height");
		attributeNames.add("length");

		CladeSystem cladeSystem = new CladeSystem();

		totalTrees = 10000;
		totalTreesUsed = 0;

		System.out.println("Reading trees (bar assumes 10,000 trees)...");
		System.out.println("0              25             50             75            100");
		System.out.println("|--------------|--------------|--------------|--------------|");

		int stepSize = totalTrees / 60;
		if (stepSize < 1) stepSize = 1;

		if (targetOption != USER_TARGET_TREE) {
			cladeSystem = new CladeSystem();
			FileReader fileReader = new FileReader(inputFileName);
			TreeImporter importer = new NexusImporter(fileReader);
			try {
				totalTrees = 0;
				while (importer.hasTree()) {
					Tree tree = importer.importNextTree();

					if (totalTrees >= burnin) {
						cladeSystem.add(tree, false);

						totalTreesUsed += 1;
					}

					if (totalTrees > 0 && totalTrees % stepSize == 0) {
						System.out.print("*");
						System.out.flush();
					}
					totalTrees++;
				}
			} catch (Importer.ImportException e) {
				System.err.println("Error Parsing Input Tree: " + e.getMessage());
				return;
			}
			fileReader.close();
			System.out.println();
			System.out.println();

			cladeSystem.calculateCladeCredibilities(totalTreesUsed);

			System.out.println("Total trees read: " + totalTrees);
			if (burnin > 0) {
				System.out.println("Ignoring first " + burnin + " trees.");
			}

			System.out.println("Total unique clades: " + cladeSystem.getCladeMap().keySet().size());
			System.out.println();
		}

		MutableTree targetTree;

		if (targetOption == USER_TARGET_TREE) {
			if (targetTreeFileName != null) {
				System.out.println("Reading user specified target tree, " + targetTreeFileName);

				NexusImporter importer = new NexusImporter(new FileReader(targetTreeFileName));
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

		System.out.println("Collecting node information...");
		System.out.println("0              25             50             75            100");
		System.out.println("|--------------|--------------|--------------|--------------|");

		stepSize = totalTrees / 60;
		if (stepSize < 1) stepSize = 1;

		FileReader fileReader = new FileReader(inputFileName);
		NexusImporter importer = new NexusImporter(fileReader);
		cladeSystem = new CladeSystem(targetTree);
		totalTreesUsed = 0;
		try {
			boolean firstTree = true;
			int counter = 0;
			while (importer.hasTree()) {
				Tree tree = importer.importNextTree();

				if (counter >= burnin) {
					if (firstTree) {
						setupAttributes(tree);
						firstTree = false;
					}

					cladeSystem.collectAttributes(tree);
					totalTreesUsed += 1;
				}
				if (counter > 0 && counter % stepSize == 0) {
					System.out.print("*");
					System.out.flush();
				}
				counter++;

			}

			cladeSystem.calculateCladeCredibilities(totalTreesUsed);
		} catch (Importer.ImportException e) {
			System.err.println("Error Parsing Input Tree: " + e.getMessage());
			return;
		}
		System.out.println();
		System.out.println();
		fileReader.close();

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
					String name = (String) iter.next();
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
		 */
		public CladeSystem(Tree targetTree) {
			this.targetTree = targetTree;
			add(targetTree, true);
		}

		/**
		 * adds all the clades in the tree
		 */
		public void add(Tree tree, boolean includeTips) {
			if (taxonList == null) {
				taxonList = tree;
			}

			// Recurse over the tree and add all the clades (or increment their
			// frequency if already present). The root clade is added too (for
			// annotation purposes).
			addClades(tree, tree.getRoot(), includeTips);
		}

		public Clade getClade(NodeRef node) {
			return null;
		}

		private BitSet addClades(Tree tree, NodeRef node, boolean includeTips) {

			BitSet bits = new BitSet();

			if (tree.isExternal(node)) {

				int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
				bits.set(index);

				if (includeTips) {
					addClade(bits);
				}

			} else {

				for (int i = 0; i < tree.getChildCount(node); i++) {

					NodeRef node1 = tree.getChild(node, i);

					bits.or(addClades(tree, node1, includeTips));
				}

				addClade(bits);
			}

			return bits;
		}

		private void addClade(BitSet bits) {
			Clade clade = cladeMap.get(bits);
			if (clade == null) {
				clade = new Clade(bits);
				cladeMap.put(bits, clade);
			}
			clade.setCount(clade.getCount() + 1);
		}

		public void collectAttributes(Tree tree) {
			collectAttributes(tree, tree.getRoot());
		}

		private BitSet collectAttributes(Tree tree, NodeRef node) {

			BitSet bits = new BitSet();

			if (tree.isExternal(node)) {

				int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
				if (index < 0) {
					throw new IllegalArgumentException("Taxon, " + tree.getNodeTaxon(node).getId() + ", not found in target tree");
				}
				bits.set(index);

			} else {

				for (int i = 0; i < tree.getChildCount(node); i++) {

					NodeRef node1 = tree.getChild(node, i);

					bits.or(collectAttributes(tree, node1));
				}
			}

			collectAttributesForClade(bits, tree, node);

			return bits;
		}

		private void collectAttributesForClade(BitSet bits, Tree tree, NodeRef node) {
			Clade clade = cladeMap.get(bits);
			if (clade != null) {

				if (clade.attributeValues == null) {
					clade.attributeValues = new ArrayList<Object[]>();
				}

				int i = 0;
				Object[] values = new Object[attributeNames.size()];
				for (String attributeName : attributeNames) {
					Object value;
					if (attributeName.equals("height")) {
						value = tree.getNodeHeight(node);
					} else if (attributeName.equals("length")) {
						value = tree.getBranchLength(node);
					} else {
						value = tree.getNodeAttribute(node, attributeName);
						if (value instanceof String && ((String) value).startsWith("\"")) {
							value = ((String) value).replaceAll("\"", "");
						}
					}

					//if (value == null) {
					//    System.out.println("attribute " + attributeNames[i] + " is null.");
					//}

					values[i] = value;
					i++;
				}
				clade.attributeValues.add(values);

				clade.setCount(clade.getCount() + 1);
			}
		}

		public Map getCladeMap() {
			return cladeMap;
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

			int i = 0;
			for (String attributeName : attributeNames) {

				if (clade.attributeValues != null && clade.attributeValues.size() > 0) {
					double[] values = new double[clade.attributeValues.size()];

					HashMap<String, Integer> hashMap = new HashMap<String, Integer>();

					Object[] v = clade.attributeValues.get(0);
					if (v[i] != null) {

						boolean isHeight = attributeName.equals("height");
						boolean isBoolean = v[i] instanceof Boolean;

						boolean isDiscrete = v[i] instanceof String;

						double minValue = Double.MAX_VALUE;
						double maxValue = -Double.MAX_VALUE;

						boolean isDoubleArray = v[i] instanceof Object[] && ((Object[]) v[i])[0] instanceof Double;
						// todo Handle other types of arrays

						double[][] valuesArray = null;
						double[] minValueArray = null;
						double[] maxValueArray = null;
						int lenArray = 0;

						if (isDoubleArray) {
							lenArray = ((Object[]) v[i]).length;

							valuesArray = new double[lenArray][clade.attributeValues.size()];
							minValueArray = new double[lenArray];
							maxValueArray = new double[lenArray];

							for (int k = 0; k < lenArray; k++) {
								minValueArray[k] = Double.MAX_VALUE;
								maxValueArray[k] = -Double.MAX_VALUE;
							}
						}

						for (int j = 0; j < clade.attributeValues.size(); j++) {
							Object value = clade.attributeValues.get(j)[i];
							if (isDiscrete) {
                                final String s = (String) value;
                                if (hashMap.containsKey(s)) {
                                    hashMap.put(s, hashMap.get(s) + 1);
								} else {
									hashMap.put(s, 1);
								}
							} else if (isBoolean) {
								values[j] = (((Boolean) value) ? 1.0 : 0.0);
							} else if (isDoubleArray) {
								Double[] array = (Double[]) value;
								for (int k = 0; k < lenArray; k++) {
									valuesArray[k][j] = array[k];
									if (valuesArray[k][j] < minValueArray[k]) minValueArray[k] = valuesArray[k][j];
									if (valuesArray[k][j] > maxValueArray[k]) maxValueArray[k] = valuesArray[k][j];
								}
							} else {
								values[j] = ((Number) value).doubleValue();
								if (values[j] < minValue) minValue = values[j];
								if (values[j] > maxValue) maxValue = values[j];
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
							if (!isDiscrete) {
								if (!isDoubleArray)
									annotateMeanAttribute(tree, node, attributeName, values);
								else {
									for (int k = 0; k < lenArray; k++) {
										annotateMeanAttribute(tree, node, attributeName + (k + 1), valuesArray[k]);
									}
								}
							} else {
								annotateModeAttribute(tree, node, attributeName, hashMap);
								annotateFrequencyAttribute(tree, node, attributeName, hashMap);
							}
							if (!isBoolean && minValue < maxValue && !isDiscrete && !isDoubleArray) {
								// Basically, if it is a boolean (0, 1) then we don't need the distribution information
								// Likewise if it doesn't vary.
								annotateMedianAttribute(tree, node, attributeName + "_median", values);
								annotateHPDAttribute(tree, node, attributeName + "_95%_HPD", 0.95, values);
								annotateRangeAttribute(tree, node, attributeName + "_range", values);
							}
							if (isDoubleArray) {
								for (int k = 0; k < lenArray; k++) {
									if (minValueArray[k] < maxValueArray[k]) {
										annotateMedianAttribute(tree, node, attributeName + (k + 1) + "_median", valuesArray[k]);
										annotateRangeAttribute(tree, node, attributeName + (k + 1) + "_range", valuesArray[k]);
										if (!processBivariateAttributes || lenArray != 2)
											annotateHPDAttribute(tree, node, attributeName + (k + 1) + "_95%_HPD", 0.95, valuesArray[k]);
									}
								}
								// 2D contours
								if (processBivariateAttributes && lenArray == 2) {

									boolean variationInFirst = (minValueArray[0] < maxValueArray[0]);
									boolean variationInSecond = (minValueArray[1] < maxValueArray[1]);

									if (variationInFirst && !variationInSecond)
										annotateHPDAttribute(tree, node, attributeName + "1" + "_95%_HPD", 0.95, valuesArray[0]);

									if (variationInSecond && !variationInFirst)
										annotateHPDAttribute(tree, node, attributeName + "2" + "_95%_HPD", 0.95, valuesArray[1]);

									if (variationInFirst && variationInSecond)
										annotate2DHPDAttribute(tree, node, attributeName + "_95%HPD", 0.95, valuesArray);
								}
							}
						}
					}
				}
				i++;
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
				if (thisCount == maxCount) {
                    // I hope this is the intention
                    mode = mode.concat("+" + key);
				} else if (thisCount > maxCount) {
					mode = key;
					maxCount = thisCount;
				}
				totalCount += thisCount;
			}
			double freq = (double) maxCount / (double) totalCount;
			tree.setNodeAttribute(node, label, mode);
			tree.setNodeAttribute(node, label + ".prob", freq);
		}

		private void annotateFrequencyAttribute(MutableTree tree, NodeRef node, String label, HashMap<String, Integer> values) {
			String mode = null;
			int maxCount = 0;
			int totalCount = 0;

			for (String key : values.keySet()) {
				int thisCount = values.get(key);
				if (thisCount == maxCount)
					mode = mode.concat("+" + key);
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

		// todo Move rEngine to outer class; create once.
		Rengine rEngine = null;

		private final String[] rArgs = new String[]{"--no-save"};

//	    private int called = 0;

		private final String[] rBootCommands = new String[]{
				"library(MASS)",
				"makeContour = function(var1, var2, prob=0.95, n=50, h=c(1,1)) {" +
						"post1 = kde2d(var1, var2, n = n, h = h); " +
						"dx = diff(post1$x[1:2]); " +
						"dy = diff(post1$y[1:2]); " +
						"sz = sort(post1$z); " +
						"c1 = cumsum(sz) * dx * dy; " +
						"levels = sapply(prob, function(x) { approx(c1, sz, xout = 1 - x)$y }); " +
						"line = contourLines(post1$x, post1$y, post1$z, level = levels); " +
						"return(line) }",
		};

		private String makeRString(double[] values) {
			StringBuffer sb = new StringBuffer("c(");
			sb.append(values[0]);
			for (int i = 1; i < values.length; i++)
				sb.append(",").append(values[i]);
			sb.append(")");
			return sb.toString();
		}

		public static final String CORDINATE = "cordinates";

		private String formattedLocation(double x, double y) {
			return String.format("%5.2f,%5.2f,0.00", x, y);
		}

		private void annotate2DHPDAttribute(MutableTree tree, NodeRef node, String label, double hpd, double[][] values) {

			// Uses R-Java interface, and the HPD routines from 'emdbook' and 'coda'

			if (rEngine == null) {

				if (!Rengine.versionCheck()) {
					throw new RuntimeException("JRI library version mismatch");
				}

				rEngine = new Rengine(rArgs, false, null);

				if (!rEngine.waitForR()) {
					throw new RuntimeException("Cannot load R");
				}

				for (String command : rBootCommands) {
					rEngine.eval(command);
				}

//			    called++;
//			    System.err.println("Init call #"+called);
			}

//			int N = values[0].length / 10;
//		    if (N > 50)
//		        N = 50;
			int N = 25;

			double prob = 0.95;

//		    System.err.println("V0 = "+makeRString(values[0]));
//		    System.err.println("V1 = "+makeRString(values[1]));
//		    System.err.println("N = "+N);

			REXP x = rEngine.eval("makeContour(" +
					makeRString(values[0]) + "," +
					makeRString(values[1]) + "," +
					prob + "," +
					N + ")");

			RVector contourList = x.asVector();
			int numberContours = contourList.size();
			StringBuffer output = new StringBuffer();
			for (int i = 0; i < numberContours; i++) {
				output.append("\n<" + CORDINATE + ">\n");
				RVector oneContour = contourList.at(i).asVector();
				double[] xList = oneContour.at(1).asDoubleArray();
				double[] yList = oneContour.at(2).asDoubleArray();
				for (int k = 0; k < xList.length; k++) {
					output.append(formattedLocation(xList[k], yList[k]) + "\n");
				}
				output.append(formattedLocation(xList[0], yList[0]) + "\n</" + CORDINATE + ">\n");
			}

			tree.setNodeAttribute(node, label, output.toString());

//          // Uses some crap that Marc was trying to write
//		    int gridPoints = 10;
//		    DensityMap densityMap = new DensityMap(gridPoints,gridPoints,values[0],values[1]);
//		    densityMap.writeAsTIFF("Node"+node.getNumber()+"-density.tiff");
//
//		    double criticalValue = 1.0 / values.length;     // Throw away all points with mass less than epsilon
//
//		    ContourGenerator contour = new ContourGenerator(densityMap.getXMidPoints(),
//				                                            densityMap.getYMidPoints(),
//				                                            densityMap.getNormalizedDensity(1.0),
//				                                            new ContourAttrib[]{ new ContourAttrib(criticalValue)});
//		    ContourPath[] paths = null;
//		    try {
//			    paths = contour.getContours();
//		    } catch (InterruptedException e) {
//			    e.printStackTrace();
//		    }
//
//		    System.err.println("Found "+paths.length+" paths.");
//		    for(ContourPath p : paths) {
//			    StringBuffer sb = new StringBuffer("Path =");
//			    double[] pathX = p.getAllX();
//			    double[] pathY = p.getAllY();
//			    for(int i=0; i<pathX.length; i++)
//			        sb.append(" {"+pathX[i]+","+pathY[i]+"}");
//			    sb.append((p.isClosed() ? " closed" : " open"));
//			    System.err.println(sb.toString());
//		    }

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
			List<Object[]> attributeValues = null;
		}

		//
		// Private stuff
		//
		TaxonList taxonList = null;
		Map<BitSet, Clade> cladeMap = new HashMap<BitSet, Clade>();

		Tree targetTree;
	}

	int totalTrees = 0;
	int totalTreesUsed = 0;
	double posteriorLimit = 0.0;

	Set<String> attributeNames = new HashSet<String>();
	TaxonList taxa = null;

	static boolean processBivariateAttributes = false;

	static {
		try {
			System.loadLibrary("jri");
			processBivariateAttributes = true;
			System.err.println("JRI loaded. Will process bivariate attributes");
		} catch (UnsatisfiedLinkError e) {

			System.err.println("JRI not available. Will not process bivariate attributes");
		}
	}

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

			/*ConsoleApplication consoleApp = */
			new ConsoleApplication(nameString, aboutString, icon, true);

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
						outputFileName);

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

		new TreeAnnotator(burnin,
				heights,
				posteriorLimit,
				target,
				targetTreeFileName,
				inputFileName,
				outputFileName);

		System.exit(0);
	}

}

