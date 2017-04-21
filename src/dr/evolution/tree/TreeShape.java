/*
 * TreeShape.java
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

package dr.evolution.tree;

import dr.evolution.io.NexusImporter;

import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;

/**
 * A collection of tree shape summary statistics.
 *
 * @version $Id: TreeShape.java,v 1.7 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class TreeShape {

	/**
	 * The NBar statistic is the average number of nodes above a terminal node.
	 * see Kirkpatrick & Slatkin (1992)
	 */
	public static double getNBarStatistic(Tree tree) {
	
		double NBar = 0.0;	
		for (int i =0; i < tree.getExternalNodeCount(); i++) {
			NodeRef node = tree.getExternalNode(i);
			while (!tree.isRoot(node)) {
				node = tree.getParent(node);
				NBar += 1.0;
			}
		} 
		return NBar / tree.getExternalNodeCount();
	}
	
	/**
	 * The varNBar statistic is the average number of nodes above a terminal node.
	 * see Kirkpatrick & Slatkin (1992)
	 */
	public static double getVarNBarStatistic(Tree tree) {
	
		double[] N = new double[tree.getExternalNodeCount()];	
		for (int i =0; i < N.length; i++) {
			
			NodeRef node = tree.getExternalNode(i);
			while (!tree.isRoot(node)) {
				node = tree.getParent(node);
				N[i] += 1.0;
			}
		} 
		return dr.stats.DiscreteStatistics.variance(N);
	}
	
	/**
	 * The C statistic is the normalized sum of differences of
	 * number of children in left and right subtrees over all
	 * internal nodes.
	 * see Kirkpatrick & Slatkin (1992)
	 */
	public static double getCStatistic(Tree tree) {
	
		double C = 0.0;
		int n = tree.getExternalNodeCount();
		for (int i =0; i < n-1; i++) {
			
			NodeRef node = tree.getInternalNode(i);
			
			int r = TreeUtils.getLeafCount(tree, tree.getChild(node, 0));
			int s = TreeUtils.getLeafCount(tree, tree.getChild(node, 1));
			
			C += Math.abs(r-s);
		} 
		C *= 2.0 / (n * (n - 3) + 2);
		return C;
	}
	
	/**
	 * The B1 statistic measures the maximum number of nodes between
	 * each interior node and a tip (Mi) for all internal nodes except the
	 * root. THe statistic is the sum of the reciprocals of this Mi
	 * statistic.
	 * see Kirkpatrick & Slatkin (1992)
	 */
	public static double getB1Statistic(Tree tree) {
	
		double B1 = 0.0;
		int n = tree.getInternalNodeCount();
		for (int i =0; i < n; i++) {
			
			NodeRef node = tree.getInternalNode(i);
			if (!tree.isRoot(node)) {
				B1 += 1.0/getMi(tree, node);
			}
		} 
		return B1;
	}
	
	/**
	 * Assumes strictly bifurcating tree 
	 */
	private static int getMi(Tree tree, NodeRef node) {
		int childCount = tree.getChildCount(node);
		if (childCount == 0) return 0;
		int Mi = 0;
		for (int i =0; i < childCount; i++) {
			int mi = getMi(tree, tree.getChild(node, i));
			if (mi > Mi) Mi = mi;
		}
		Mi += 1;
		return Mi;
	}
	
	public static double getTreeness(Tree tree) {
		double noise = 0.0;
		double signal = 0.0;
		for (int i = 0; i < tree.getExternalNodeCount(); i++) {
			NodeRef node = tree.getExternalNode(i);
			NodeRef parent = tree.getParent(node);
			noise += tree.getNodeHeight(parent) - tree.getNodeHeight(node);
		}
		
		for (int i = 0; i < tree.getInternalNodeCount(); i++) {
			NodeRef node = tree.getInternalNode(i);
			if (!tree.isRoot(node)) {
				NodeRef parent = tree.getParent(node);
				signal += tree.getNodeHeight(parent) - tree.getNodeHeight(node);
			}
		}
		return signal/noise;
	}
	
	/**
	 * @return the gamma statistic of tree shape (see Pybus & Harvey, 2000).
	 * Assumes a strictly bifurcating tree and contemporaneous sequences.
	 */
	public static double getGammaStatistic(Tree tree) {

		int n = tree.getExternalNodeCount();
		double[] g = getIntervals(tree);
		
		double T = 0; // total branch length
		for (int j = 2; j <= n; j++) {
			T += j * g[j-2];
		}
		
		double gamma = 0.0;
		for (int i = 2; i < n; i++) {
			for (int k = 2; k <= i; k++) {
				gamma += k * g[k-2];
			}
		}
		
		gamma *= (1.0 / n-2.0);
		gamma -= T/2.0;
		gamma /= T * Math.sqrt(1.0/(12.0*(n-2.0)));
		
		return gamma;
	}

	/**
	 * @return the intervals in an ultrametric tree in order from root to tips.
	 */
	private static double[] getIntervals(Tree tree) {
	
		ArrayList<Double> heights = new ArrayList<Double>();
		if (TreeUtils.isUltrametric(tree)) {
			for (int i = 0; i < tree.getInternalNodeCount(); i++) {
				heights.add(tree.getNodeHeight(tree.getInternalNode(i)));
			}
			Collections.sort(heights);
//			for (int i = 0; i < heights.size(); i++) {
//				System.out.print(heights.get(i)+" ");
//			}
//			System.out.println();
			
			double[] intervals = new double[heights.size()];
			for (int i = 0; i < intervals.length-1; i++) {
				double height1 = heights.get(i);
				double height2 = heights.get(i + 1);
				
				intervals[i] = height1 - height2;
			}
			intervals[intervals.length - 1] = heights.get(intervals.length - 1);
		
			return intervals;
		} else 
			throw new IllegalArgumentException("Expecting ultrametric tree.");
		
	}

	public static void main(String[] args) throws Exception {
		NexusImporter importer = new NexusImporter(new FileReader(args[0]));
		
		
		Tree[] trees = importer.importTrees(null);
		
		System.out.println("File = " + args[0]);
		double[] treeness = new double[trees.length];
		for (int i = 0; i < treeness.length; i++) {
			treeness[i] =  getTreeness(trees[i]);
		}
		System.out.println("Mean treeness = " + dr.stats.DiscreteStatistics.mean(treeness));
		System.out.println("Lower (95%) treeness = " + dr.stats.DiscreteStatistics.quantile(0.025, treeness));
		System.out.println("Upper (95%) treeness = " + dr.stats.DiscreteStatistics.quantile(0.975, treeness));
		
	}

}
