/*
 * TreeTraceAnalysis.java
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

package dr.evomodel.tree;

import dr.evolution.io.Importer;
import dr.evolution.tree.CladeSet;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.util.FrequencySet;
import dr.util.NumberFormatter;

import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Set;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: TreeTraceAnalysis.java,v 1.20 2005/06/07 16:28:18 alexei Exp $
 */
public class TreeTraceAnalysis {
	
	private TreeTraceAnalysis(TreeTrace trace, int burnIn, boolean verbose) {

        this(new TreeTrace[] {trace}, burnIn, verbose);
	}

    public TreeTraceAnalysis(TreeTrace[] trace, int burnIn, boolean verbose) {

        this.trace = trace;

        int minMaxState = Integer.MAX_VALUE;
        for (int i = 0; i < trace.length; i++) {
            if (trace[i].getMaximumState() < minMaxState) {
                minMaxState = trace[i].getMaximumState();
            }
        }

        if (burnIn < 0 || burnIn >= minMaxState) {
            this.burnin = minMaxState / (10 * trace[0].getStepSize());
            if (verbose) System.out.println("WARNING: Burn-in larger than total number of states - using 10% of smallest trace");
        } else {
            this.burnin = burnIn;
        }

        analyze(verbose);
    }


	/**
	 * Actually analyzes the trace given the burnin
	 */
	public void analyze(boolean verbose) {

		if (verbose) {
            if (trace.length > 1) System.out.println("Combining " + trace.length + " traces.");
        }

        Tree tree = getTree(0);
		
		cladeSet = new CladeSet(tree);
		treeSet = new FrequencySet();
		treeSet.add(Tree.Utils.uniqueNewick(tree,tree.getRoot()));
		
		for (int t = 0; t < trace.length; t++) {
            int treeCount = trace[t].getTreeCount(burnin*trace[t].getStepSize());
            double stepSize = treeCount/60.0;
            int counter = 1;

            if (verbose) {
                System.out.println("Analyzing " + treeCount + " trees...");
                System.out.println("0              25             50             75            100");
                System.out.println("|--------------|--------------|--------------|--------------|");
                System.out.print(  "*");
            }
            for (int i = 1; i < treeCount; i++) {
                tree = (Tree)trace[t].getTree(i, burnin*trace[t].getStepSize());
                cladeSet.add(tree);
                treeSet.add(Tree.Utils.uniqueNewick(tree, tree.getRoot()));
                if (i >= (int)Math.round(counter*stepSize) && counter <= 60) {
                    if (verbose) {
                        System.out.print("*");
                        System.out.flush();
                    }
                    counter += 1;
                }
            }
            if (verbose) { System.out.println("*"); }
        }
	}
	
	/**
	 * Actually analyzes a particular tree using the trace given the burnin
	 */
	public final Tree analyzeTree(String target) {
			
		int n = getTreeCount();
		
		FlexibleTree meanTree = null;
		
		for (int i = 0; i < n; i++) {
			Tree tree = getTree(i);
				
			if (Tree.Utils.uniqueNewick(tree, tree.getRoot()).equals(target)) {
				meanTree = new FlexibleTree(tree);		
				break;
			}				
		}
		if (meanTree == null) throw new RuntimeException("No target tree in trace");
		
		int m = meanTree.getInternalNodeCount();
		for (int j = 0; j < m; j++) {
			double[] heights = new double[n];
			NodeRef node1 = meanTree.getInternalNode(j);
			Set leafSet = Tree.Utils.getDescendantLeaves(meanTree, node1);	
			
			for (int i = 0; i < n; i++) {
				Tree tree = getTree(i);
				
				NodeRef node2 = Tree.Utils.getCommonAncestorNode(tree, leafSet);
				heights[i] =  tree.getNodeHeight(node2);						
			}
			meanTree.setNodeHeight(node1,dr.stats.DiscreteStatistics.mean(heights));
			meanTree.setNodeAttribute(node1, "upper", new Double(dr.stats.DiscreteStatistics.quantile(0.975, heights)));
			meanTree.setNodeAttribute(node1, "lower", new Double(dr.stats.DiscreteStatistics.quantile(0.025, heights)));
			
		}
		
		return meanTree;
	}

    public final int getTreeCount() {

        int treeCount = 0;
        for (int i = 0; i < trace.length; i++) {
            treeCount += trace[i].getTreeCount(burnin*trace[i].getStepSize());
        }
        return treeCount;
    }

    public final Tree getTree(int index) {

        int oldTreeCount = 0;
        int newTreeCount = 0;
        for (int i = 0; i < trace.length; i++) {
            newTreeCount += trace[i].getTreeCount(burnin*trace[i].getStepSize());

            if (index < newTreeCount) {
                return trace[i].getTree(index-oldTreeCount, burnin*trace[i].getStepSize());
            }
            oldTreeCount = newTreeCount;
        }
        throw new RuntimeException("Couldn't find tree " + index);
    }



	public void report() throws IOException {
	
		int fieldWidth = 14;
		NumberFormatter formatter = new NumberFormatter(6);
		formatter.setPadding(true);
		formatter.setFieldWidth(fieldWidth);

        int n = treeSet.size();
        int totalTrees = treeSet.getSumFrequency();

		System.out.println();
		System.out.println("burnIn=" + burnin);
		System.out.println("total trees used =" + totalTrees);
		System.out.println();
										

		System.out.println("95% credible set (" + n + " unique trees, " + totalTrees + " total):");
		System.out.println("Count\tPercent\tTree");
		int credSet = (95 * totalTrees) / 100;
		int sumFreq = 0;
		
		NumberFormatter nf = new NumberFormatter(8);
		
		for (int i = 0; i < n; i++) {
			int freq = treeSet.getFrequency(i);
			double prop = ((double)freq) / totalTrees;
			System.out.print(freq);
			System.out.print("\t" + nf.formatDecimal(prop * 100.0, 2) + "%");
			
			sumFreq += freq;
			double sumProp = ((double)sumFreq) / totalTrees;
			System.out.print("\t" + nf.formatDecimal(sumProp * 100.0, 2) + "%");
			
			String newickTree = (String)treeSet.get(i);
			
			if (freq > 100) {
				// calculate conditional average node heights
				Tree meanTree = analyzeTree(newickTree);
				System.out.println("\t" + Tree.Utils.newick(meanTree));
				/*for (int k = 0; k < meanTree.getInternalNodeCount(); k++) {
					NodeRef node = meanTree.getInternalNode(k);
					System.out.println("node " + k + "\t" +
						meanTree.getNodeHeight(node) + "\t" +
						meanTree.getNodeAttribute(node, "lower") + "\t" +
						meanTree.getNodeAttribute(node, "upper"));
				}*/

			} else {
				System.out.println("\t" + newickTree);
			}
			
			
			if (sumFreq >= credSet) {
                System.out.println();
                System.out.println("95% credible set has " + (i + 1) + " trees.");                
                break;
            }
		}
			
		System.out.println();
		System.out.println("Majority rule clades (" + cladeSet.size() + " unique clades):");
		n = cladeSet.size();
		for (int i = 0; i < n; i++) {
			int freq = cladeSet.getFrequency(i);
			double prop = ((double)freq) / totalTrees;
			if (prop > 0.5) {
				System.out.print(freq);
				System.out.print("\t" + nf.formatDecimal(prop * 100.0, 2) + "%");
				System.out.println("\t" + cladeSet.getClade(i));
			}
		}
							
		System.out.flush();
	}

    public void shortReport(String name, Tree tree, boolean drawHeader) throws IOException {

        String targetTree = "";
        if (tree != null) targetTree = Tree.Utils.uniqueNewick(tree, tree.getRoot());

        int n = treeSet.size();
        int totalTrees = treeSet.getSumFrequency();
        double highestProp = ((double)treeSet.getFrequency(0)) / totalTrees;
        String mapTree = (String)treeSet.get(0);

        if (drawHeader) {
            System.out.println("file\ttrees\tuniqueTrees\tp(MAP)\tMAP tree\t95credSize\ttrue_I\tp(true)\tcum(true)");
        }

        System.out.print(name+"\t");
        System.out.print(totalTrees+"\t");
        System.out.print(n+"\t");
        System.out.print(highestProp +"\t");
        System.out.print(mapTree +"\t");

        int credSet = (95 * totalTrees) / 100;
        int sumFreq = 0;

        int credSetSize = -1;
        int targetTreeIndex = -1;
        double targetTreeProb = 0.0;
        double targetTreeCum = 1.0;
        for (int i = 0; i < n; i++) {
            int freq = treeSet.getFrequency(i);
            double prop = ((double)freq) / totalTrees;

            sumFreq += freq;
            double sumProp = ((double)sumFreq) / totalTrees;

            String newickTree = (String)treeSet.get(i);

            if (newickTree.equals(targetTree)) {
                targetTreeIndex = i + 1;
                targetTreeProb = prop;
                targetTreeCum = sumProp;
            }

            if (sumFreq >= credSet) {
                if (credSetSize == -1) credSetSize = i + 1;
            }
        }

        System.out.print(credSetSize +"\t");
        System.out.print(targetTreeIndex +"\t");
        System.out.print(targetTreeProb +"\t");
        System.out.println(targetTreeCum);
    }

	public int getBurnin() { return burnin; }
	
	/**
	 * @return an analyses of the trees in a log file.
	 */
	private static TreeTraceAnalysis analyzeLogFile(String fileName, int burnin) throws IOException {
		return analyzeLogFile(new Reader[] {new FileReader(fileName)}, burnin, true);
	}
	
	/**
	 * @return an analyses of the trees in a log file.
	 */
	public static TreeTraceAnalysis analyzeLogFile(Reader[] reader, int burnin, boolean verbose) throws IOException {
		
		TreeTrace[] trace = new TreeTrace[reader.length];
		for (int i = 0; i < reader.length; i++) {
            try {
                trace[i] = TreeTrace.loadTreeTrace(reader[i]);
            } catch (Importer.ImportException ie) {
                throw new RuntimeException(ie.toString());
            }
            reader[i].close();

        }

		TreeTraceAnalysis analysis = new TreeTraceAnalysis(trace, burnin, verbose);
		
		return analysis;
	}
		
	private int burnin = -1;
	private TreeTrace[] trace;
	
	private CladeSet cladeSet;
	private FrequencySet treeSet; 
}