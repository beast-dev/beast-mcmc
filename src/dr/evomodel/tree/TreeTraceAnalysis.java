/*
 * TreeTraceAnalysis.java
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

package dr.evomodel.tree;

import dr.app.tools.NexusExporter;
import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.TreeTrace;
import dr.evolution.tree.*;
import dr.util.FrequencySet;
import dr.util.NumberFormatter;
import jebl.evolution.treemetrics.RobinsonsFouldMetric;

import java.io.IOException;
import java.io.PrintStream;
import java.io.Reader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: TreeTraceAnalysis.java,v 1.20 2005/06/07 16:28:18 alexei Exp $
 */
public class TreeTraceAnalysis {

    private TreeTraceAnalysis(TreeTrace[] traces, int burnIn, boolean verbose) {

        this.traces = traces;

        int minMaxState = Integer.MAX_VALUE;
        for (TreeTrace trace : traces) {
            if (trace.getMaximumState() < minMaxState) {
                minMaxState = trace.getMaximumState();
            }
        }

        if (burnIn < 0 || burnIn >= minMaxState) {
            this.burnin = minMaxState / (10 * traces[0].getStepSize());
            if (verbose) {
                String reason = burnIn < 0 ? "Defalt burn-in" : "WARNING: Burn-in larger than total number of states";
                System.out.println(reason + " - using 10% of smallest trace");
            }
        } else {
            this.burnin = burnIn;
        }

        analyze(verbose);
    }

    public static double[] getSymmetricTreeDistanceTrace(TreeTrace treeTrace, Tree targetTree) {
        double[] symDistance = new double[treeTrace.getTreeCount(0)];

        RobinsonsFouldMetric metric = new RobinsonsFouldMetric();
        jebl.evolution.trees.RootedTree jreference = TreeUtils.asJeblTree(targetTree);
        for (int i = 0; i < symDistance.length; i++) {
            final jebl.evolution.trees.RootedTree tree = TreeUtils.asJeblTree(treeTrace.getTree(i, 0));

            symDistance[i] = metric.getMetric(jreference, tree);
        }
        return symDistance;
    }


    /**
     * Actually analyzes the trace given the burnin
     *
     * @param verbose if true then progress is logged to stdout
     */
    void analyze(boolean verbose) {

        if (verbose) {
            if (traces.length > 1) System.out.println("Combining " + traces.length + " traces.");
        }

        final Tree tree0 = getTree(0);

        double[][] changed = new double[tree0.getNodeCount()][tree0.getNodeCount()];
        double[] rateConditionalOnChange = new double[tree0.getNodeCount()];
        boolean changesFound = false;

        cladeSet = new CladeSet(tree0);
        treeSet = new FrequencySet<String>();
        treeSet.add(TreeUtils.uniqueNewick(tree0, tree0.getRoot()));

        final int reportRate = 60;

        for (TreeTrace trace : traces) {
            final int treeCount = trace.getTreeCount(burnin * trace.getStepSize());
            final double stepSize = treeCount / (double) reportRate;
            int counter = 1;

            if (verbose) {
                System.out.println("Analyzing " + treeCount + " trees...");
                System.out.println("0              25             50             75            100");
                System.out.println("|--------------|--------------|--------------|--------------|");
                System.out.print("*");
            }

            for (int i = 1; i < treeCount; i++) {
                Tree tree = trace.getTree(i, burnin * trace.getStepSize());
                for (int j = 0; j < tree.getNodeCount(); j++) {
                    if (tree.getNode(j) != tree.getRoot() && tree.getNodeAttribute(tree.getNode(j), "changed") != null) {
                        changesFound = true;
                        final Object o = tree.getNodeAttribute(tree.getNode(j), "changed");
                        if (o != null) {
                            boolean ch = getChanged(tree, j);
                            if (ch) {
                                rateConditionalOnChange[j] += (Double) tree.getNodeAttribute(tree.getNode(j), "rate");
                            }
                            for (int k = 0; k < tree.getNodeCount(); k++) {
                                if (tree.getNode(k) != tree.getRoot()) {

                                    changed[j][k] += (ch && getChanged(tree, k)) ? 1 : 0;
                                }
                            }
                        }
                    }
                }

                cladeSet.add(tree);
                treeSet.add(TreeUtils.uniqueNewick(tree, tree.getRoot()));

                if (verbose && i >= (int) Math.round(counter * stepSize) && counter <= reportRate) {
                    System.out.print("*");
                    System.out.flush();

                    counter += 1;
                }
            }
            if (verbose) {
                System.out.println("*");
            }
        }
        if (changesFound) {
            for (int j = 0; j < tree0.getNodeCount(); j++) {
                System.out.println(j + "\t" + rateConditionalOnChange[j]);
            }
            System.out.println();
            for (int j = 0; j < tree0.getNodeCount(); j++) {
                for (int k = 0; k < tree0.getNodeCount(); k++) {
                    System.out.print(changed[j][k] + "\t");
                }
                System.out.println();
            }

        }
    }

    private boolean getChanged(Tree tree, int j) {
        final Object o = tree.getNodeAttribute(tree.getNode(j), "changed");
        if (o instanceof Integer) return (Integer) o == 1;
        return (Boolean) o;
    }

    /**
     * Actually analyzes a particular tree using the trace given the burnin
     *
     * @param target a tree in uniqueNewick format
     * @return a tree with mean node heights
     */
    final MutableTree analyzeTree(String target) {

        final int n = getTreeCount();

        FlexibleTree meanTree = null;

        // todo using CladeSet may probably speed this a lot
        for (int i = 0; i < n; i++) {
            final Tree tree = getTree(i);

            if (TreeUtils.uniqueNewick(tree, tree.getRoot()).equals(target)) {
                meanTree = new FlexibleTree(tree);
                break;
            }
        }
        if (meanTree == null) {
            throw new RuntimeException("No target tree in trace");
        }

        final int inc = meanTree.getInternalNodeCount();
        for (int j = 0; j < inc; j++) {
            double[] heights = new double[n];
            NodeRef nodej = meanTree.getInternalNode(j);
            Set<String> leafSet = TreeUtils.getDescendantLeaves(meanTree, nodej);

            for (int i = 0; i < n; i++) {
                final Tree tree = getTree(i);

                NodeRef can = TreeUtils.getCommonAncestorNode(tree, leafSet);
                heights[i] = tree.getNodeHeight(can);
            }

            meanTree.setNodeHeight(nodej, dr.stats.DiscreteStatistics.mean(heights));
            final double upper = dr.stats.DiscreteStatistics.quantile(0.975, heights);
            meanTree.setNodeAttribute(nodej, "upper", upper);
            final double lower = dr.stats.DiscreteStatistics.quantile(0.025, heights);
            meanTree.setNodeAttribute(nodej, "lower", lower);
            // Make it possible to display bars in figtree
            meanTree.setNodeAttribute(nodej, "range", new Double[]{lower, upper});
        }
        return meanTree;
    }

    final int getTreeCount() {

        int treeCount = 0;
        for (TreeTrace trace : traces) {
            treeCount += trace.getTreeCount(burnin * trace.getStepSize());
        }
        return treeCount;
    }

    final Tree getTree(int index) {

        int oldTreeCount = 0;
        int newTreeCount = 0;
        for (TreeTrace trace : traces) {
            final int br = burnin * trace.getStepSize();
            newTreeCount += trace.getTreeCount(br);

            if (index < newTreeCount) {
                return trace.getTree(index - oldTreeCount, br);
            }
            oldTreeCount = newTreeCount;
        }
        throw new RuntimeException("Couldn't find tree " + index);
    }

    public void report(int minNT) throws IOException {
        report(0.5, 0.95, minNT);
    }

    public void report(double minCladeProbability, int minNT) throws IOException {
        report(minCladeProbability, 0.95, minNT);
    }

    /**
     * @param minCladeProbability clades with at least this posterior probability will be included in report.
     * @throws IOException if general I/O error occurs
     */
    public void report(double minCladeProbability, double credSetProbability, int minNT) throws IOException {

        System.err.println("making report");

        final int fieldWidth = 14;
        NumberFormatter formatter = new NumberFormatter(6);
        formatter.setPadding(true);
        formatter.setFieldWidth(fieldWidth);

        final int nTreeSet = treeSet.size();
        int totalTrees = treeSet.getSumFrequency();

        System.out.println();
        System.out.println("burnIn=" + burnin);
        System.out.println("total trees used =" + totalTrees);
        System.out.println();


        System.out.println((Math.round(credSetProbability * 100.0))
                + "% credible set (" + nTreeSet + " unique trees, " + totalTrees + " total):");

        System.out.println("Count\tPercent\tTree");
        int credSet = (int) (credSetProbability * totalTrees);
        int sumFreq = 0;
        int skipped = 0;

        NumberFormatter nf = new NumberFormatter(8);

        for (int i = 0; i < nTreeSet; i++) {
            final int freq = treeSet.getFrequency(i);
            boolean show = true;
            if( minNT > 0 && freq <= minNT ) {
                show = false;
                skipped += 1;
            }
            final double prop = ((double) freq) / totalTrees;
            if( show ) {
                System.out.print(freq);
                System.out.print("\t" + nf.formatDecimal(prop * 100.0, 2) + "%");
            }

            sumFreq += freq;
            final double sumProp = ((double) sumFreq) / totalTrees;
            if( show ) {
                System.out.print("\t" + nf.formatDecimal(sumProp * 100.0, 2) + "%");

                String newickTree = treeSet.get(i);

                if (freq > 100) {
                    // calculate conditional average node heights
                    Tree meanTree = analyzeTree(newickTree);
                    System.out.println("\t" + TreeUtils.newick(meanTree));

                } else {
                    System.out.println("\t" + newickTree);
                }
            }

            if (sumFreq >= credSet) {
                if( skipped > 0 ) {
                   System.out.println();
                   System.out.println("... (" + skipped + ") trees.");  
                }
                System.out.println();
                System.out.println("95% credible set has " + (i + 1) + " trees.");
                break;
            }
        }

        System.out.println();
        System.out.println(Math.round(minCladeProbability * 100.0) +
                "%-rule clades (" + cladeSet.size() + " unique clades):");
        final int nCladeSet = cladeSet.size();
        for (int i = 0; i < nCladeSet; i++) {
            final int freq = cladeSet.getFrequency(i);
            final double prop = ((double) freq) / totalTrees;
            if (prop >= minCladeProbability) {
                System.out.print(freq);
                System.out.print("\t" + nf.formatDecimal(prop * 100.0, 2) + "%");
                System.out.print("\t" + cladeSet.getMeanNodeHeight(i));
                System.out.println("\t" + cladeSet.getClade(i));
            }
        }

        System.out.flush();

        System.out.println("Clade credible sets:");


        int fiveCredSet = (5 * totalTrees) / 100;
        int halfCredSet = (50 * totalTrees) / 100;
        sumFreq = 0;
        assert nTreeSet == treeSet.size();

        final CladeSet tempCladeSet = new CladeSet();
        for (int nt = 0; nt < nTreeSet; nt++) {

            sumFreq += treeSet.getFrequency(nt);

            String newickTree = treeSet.get(nt);
            NewickImporter importer = new NewickImporter(new StringReader(newickTree));

            try {
                Tree tree = importer.importNextTree();

                tempCladeSet.add(tree);
            } catch (Importer.ImportException e) {
                System.err.println("Err");
            }

            if (sumFreq >= fiveCredSet) {
                System.out.println();
                System.out.println("5% credible set has " + tempCladeSet.getCladeCount() + " clades.");
                // don't do it more than once
                fiveCredSet = totalTrees + 1;
            }

            if (sumFreq >= halfCredSet) {
                System.out.println();
                System.out.println("50% credible set has " + tempCladeSet.getCladeCount() + " clades.");
                // don't do it more than once
                halfCredSet = totalTrees + 1;
            }
        }

        System.out.flush();
    }

    public void shortReport(String name, Tree tree, boolean drawHeader) {
        shortReport(name, tree, drawHeader, 0.95);
    }

    public void shortReport(String name, Tree tree, boolean drawHeader, double credSetProbability) {

        String targetTree = "";
        if (tree != null) targetTree = TreeUtils.uniqueNewick(tree, tree.getRoot());

        final int n = treeSet.size();
        final int totalTrees = treeSet.getSumFrequency();
        final double highestProp = ((double) treeSet.getFrequency(0)) / totalTrees;
        String mapTree = treeSet.get(0);

        if (drawHeader) {
            System.out.println("file\ttrees\tuniqueTrees\tp(MAP)\tMAP tree\t" + (int) credSetProbability * 100 + "credSize\ttrue_I\tp(true)\tcum(true)");
        }

        System.out.print(name + "\t");
        System.out.print(totalTrees + "\t");
        System.out.print(n + "\t");
        System.out.print(highestProp + "\t");
        System.out.print(mapTree + "\t");

        int credSet = (int) (credSetProbability * totalTrees);
        int sumFreq = 0;

        int credSetSize = -1;
        int targetTreeIndex = -1;
        double targetTreeProb = 0.0;
        double targetTreeCum = 1.0;
        for (int i = 0; i < n; i++) {
            final int freq = treeSet.getFrequency(i);
            final double prop = ((double) freq) / totalTrees;

            sumFreq += freq;
            final double sumProp = ((double) sumFreq) / totalTrees;

            String newickTree = treeSet.get(i);

            if (newickTree.equals(targetTree)) {
                targetTreeIndex = i + 1;
                targetTreeProb = prop;
                targetTreeCum = sumProp;
            }

            if (sumFreq >= credSet) {
                if (credSetSize == -1) credSetSize = i + 1;
            }
        }

        System.out.print(credSetSize + "\t");
        System.out.print(targetTreeIndex + "\t");
        System.out.print(targetTreeProb + "\t");
        System.out.println(targetTreeCum);
    }

    public void export(PrintStream out, double minTreeProbability, int max, boolean verbose) {
        NexusExporter exporter = new NexusExporter(out);
        int n = treeSet.size();
        if (max < 0) max = n;

        final int totalTrees = treeSet.getSumFrequency();

        List<Tree> trees = new ArrayList<Tree>();

        final int totExport = Math.min(max, n);
        final boolean progress = verbose && totExport > 60;
        if (progress) {
            System.out.println("Exporting " + totExport + " trees...");
            System.out.println("0              25             50             75            100");
            System.out.println("|--------------|--------------|--------------|--------------|");
            System.out.print("*");
        }
        // todo have an option for threshold and sort by NCP. 
        for (int i = 0; i < n; i++) {
            int freq = treeSet.getFrequency(i);
            double prop = ((double) freq) / totalTrees;
            if (prop < minTreeProbability) {
                continue;
            }
            final String newickTree = treeSet.get(i);
            // calculate conditional average node heights
            final MutableTree tree = analyzeTree(newickTree);
            tree.setAttribute("weight", prop);
            double p = cladeSet.annotate(tree, "posterior");
            tree.setNodeAttribute(tree.getRoot(), "posterior", Math.exp(p / tree.getInternalNodeCount()));
            trees.add(tree);

            if (progress && ((i + 1) % (totExport / 60)) == 0) {
                System.out.print("*");
            }

            if (trees.size() == max) {
                break;
            }
        }

        if (trees.size() > 0) {
            exporter.exportTrees(trees.toArray(new Tree[trees.size()]), true);
        }
    }

    public int getBurnin() {
        return burnin;
    }

    /**
     * @param reader  the readers to be analyzed
     * @param burnin  the burnin in states
     * @param verbose true if progress should be logged to stdout
     * @return an analyses of the trees in a log file.
     * @throws java.io.IOException if general I/O error occurs
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

        return new TreeTraceAnalysis(trace, burnin, verbose);
    }

    private int burnin = -1;
    private final TreeTrace[] traces;

    private CladeSet cladeSet;
    private FrequencySet<String> treeSet;
}