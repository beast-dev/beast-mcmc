/*
 * TreeAnnotator.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.tools.treeannotator;

import dr.app.beast.BeastVersion;
import dr.app.tools.BaseTreeTool;
import dr.app.tools.NexusExporter;
import dr.app.tools.logcombiner.LogCombiner;
import dr.app.util.Arguments;
import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.treemetrics.*;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.util.Version;
import jam.console.ConsoleApplication;

import javax.swing.*;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @author Philippe Lemey
 * @author Guy Baele
 * @author Alexei Drummond
 */
public class TreeAnnotator extends BaseTreeTool {

    private static final Version VERSION = new BeastVersion();

    private static final HeightsSummary DEFAULT_HEIGHTS_SUMMARY = HeightsSummary.MEAN_HEIGHTS;
    private static final boolean COUNT_TREES = true;

    private static final boolean THREADED_READING = true;

    // Messages to stderr, output to stdout
    private static PrintStream progressStream = System.err;
    private static final boolean extendedMetrics = false;

    private final CollectionAction collectionAction;
    private final AnnotationAction annotationAction;
    private final int threadCount;

    private TaxonList taxa = null;
    private int totalTrees;
    private int totalTreesUsed;

    List<Double> rootHeights = new ArrayList<>();

    enum Target {
        HIPSTR("Highest independent posterior subtree reconstruction (HIPSTR)"),
        MRHIPSTR("Majority rule highest independent posterior subtree reconstruction (MrHIPSTR)"),
        MAX_CLADE_CREDIBILITY("Maximum clade credibility tree"),
        MAJORITY_RULE("Majority-rule consensus tree"),
        USER_TARGET_TREE("User target tree");

        String desc;

        Target(String s) {
            desc = s;
        }

        public String toString() {
            return desc;
        }
    }

    enum HeightsSummary {
        MEAN_HEIGHTS("Mean heights"),
        MEDIAN_HEIGHTS("Median heights"),
        KEEP_HEIGHTS("Keep target heights");

        String desc;

        HeightsSummary(String s) {
            desc = s;
        }

        public String toString() {
            return desc;
        }
    }

    /**
     * Burnin can be specified as the number of trees or the number of states
     * (one or other should be zero).
     */
    public TreeAnnotator(final int burninTrees,
                         final long burninStates,
                         final HeightsSummary heightsOption,
                         final double posteriorLimit,
                         final int countLimit,
                         final double[] hpd2D,
                         final boolean computeESS,
                         final int threadCount,
                         final Target targetOption,
                         final String targetTreeFileName,
                         final String referenceTreeFileName,
                         final String treeMetricFileName,
                         final String inputFileName,
                         final String outputFileName
    ) throws IOException {

        long totalStartTime = System.currentTimeMillis();

        collectionAction = new CollectionAction();

        collectionAction.addAttributeName("height");
        collectionAction.addAttributeName("length");

        annotationAction = new AnnotationAction(heightsOption, posteriorLimit, countLimit, hpd2D, computeESS, true);

        annotationAction.addAttributeName("height");
        annotationAction.addAttributeName("length");

        int burnin = -1;

        totalTrees = 10000;
        totalTreesUsed = 0;

        this.threadCount = threadCount;

        CladeSystem cladeSystem = new CladeSystem(targetOption == Target.HIPSTR || targetOption == Target.MRHIPSTR || targetOption == Target.MAJORITY_RULE, targetOption == Target.MAJORITY_RULE);

        if (COUNT_TREES) {
            countTrees(inputFileName);
            progressStream.println("Reading trees...");
        } else {
            totalTrees = 10000;
            progressStream.println("Reading trees (assuming 10,000 trees)...");
        }

        burnin = readTrees(inputFileName, burninTrees, burninStates, cladeSystem);

        cladeSystem.calculateCladeCredibilities(totalTreesUsed);

        progressStream.println("Total trees read: " + totalTrees);
        progressStream.println("Size of trees: " + taxa.getTaxonCount() + " tips");
        if (burninTrees > 0) {
            progressStream.println("Ignoring first " + burninTrees + " trees" +
                    (burninStates > 0 ? " (" + burninStates + " states)." : "." ));
        } else if (burninStates > 0) {
            progressStream.println("Ignoring first " + burninStates + " states (" + burnin + " trees).");
        }

        printCladeInformation(cladeSystem);
//        }

        MutableTree targetTree = null;

        switch (targetOption) {
            case USER_TARGET_TREE: {
                if (targetTreeFileName != null) {
                    targetTree = readUserTargetTree(targetTreeFileName, cladeSystem);
                } else {
                    System.err.println("No user target tree specified.");
                    System.exit(1);
                }
                break;
            }
            case MAX_CLADE_CREDIBILITY: {
                progressStream.println("Finding maximum credibility tree...");
                targetTree = new FlexibleTree(getMCCTree(burnin, cladeSystem, inputFileName));
                break;
            }
            case HIPSTR: {
                progressStream.println("Finding highest independent posterior subtree reconstruction (HIPSTR) tree...");
                targetTree = getHIPSTRTree(cladeSystem, false);
                break;
            }
            case MRHIPSTR: {
                progressStream.println("Finding majority rule highest independent posterior subtree reconstruction (MrHIPSTR) tree...");
                targetTree = getHIPSTRTree(cladeSystem, true);
                break;
            }
            case MAJORITY_RULE: {
                progressStream.println("Finding majority-rule consensus tree...");
                targetTree = getMajorityRuleConsensusTree(cladeSystem);
                break;
            }
            default: throw new IllegalArgumentException("Unknown targetOption");
        }


        if (referenceTreeFileName != null) {
            CladeSystem targetCladeSystem = new CladeSystem(targetTree);

            progressStream.println("Reading reference tree: " + referenceTreeFileName);

            MutableTree referenceTree = readTreeFile(referenceTreeFileName);
            CladeSystem referenceCladeSystem = new CladeSystem(referenceTree);

            int commonCladeCount = targetCladeSystem.getCommonCladeCount(referenceCladeSystem);
            progressStream.println("Clades in common with reference tree: " + commonCladeCount +
                    " (out of " + referenceCladeSystem.getCladeCount() + ")");
            progressStream.println();
        }

        collectNodeAttributes(cladeSystem, inputFileName, burnin);

        annotateTargetTree(cladeSystem, heightsOption, countLimit, targetTree);

        writeAnnotatedTree(outputFileName, targetTree);

        if (treeMetricFileName != null && !treeMetricFileName.isEmpty()) {
            writeTreeMetrics(burnin, targetTree, inputFileName, treeMetricFileName);
        }

        long timeElapsed =  (System.currentTimeMillis() - totalStartTime) / 1000;
        progressStream.println("Total time: " + timeElapsed + " secs");
        progressStream.println();

    }

    private static void printCladeInformation(CladeSystem cladeSystem) {
        int n = cladeSystem.getCladeCount();
        progressStream.println("Total unique clades: " + n);
        progressStream.println("Total clades in more than one tree: " + (n - cladeSystem.getCladeFrequencyCount(1)));
        progressStream.println();
    }

    private void countTrees(String inputFileName) throws IOException {
        progressStream.println("Counting trees...");
        Reader reader = new BufferedReader(new FileReader(inputFileName));
        TreeImporter importer = new BEASTTreesImporter(reader);
        totalTrees = importer.countTrees();
        if (totalTrees == 0) {
            totalTrees = 10000;
        }
        reader.close();
        progressStream.println("Total number of trees: " + totalTrees);
        progressStream.println();
    }

    private int readTrees(String inputFileName, int burninTrees, long burninStates, CladeSystem cladeSystem) throws IOException {
        long timeElapsed;
        long startTime;

        int burnin = -1;

        long stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        startTime = System.currentTimeMillis();

        try {
            // read the first tree using NexusImport to get the taxon list and tip number to taxon mapping
//            Reader reader = new BufferedReader(new FileReader(inputFileName));
//            NexusImporter nexusImporter = new NexusImporter(reader, true);
//            taxa = nexusImporter.importTree(null);
//
//            reader = new BufferedReader(new FileReader(inputFileName));
//            BEASTTreesImporter importer = new BEASTTreesImporter(reader, false);
//            importer.setTaxonList(taxa);

            Reader reader = new BufferedReader(new FileReader(inputFileName));
            NexusImporter importer = new NexusImporter(reader, true);

            ExecutorService pool;
            List<Future<?>> futures;
            if (THREADED_READING) {
                if (threadCount <= 0) {
                    pool = Executors.newCachedThreadPool();
                } else {
                    pool = Executors.newFixedThreadPool(threadCount);
                }
                futures = new ArrayList<>();
            }
            totalTrees = 0;
            boolean firstTree = true;
            while (importer.hasTree()) {
                final Tree tree = importer.importNextTree();
                long state = 0;

                if (taxa == null) {
                    taxa = new Taxa(tree);
                }

                if (burninStates > 0) {
                    // if burnin has been specified in states, try to parse it out...
                    String name = tree.getId().trim();

                    long maxState;
                    if (name.startsWith("STATE_")) {
                        state = Long.parseLong(name.split("_")[1]);
                        maxState = state;
                    } else {
                        maxState = state;
                        state += 1;
                    }
                }

                if (totalTrees >= burninTrees && state >= burninStates) {
                    // if either of the two burnin thresholds have been reached...

                    if (burnin < 0) {
                        // if this is the first time this point has been reached,
                        // record the number of trees this represents for future use...
                        burnin = totalTrees;
                    }

                    if (firstTree) {
                        // for the first tree do it outside a thread
                        cladeSystem.add(tree);
                        firstTree = false;
                    } else {
                        if (THREADED_READING) {
                            futures.add(pool.submit(() -> {
                                cladeSystem.add(tree);
                            }));
                        } else {
                            cladeSystem.add(tree);
                        }
                    }
                    totalTreesUsed += 1;
                }

                if (totalTrees > 0 && totalTrees % stepSize == 0) {
                    progressStream.print("*");
                    progressStream.flush();
                }
                totalTrees++;
            }

            if (THREADED_READING) {
                try {
                    // wait for all the threads to run to completion
                    for (Future<?> f : futures) {
                        f.get();
                    }
                } catch (InterruptedException | ExecutionException e) {
                    throw new RuntimeException(e);
                }
            }

            reader.close();
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            System.exit(1);
        }
        timeElapsed =  (System.currentTimeMillis() - startTime) / 1000;
        progressStream.println("* [" + timeElapsed + " secs]");
        progressStream.println();

        if (totalTrees < 1) {
            System.err.println("No trees in input file");
            System.exit(1);
        }
        if (totalTreesUsed < 1) {
            if (burninTrees > 0 || burninStates > 0) {
                System.err.println("No trees to use: burnin greater than number of trees in input file");
                System.exit(1);
            }
        }
        return burnin;
    }

    private void collectNodeAttributes(CladeSystem cladeSystem, String inputFileName, int burnin) throws IOException {
        progressStream.println("Collecting node information...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        Reader reader = new BufferedReader(new FileReader(inputFileName));
//         TreeImporter importer = new BEASTTreesImporter(reader, true);
        TreeImporter importer = new NexusImporter(reader, false);

        long startTime = System.currentTimeMillis();

        totalTreesUsed = 0;
        try {
            ExecutorService pool;
            List<Future<?>> futures;
            if (THREADED_READING) {
                if (threadCount <= 0) {
                    pool = Executors.newCachedThreadPool();
                } else {
                    pool = Executors.newFixedThreadPool(threadCount);
                }

                futures = new ArrayList<>();
            }

            boolean firstTree = true;
            int counter = 0;

            while (importer.hasTree()) {
                final Tree tree = importer.importNextTree();

                if (counter >= burnin) {
                    if (firstTree) {
                        setupAttributes(tree);
                        firstTree = false;
                    }

                    if (THREADED_READING) {
                        futures.add(pool.submit(() -> {
                            cladeSystem.collectCladeHeights(tree);
                            rootHeights.add(tree.getNodeHeight(tree.getRoot()));
                            cladeSystem.traverseTree(tree, collectionAction);
                        }));
                    } else {
                        cladeSystem.collectCladeHeights(tree);
                        rootHeights.add(tree.getNodeHeight(tree.getRoot()));
                        cladeSystem.traverseTree(tree, collectionAction);
                    }
                    totalTreesUsed += 1;
                }
                if (counter > 0 && counter % stepSize == 0) {
                    progressStream.print("*");
                    progressStream.flush();
                }
                counter++;

            }

            if (THREADED_READING) {
                try {
                    // wait for all the threads to run to completion
                    for (Future<?> f : futures) {
                        f.get();
                    }
                } catch (InterruptedException | ExecutionException e) {
                    throw new RuntimeException(e);
                }
            }

            cladeSystem.calculateCladeCredibilities(totalTreesUsed);
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            System.exit(1);
        }
        long timeElapsed =  (System.currentTimeMillis() - startTime) / 1000;
        progressStream.println("* [" + timeElapsed + " secs]");
        progressStream.println();
        reader.close();
    }

    public void setupAttributes(Tree tree) {
        Set<String> attributeNames = new TreeSet<>(String.CASE_INSENSITIVE_ORDER);
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
        collectionAction.addAttributeNames(attributeNames);
        annotationAction.addAttributeNames(attributeNames);
    }

    private MutableTree readUserTargetTree(String targetTreeFileName, CladeSystem cladeSystem) throws IOException {
        progressStream.println("Reading user specified target tree, " + targetTreeFileName + ", ...");

        MutableTree targetTree = readTreeFile(targetTreeFileName);

        progressStream.println();
        double score = scoreTree(targetTree, cladeSystem);
        progressStream.println("Target tree's log clade credibility: " + String.format("%.4f", score));
        reportStatistics(cladeSystem, targetTree);
//        reportStatisticTables(cladeSystem, targetTree);

        progressStream.println();
        return targetTree;
    }


    private static MutableTree readTreeFile(String treeFileName) throws IOException {
        NexusImporter importer = new NexusImporter(new FileReader(treeFileName));
        Tree tree = null;
        try {
            tree = importer.importNextTree();
            if (tree == null) {
                NewickImporter x = new NewickImporter(new FileReader(treeFileName));
                tree = x.importNextTree();
            }
            if (tree == null) {
                System.err.println("No tree in nexus or newick file " + treeFileName);
                System.exit(1);
            }
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Target Tree: " + e.getMessage());
            System.exit(1);
        }
        return new FlexibleTree(tree);
    }

    private Tree getMCCTree(int burnin, CladeSystem cladeSystem, String inputFileName)
            throws IOException {

        long startTime = System.currentTimeMillis();

        Tree bestTree = null;
        double bestScore = Double.NEGATIVE_INFINITY;

        progressStream.println("Analyzing " + totalTreesUsed + " trees...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        int counter = 0;
        int bestTreeNumber = 0;
//        TreeImporter importer = new BEASTTreesImporter(new FileReader(inputFileName), false);
        TreeImporter importer = new NexusImporter(new FileReader(inputFileName), true);
        try {
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (counter >= burnin) {
                    double score = scoreTree(tree, cladeSystem);
//                    progressStream.println(score);
                    if (score > bestScore) {
                        bestTree = tree;
                        bestScore = score;
                        bestTreeNumber = counter + 1;
                    }
                }
                if (counter > 0 && counter % stepSize == 0) {
                    progressStream.print("*");
                    progressStream.flush();
                }
                counter++;
            }
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            System.exit(1);
        }

        long timeElapsed =  (System.currentTimeMillis() - startTime) / 1000;
        progressStream.println("* [" + timeElapsed + " secs]");
        progressStream.println();
        progressStream.println("Best tree: " + bestTree.getId() + " (tree number " + bestTreeNumber + ")");
        progressStream.println("Best tree's log clade credibility: " + String.format("%.4f", bestScore));
        reportStatistics(cladeSystem, bestTree);
//        reportStatisticTables(cladeSystem, bestTree);
        progressStream.println();

        return bestTree;
    }

    private MutableTree getMajorityRuleConsensusTree(CladeSystem cladeSystem) {

        long startTime = System.currentTimeMillis();

        MajorityRuleTreeBuilder treeBuilder = new MajorityRuleTreeBuilder();
        MutableTree tree = treeBuilder.getMajorityRuleConsensusTree(cladeSystem, taxa);

        // majority rule tree may be non-bifurcating
        // double score = scoreTree(tree, cladeSystem);

        double timeElapsed =  (double)(System.currentTimeMillis() - startTime) / 1000;
        progressStream.println("[" + timeElapsed + " secs]");
        progressStream.println();
//        progressStream.println("Majority rule consensus tree's log clade credibility: " + String.format("%.4f", score));
        reportStatistics(cladeSystem, tree);
//        reportStatisticTables(cladeSystem, tree);
        progressStream.println();

        return tree;
    }

    private MutableTree getHIPSTRTree(CladeSystem cladeSystem, boolean majorityRule) {

        long startTime = System.currentTimeMillis();

        HIPSTRTreeBuilder treeBuilder = new HIPSTRTreeBuilder();
        MutableTree tree = treeBuilder.getHIPSTRTree(cladeSystem, taxa, majorityRule);

        double score = scoreTree(tree, cladeSystem);

        double timeElapsed =  (double)(System.currentTimeMillis() - startTime) / 1000;
        progressStream.println("[" + timeElapsed + " secs]");
        progressStream.println();
        if (majorityRule) {
            progressStream.println("MrHIPSTR tree's log clade credibility: " + String.format("%.4f", score));
        } else {
            progressStream.println("HIPSTR tree's log clade credibility: " + String.format("%.4f", score));
        }
        reportStatistics(cladeSystem, tree);
//        reportStatisticTables(cladeSystem, tree);
        progressStream.println();

        return tree;
    }

    private static void reportStatistics(CladeSystem cladeSystem, Tree tree) {
        progressStream.println("Lowest individual clade credibility: " + String.format("%.4f", cladeSystem.getMinimumCladeCredibility(tree)));
        progressStream.println("Mean individual clade credibility: " + String.format("%.4f", cladeSystem.getMeanCladeCredibility(tree)));
        progressStream.println("Median individual clade credibility: " + String.format("%.4f", cladeSystem.getMedianCladeCredibility(tree)));
        progressStream.println("Number of clades with credibility 1.0: " + cladeSystem.getTopCladeCount(tree, 1.0));
        reportCladeCredibilityCount(cladeSystem, tree, 0.99, extendedMetrics);
        reportCladeCredibilityCount(cladeSystem, tree, 0.95, extendedMetrics);
        if (extendedMetrics) {
            reportCladeCredibilityCount(cladeSystem, tree, 0.9, extendedMetrics);
            reportCladeCredibilityCount(cladeSystem, tree, 0.8, extendedMetrics);
            reportCladeCredibilityCount(cladeSystem, tree, 0.7, extendedMetrics);
            reportCladeCredibilityCount(cladeSystem, tree, 0.6, extendedMetrics);
        }
        reportCladeCredibilityCount(cladeSystem, tree, 0.5, extendedMetrics, true);
        if (extendedMetrics) {
            reportCladeCredibilityCount(cladeSystem, tree, 0.25, extendedMetrics);
            reportCladeCredibilityCount(cladeSystem, tree, 0.1, extendedMetrics);
            reportCladeCredibilityCount(cladeSystem, tree, 0.05, extendedMetrics);
        }
    }

    private static void reportCladeCredibilityCount(CladeSystem cladeSystem, Tree tree, double threshold, boolean extendedMetrics) {
        reportCladeCredibilityCount(cladeSystem, tree, threshold, extendedMetrics, false);
    }

    private static void reportCladeCredibilityCount(CladeSystem cladeSystem, Tree tree, double threshold, boolean extendedMetrics, boolean showMissingClades) {
        int treeCladeCount = cladeSystem.getTopCladeCount(tree, threshold);
        int allCladeCount = cladeSystem.getTopCladeCount(threshold);
        progressStream.print("Number of clades with credibility > " + threshold + ": " +
                treeCladeCount);
        if (treeCladeCount < allCladeCount) {
            if (extendedMetrics) {
                Set<BiClade> treeClades = cladeSystem.getTopClades(tree, threshold);
                Set<BiClade> allClades = cladeSystem.getTopClades(threshold);

                Set<BiClade> missingClades = new HashSet<>(allClades);
                missingClades.removeAll(treeClades);
                double sum = 0;
                int max = 0;
                for (BiClade missing : missingClades) {
                    sum += missing.getSize();
                    if (missing.getSize() > max) {
                        max = missing.getSize();
                    }
                }
                progressStream.print(" / " + allCladeCount + " | missing: " + missingClades.size());
                progressStream.printf(",  mean size: %.2f", (sum / missingClades.size()));
                progressStream.printf(",  max size: %d", max);
                if (showMissingClades) {
                    progressStream.println();
                    for (BiClade missing : missingClades) {
                        progressStream.println(" (" + missing.getSize() + ", " + missing.getCredibility() + ", {" + missing + "} )");
                    }
                }

            } else {
                progressStream.print(" / " + allCladeCount + " (in all trees)");
            }
        }
        progressStream.println();
    }

    private static void reportStatisticTables(CladeSystem cladeSystem, Tree tree) {
        int count = 100;
//        double[] table = new double[count + 1];
//        for (int i = 0; i <= count; i++) {
//            double threshold = ((double) (i)) / count;
//            table[i] = cladeSystem.getTopCladeCredibility(tree, threshold);
//        }

        progressStream.println("threshold, #clades");
        for (int i = 0; i <= count; i++) {
            double threshold = ((double) (i)) / count;
            progressStream.print(threshold);
            progressStream.print(",");
            progressStream.print(cladeSystem.getTopCladeCount(tree, threshold));
            progressStream.print(",");
            progressStream.println(cladeSystem.getTopCladeCount(threshold));
        }

    }

    private void writeTreeMetrics(int burnin, Tree referenceTree, String inputFileName, String outputFileName)
            throws IOException {

        long startTime = System.currentTimeMillis();

        TreeMetric[] treeMetrics = new TreeMetric[] {
                new BranchScoreMetric(),
                new CladeHeightMetric(),
                new KendallColijnPathDifferenceMetric(0.0),
                new KendallColijnPathDifferenceMetric(0.5),
                new KendallColijnPathDifferenceMetric(1.0),
                new RobinsonFouldsMetric(),
                new RootedBranchScoreMetric(),
                new SteelPennyPathDifferenceMetric()
        };

        try {
            final PrintStream stream = outputFileName != null ?
                    new PrintStream(Files.newOutputStream(Paths.get(outputFileName))) :
                    System.out;


            progressStream.println("Writing tree metrics for " + totalTreesUsed + " trees...");
            progressStream.println("0              25             50             75            100");
            progressStream.println("|--------------|--------------|--------------|--------------|");

            int stepSize = totalTrees / 60;
            if (stepSize < 1) stepSize = 1;

            int counter = 0;
            TreeImporter importer = new NexusImporter(new FileReader(inputFileName), true);
            try {
                stream.print("tree");
                for (TreeMetric treeMetric : treeMetrics) {
                    stream.print("\t" + treeMetric.getType().getName());
                }

                stream.println();

                while (importer.hasTree()) {
                    Tree tree = importer.importNextTree();

                    if (counter >= burnin) {
                        stream.print(counter);

                        for (TreeMetric treeMetric : treeMetrics) {
                            double score = treeMetric.getMetric(tree, referenceTree);

                            stream.print("\t" + score);
                        }

                        stream.println();
                    }

                    if (counter > 0 && counter % stepSize == 0) {
                        progressStream.print("*");
                        progressStream.flush();
                    }
                    counter++;
                }
            } catch (Importer.ImportException e) {
                System.err.println("Error Parsing Input Tree: " + e.getMessage());
                System.exit(1);
            }

            long timeElapsed =  (System.currentTimeMillis() - startTime) / 1000;
            progressStream.println("* [" + timeElapsed + " secs]");
            progressStream.println();
            progressStream.println("Tree metric comparisons to target tree written to file: " + outputFileName);
            progressStream.println();
        } catch (Exception e) {
            System.err.println("Error writing tree metric file: " + e.getMessage());
            System.exit(1);
        }
    }

    private void annotateTargetTree(CladeSystem cladeSystem, HeightsSummary heightsOption, int countLimit, MutableTree targetTree) {
        progressStream.println("Annotating target tree...");

        try {
            cladeSystem.traverseTree(targetTree, new SetHeightsAction(rootHeights, countLimit));

            cladeSystem.traverseTree(targetTree, annotationAction);
        } catch (Exception e) {
            System.err.println("Error annotating tree: " + e.getMessage() + "\nPlease check the tree log file format.");
            System.exit(1);
        }
        progressStream.println();
    }

    private void writeAnnotatedTree(String outputFileName, MutableTree targetTree) {
        progressStream.println("Writing annotated tree...");

        try {
            final PrintStream stream = outputFileName != null ?
                    new PrintStream(Files.newOutputStream(Paths.get(outputFileName))) :
                    System.out;

            new NexusExporter(stream).exportTree(targetTree);
        } catch (Exception e) {
            System.err.println("Error writing annotated tree file: " + e.getMessage());
            System.exit(1);
        }
        progressStream.println();
        progressStream.println("Written to file: " + outputFileName);
        progressStream.println();
    }


    private double scoreTree(Tree tree, CladeSystem cladeSystem) {
        return cladeSystem.getLogCladeCredibility(tree);
    }

    public static void printTitle() {
        progressStream.println();
        centreLine("TreeAnnotator " + VERSION.getVersionString() + ", " + VERSION.getDateString(), 60);
        centreLine("MCMC Output analysis", 60);
        centreLine("by", 60);
        centreLine("Andrew Rambaut, Marc A. Suchard and Alexei J. Drummond", 60);
        progressStream.println();
        centreLine("Institute of Ecology and Evolution", 60);
        centreLine("University of Edinburgh", 60);
        centreLine("a.rambaut@ed.ac.uk", 60);
        progressStream.println();
        centreLine("David Geffen School of Medicine", 60);
        centreLine("University of California, Los Angeles", 60);
        centreLine("msuchard@ucla.edu", 60);
        progressStream.println();
        centreLine("Department of Computer Science", 60);
        centreLine("University of Auckland", 60);
        centreLine("alexei@cs.auckland.ac.nz", 60);
        progressStream.println();
        progressStream.println();
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("treeannotator", "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: treeannotator test.trees out.tree");
        progressStream.println("  Example: treeannotator -burnin 100 -heights mean test.trees out.tree");
        progressStream.println("  Example: treeannotator -type hipstr -burnin 100 -heights mean test.trees out.tree");
        progressStream.println();
    }

    public static double[] parseVariableLengthDoubleArray(String inString) throws Arguments.ArgumentException {

        List<Double> returnList = new ArrayList<Double>();
        StringTokenizer st = new StringTokenizer(inString,",");
        while(st.hasMoreTokens()) {
            try {
                returnList.add(Double.parseDouble(st.nextToken()));
            } catch (NumberFormatException e) {
                throw new Arguments.ArgumentException();
            }

        }

        if (!returnList.isEmpty()) {
            double[] doubleArray = new double[returnList.size()];
            for(int i=0; i<doubleArray.length; i++)
                doubleArray[i] = returnList.get(i);
            return doubleArray;
        }
        return null;
    }

    //Main method
    public static void main(String[] args) throws IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        String targetTreeFileName = null;
        String referenceTreeFileName = null;
        String treeMetricFileName = null;
        String inputFileName = null;
        String outputFileName = null;

        boolean forceIntegerToDiscrete = false;
        boolean computeESS = false;

        if (args.length == 0) {
            System.setProperty("com.apple.macos.useScreenMenuBar", "true");
            System.setProperty("apple.laf.useScreenMenuBar", "true");
            System.setProperty("apple.awt.showGrowBox", "true");

            java.net.URL url = LogCombiner.class.getResource("/images/utility.png");
            Icon icon = null;

            if (url != null) {
                icon = new ImageIcon(url);
            }

            final String versionString = VERSION.getVersionString();
            String nameString = "TreeAnnotatorX " + versionString;
            String aboutString = "<html><center><p>" + versionString + ", " + VERSION.getDateString() + "</p>" +
                    "<p>by<br>" +
                    "Andrew Rambaut and Alexei J. Drummond</p>" +
                    "<p>Institute of Ecology and Evolution, University of Edinburgh<br>" +
                    "<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
                    "<p>Department of Computer Science, University of Auckland<br>" +
                    "<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +
                    "<p>Part of the BEAST package:<br>" +
                    "<a href=\"http://beast.community\">http://beast.community</a></p>" +
                    "</center></html>";

            new ConsoleApplication(nameString, aboutString, icon, true);

            // The ConsoleApplication will have overridden System.out so set progressStream
            // to capture the output to the window:
            progressStream = System.out;

            printTitle();

            TreeAnnotatorDialog dialog = new TreeAnnotatorDialog(new JFrame());

            if (!dialog.showDialog("TreeAnnotatorX " + versionString)) {
                return;
            }

            long burninStates = dialog.getBurninStates();
            int burninTrees = dialog.getBurninTrees();
            double posteriorLimit = dialog.getPosteriorLimit();
            double[] hpd2D = {0.80};
            Target targetOption = dialog.getTargetOption();
            HeightsSummary heightsOption = dialog.getHeightsOption();

            targetTreeFileName = dialog.getTargetFileName();
            if (targetOption == Target.USER_TARGET_TREE && targetTreeFileName == null) {
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
                new TreeAnnotator(
                        burninTrees,
                        burninStates,
                        heightsOption,
                        posteriorLimit,
                        5,
                        hpd2D,
                        computeESS,
                        -1,
                        targetOption,
                        targetTreeFileName,
                        referenceTreeFileName,
                        null,
                        inputFileName,
                        outputFileName);

            } catch (Exception ex) {
                System.err.println("Exception: " + ex.getMessage());
            }

            progressStream.println("Finished - Quit program to exit.");
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
                        new Arguments.StringOption("type", new String[] {"hipstr", "mrhipstr", "mcc", "mrc"}, false, "an option of 'hipstr' (default), 'mrhipstr', 'mcc' or 'mrc'"),
                        new Arguments.StringOption("heights", new String[] {"keep", "median", "mean", "ca"}, false,
                                "an option of 'keep', 'median' or 'mean' (default)"),
                        new Arguments.LongOption("burnin", "the number of states to be considered as 'burn-in'"),
                        new Arguments.IntegerOption("burninTrees", "the number of trees to be considered as 'burn-in'"),
                        new Arguments.RealOption("limit", "the minimum posterior probability for a node to be annotated"),
                        new Arguments.IntegerOption("limitCount", "the minimum sample count for a node to be annotated (default 5)"),
                        new Arguments.StringOption("target", "target_file_name", "specifies a user target tree to be annotated"),
                        new Arguments.StringOption("reference", "tree_file_name", "specifies a reference tree for sampled trees to be compared with"),
                        new Arguments.StringOption("metrics", "output_file_name", "file name to write tree metrics for each tree compared to the target"),
                        new Arguments.IntegerOption("threads", "max number of threads (default automatic)"),
                        new Arguments.Option("help", "option to print this message"),
                        new Arguments.Option("forceDiscrete", "forces integer traits to be treated as discrete traits."),
                        new Arguments.StringOption("hpd2D", "the HPD interval to be used for the bivariate traits", "specifies a (vector of comma separated) HPD proportion(s)"),
                        new Arguments.Option("ess", "compute ess for branch parameters")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            progressStream.println(ae);
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("forceDiscrete")) {
            progressStream.println("  Forcing integer traits to be treated as discrete traits.");
            forceIntegerToDiscrete = true;
        }

        if (arguments.hasOption("help")) {
            printUsage(arguments);
            System.exit(0);
        }

        HeightsSummary heights = DEFAULT_HEIGHTS_SUMMARY;
        if (arguments.hasOption("heights")) {
            String value = arguments.getStringOption("heights");
            if (value.equalsIgnoreCase("mean")) {
                heights = HeightsSummary.MEAN_HEIGHTS;
            } else if (value.equalsIgnoreCase("median")) {
                heights = HeightsSummary.MEDIAN_HEIGHTS;
            } else if (value.equalsIgnoreCase("ca")) {
                progressStream.println("CA heights are not supported - this has been superseded by the HIPSTR tree (--type hipstr)");
                printUsage(arguments);
                System.exit(1);
            }
        }

        long burninStates = -1;
        int burninTrees = -1;
        if (arguments.hasOption("burnin")) {
            burninStates = arguments.getLongOption("burnin");
        }
        if (arguments.hasOption("burninTrees")) {
            burninTrees = arguments.getIntegerOption("burninTrees");
        }

        if (arguments.hasOption("ess")) {
            if (burninStates != -1) {
                progressStream.println(" Calculating ESS for branch parameters.");
                computeESS = true;
            } else {
                throw new RuntimeException("Specify burnin as states to use 'ess' option.");
            }
        }

        double posteriorLimit = 0.0;
        if (arguments.hasOption("limit")) {
            posteriorLimit = arguments.getRealOption("limit");
        }

        if (arguments.hasOption("limitFrequency")) {
            posteriorLimit = arguments.getRealOption("limitFrequency");
        }

        int countLimit = 5;
        if (arguments.hasOption("limitCount")) {
            countLimit = arguments.getIntegerOption("limitCount");
        }

        double[] hpd2D = {80};
        if (arguments.hasOption("hpd2D")) {
            try {
                hpd2D = parseVariableLengthDoubleArray(arguments.getStringOption("hpd2D"));
            } catch (Arguments.ArgumentException e) {
                System.err.println("Error reading " + arguments.getStringOption("hpd2D"));
            }
        }

        Target target = Target.HIPSTR;
        if (arguments.hasOption("type")) {
            if (arguments.getStringOption("type").equalsIgnoreCase("MRHIPSTR")) {
                target = Target.MRHIPSTR;
            } else if (arguments.getStringOption("type").equalsIgnoreCase("MCC")) {
                target = Target.MAX_CLADE_CREDIBILITY;
            } else if (arguments.getStringOption("type").equalsIgnoreCase("MRC")) {
                target = Target.MAJORITY_RULE;
            }
        }

        if (arguments.hasOption("target")) {
            target = Target.USER_TARGET_TREE;
            targetTreeFileName = arguments.getStringOption("target");
        }

        if (arguments.hasOption("reference")) {
            referenceTreeFileName = arguments.getStringOption("reference");
        }

        if (arguments.hasOption("metrics")) {
            treeMetricFileName = arguments.getStringOption("metrics");
        }

        int threadCount = -1;
        if (arguments.hasOption("threads")) {
            threadCount = arguments.getIntegerOption("threads");
        }

        final String[] args2 = arguments.getLeftoverArguments();

        switch (args2.length) {
            case 2:
                outputFileName = args2[1];
                // fall to
            case 1:
                inputFileName = args2[0];
                break;
            default: {
                System.err.println("Unknown option: " + args2[2]);
                System.err.println();
                printUsage(arguments);
                System.exit(1);
            }
        }

        new TreeAnnotator(
                burninTrees,
                burninStates,
                heights,
                posteriorLimit,
                countLimit,
                hpd2D,
                computeESS,
                threadCount,
                target,
                targetTreeFileName,
                referenceTreeFileName,
                treeMetricFileName,
                inputFileName,
                outputFileName);

        if (target == Target.MAX_CLADE_CREDIBILITY) {
            progressStream.println("Found Maximum Clade Credibility (MCC) tree - citation: " +
                    "Drummond and Rambaut: 'BEAST: Bayesian evolutionary analysis by sampling trees', BMC Ecology and Evolution 2007, 7: 214.");
        } else if (target == Target.HIPSTR) {
            progressStream.println("Constructed Highest Independent Posterior Sub-Tree Reconstruction (HIPSTR) tree - citation: In prep.");
        } else if (target == Target.MRHIPSTR) {
            progressStream.println("Constructed Majority Rule Highest Independent Posterior Sub-Tree Reconstruction (MrHIPSTR) tree - citation: In prep.");
        } else if (target == Target.MAJORITY_RULE) {
            progressStream.println("Constructed majority-rule consensus tree");
        } else if (target == Target.USER_TARGET_TREE) {
//            progressStream.println("Loaded user target tree.");
        } else {
            throw new IllegalArgumentException("Unknown target option: " + target);
        }

        System.exit(0);
    }


}

