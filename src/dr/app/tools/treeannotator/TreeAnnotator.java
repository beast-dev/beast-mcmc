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
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.util.Version;
import jam.console.ConsoleApplication;

import javax.swing.*;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

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

    // Messages to stderr, output to stdout
    private static PrintStream progressStream = System.err;

    private final CollectionAction collectionAction;
    private final AnnotationAction annotationAction;

    private TaxonList taxa = null;
    private int totalTrees;
    private int totalTreesUsed;
    private long maxState;

    enum Target {
        MAX_CLADE_CREDIBILITY("Maximum clade credibility tree"),
        HIPSTR("Highest independent posterior subtree reconstruction (HIPSTR)"),
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
        MEDIAN_HEIGHTS("Median heights"),
        MEAN_HEIGHTS("Mean heights"),
        KEEP_HEIGHTS("Keep target heights"),
        CA_HEIGHTS("Common Ancestor heights");

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
                         HeightsSummary heightsOption,
                         double posteriorLimit,
                         double[] hpd2D,
                         boolean computeESS,
                         Target targetOption,
                         String targetTreeFileName,
                         String inputFileName,
                         String outputFileName
    ) throws IOException {

        long totalStartTime = System.currentTimeMillis();

        collectionAction = new CollectionAction();

        collectionAction.addAttributeName("height");
        collectionAction.addAttributeName("length");

        annotationAction = new AnnotationAction(heightsOption, posteriorLimit, hpd2D, computeESS, true);

        annotationAction.addAttributeName("height");
        annotationAction.addAttributeName("length");

        int burnin = -1;

        totalTrees = 10000;
        totalTreesUsed = 0;

        CladeSystem cladeSystem = new CladeSystem(targetOption == Target.HIPSTR);

        // read the clades in even if a target tree so it can have its stats reported
//        if (targetOption != Target.USER_TARGET_TREE) {
            // if we are not just annotating a specific target tree
            // then we need to read all the trees into a CladeSystem
            // to get Clade and SubTree frequencies.
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

            progressStream.println("Total unique clades: " + cladeSystem.getCladeCount());
            progressStream.println();
//        }

        MutableTree targetTree = null;

        switch (targetOption) {
            case USER_TARGET_TREE: {
                if (targetTreeFileName != null) {
                    targetTree = readUserTargetTree(targetTreeFileName, targetTree, cladeSystem);
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
                targetTree = getHIPSTRTree(cladeSystem);
                break;
            }
            default: throw new IllegalArgumentException("Unknown targetOption");
        }

        // Help garbage collector
        cladeSystem = null;

        CladeSystem targetCladeSystem = new CladeSystem(targetTree);

        collectNodeAttributes(targetCladeSystem, inputFileName, burnin);

        annotateTargetTree(targetCladeSystem, heightsOption, targetTree);

        if (heightsOption == HeightsSummary.CA_HEIGHTS) {
            setNodeHeightsCA(targetCladeSystem, targetTree);
        }

        writeAnnotatedTree(outputFileName, targetTree);

        long timeElapsed =  (System.currentTimeMillis() - totalStartTime) / 1000;
        progressStream.println("Total time: " + timeElapsed + " secs");
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

            totalTrees = 0;
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();
                long state = 0;

                if (taxa == null) {
                    taxa = new Taxa(tree);
                }

                if (burninStates > 0) {
                    // if burnin has been specified in states, try to parse it out...
                    String name = tree.getId().trim();

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

                    cladeSystem.add(tree);

                    totalTreesUsed += 1;
                }

                if (totalTrees > 0 && totalTrees % stepSize == 0) {
                    progressStream.print("*");
                    progressStream.flush();
                }
                totalTrees++;
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
        TreeImporter importer = new NexusImporter(reader, true);

        long startTime = System.currentTimeMillis();

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

                    cladeSystem.traverseTree(tree, collectionAction);
                    totalTreesUsed += 1;
                }
                if (counter > 0 && counter % stepSize == 0) {
                    progressStream.print("*");
                    progressStream.flush();
                }
                counter++;

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
        Set<String> attributeNames = new HashSet<>();
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
        // plugins now live in Annotation action and I am not sure this is need...
//        for (AnnotationAction.TreeAnnotationPlugin plugin : plugins) {
//            Set<String> claimed = plugin.setAttributeNames(attributeNames);
//            attributeNames.removeAll(claimed);
//        }
        collectionAction.addAttributeNames(attributeNames);
        annotationAction.addAttributeNames(attributeNames);
    }

    private MutableTree readUserTargetTree(String targetTreeFileName, MutableTree targetTree, CladeSystem cladeSystem) throws IOException {
        progressStream.println("Reading user specified target tree, " + targetTreeFileName + ", ...");

        NexusImporter importer = new NexusImporter(new FileReader(targetTreeFileName));
        try {
            Tree tree = importer.importNextTree();
            if (tree == null) {
                NewickImporter x = new NewickImporter(new FileReader(targetTreeFileName));
                tree = x.importNextTree();
            }
            if (tree == null) {
                System.err.println("No tree in target nexus or newick file " + targetTreeFileName);
                System.exit(1);
            }
            targetTree = new FlexibleTree(tree);
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Target Tree: " + e.getMessage());
            System.exit(1);
        }

        progressStream.println();
        double score = scoreTree(targetTree, cladeSystem);
        progressStream.println("Target tree's log clade credibility: " + String.format("%.4f", score));
        reportStatistics(cladeSystem, targetTree);

        progressStream.println();
        return targetTree;
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
        progressStream.println();

        return bestTree;
    }

    private MutableTree getHIPSTRTree(CladeSystem cladeSystem) {

        long startTime = System.currentTimeMillis();

        HIPSTRTreeBuilder treeBuilder = new HIPSTRTreeBuilder();
        MutableTree tree = treeBuilder.getHIPSTRTree(cladeSystem, taxa);
        double score = treeBuilder.getScore();

        long timeElapsed =  (System.currentTimeMillis() - startTime) / 1000;
        progressStream.println("[" + timeElapsed + " secs]");
        progressStream.println();
        progressStream.println("HIPSTR tree's log clade credibility: " + String.format("%.4f", score));
        reportStatistics(cladeSystem, tree);
        progressStream.println();

        return tree;
    }

    private static void reportStatistics(CladeSystem cladeSystem, Tree tree) {
        progressStream.println("Lowest individual clade credibility: " + String.format("%.4f", cladeSystem.getMinimumCladeCredibility(tree)));
        progressStream.println("Mean individual clade credibility: " + String.format("%.4f", cladeSystem.getMeanCladeCredibility(tree)));
        progressStream.println("Median individual clade credibility: " + String.format("%.4f", cladeSystem.getMedianCladeCredibility(tree)));
        progressStream.println("Number of clades with credibility 1.0: " + cladeSystem.getTopCladeCredibility(tree, 1.0));
        progressStream.println("Number of clades with credibility > 0.99: " + cladeSystem.getTopCladeCredibility(tree, 0.99) +
                " (out of " + cladeSystem.getTopCladeCredibility(0.99) + " in all trees)");
        progressStream.println("Number of clades with credibility > 0.95: " + cladeSystem.getTopCladeCredibility(tree, 0.95) +
                " (out of " + cladeSystem.getTopCladeCredibility(0.95) + " in all trees)");
        progressStream.println("Number of clades with credibility > 0.5: " + cladeSystem.getTopCladeCredibility(tree, 0.5) +
        " (out of " + cladeSystem.getTopCladeCredibility(0.5) + " in all trees)");
    }

    private void annotateTargetTree(CladeSystem cladeSystem, HeightsSummary heightsOption, MutableTree targetTree) {
        progressStream.println("Annotating target tree...");

        try {
            cladeSystem.traverseTree(targetTree, annotationAction);
        } catch (Exception e) {
            System.err.println("Error annotating tree: " + e.getMessage() + "\nPlease check the tree log file format.");
            System.exit(1);
        }
        progressStream.println();
    }

    private void setNodeHeightsCA(CladeSystem cladeSystem, MutableTree targetTree) {
        assert false : "Implement this";

        long startTime = System.currentTimeMillis();

        progressStream.println("Setting node heights...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

//        CAHeights caHeights = new CAHeights(this);
//        caHeights.setTreeHeightsByCA(targetTree, inputFileName, burnin);

        long timeElapsed = (System.currentTimeMillis() - startTime) / 1000;
        progressStream.println("* [" + timeElapsed + " secs]");
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
        progressStream.println("  Example: treeannotator test.trees out.txt");
        progressStream.println("  Example: treeannotator -burnin 100 -heights mean test.trees out.txt");
        progressStream.println("  Example: treeannotator -burnin 100 -target map.tree test.trees out.txt");
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
            String nameString = "TreeAnnotator " + versionString;
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

            if (!dialog.showDialog("TreeAnnotator " + versionString)) {
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
                        hpd2D,
                        computeESS,
                        targetOption,
                        targetTreeFileName,
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
                        new Arguments.StringOption("type", new String[]{"mcc", "hipstr"}, false, "an option of 'mcc' or 'hipstr'"),
                        new Arguments.StringOption("heights", new String[]{"keep", "median", "mean", "ca"}, false,
                                "an option of 'keep', 'median' or 'mean' (default)"),
                        //"an option of 'keep', 'median', 'mean' or 'ca' (default)"),
                        new Arguments.LongOption("burnin", "the number of states to be considered as 'burn-in'"),
                        new Arguments.IntegerOption("burninTrees", "the number of trees to be considered as 'burn-in'"),
                        new Arguments.RealOption("limit", "the minimum posterior probability for a node to be annotated"),
                        new Arguments.StringOption("target", "target_file_name", "specifies a user target tree to be annotated"),
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
            System.out.println("  Forcing integer traits to be treated as discrete traits.");
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
                progressStream.println("CA heights are not supported - to avoid negative branch lengths, construct a HIPSTR tree");
                printUsage(arguments);
                System.exit(1);
//                heights = HeightsSummary.CA_HEIGHTS;
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
                System.out.println(" Calculating ESS for branch parameters.");
                computeESS = true;
            } else {
                throw new RuntimeException("Specify burnin as states to use 'ess' option.");
            }
        }

        double posteriorLimit = 0.0;
        if (arguments.hasOption("limit")) {
            posteriorLimit = arguments.getRealOption("limit");
        }

        double[] hpd2D = {80};
        if (arguments.hasOption("hpd2D")) {
            try {
                hpd2D = parseVariableLengthDoubleArray(arguments.getStringOption("hpd2D"));
            } catch (Arguments.ArgumentException e) {
                System.err.println("Error reading " + arguments.getStringOption("hpd2D"));
            }
        }

        Target target = Target.MAX_CLADE_CREDIBILITY;
        if (arguments.hasOption("type") && arguments.getStringOption("type").equalsIgnoreCase("HIPSTR")) {
            target = Target.HIPSTR;
        }

        if (arguments.hasOption("target")) {
            target = Target.USER_TARGET_TREE;
            targetTreeFileName = arguments.getStringOption("target");
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

        new TreeAnnotator(burninTrees, burninStates, heights, posteriorLimit, hpd2D, computeESS, target, targetTreeFileName, inputFileName, outputFileName);

        if (target == Target.MAX_CLADE_CREDIBILITY) {
            progressStream.println("Found Maximum Clade Credibility (MCC) tree - citation: " +
                    "Drummond and Rambaut: 'BEAST: Bayesian evolutionary analysis by sampling trees', BMC Ecology and Evolution 2007, 7: 214.");
        } else if (target == Target.HIPSTR) {
            progressStream.println("Constructed Highest Independent Posterior Sub-Tree Reconstruction (HIPSTR) tree - citation: In prep.");
        } else if (target == Target.USER_TARGET_TREE) {
            progressStream.println("Loaded user target tree.");
        }

        if (heights == HeightsSummary.CA_HEIGHTS) {
            progressStream.println("\nUsed Clade Height option - citation: " +
                    "Heled and Bouckaert: 'Looking for trees in the forest: " +
                    "summary tree from posterior samples'. BMC Evolutionary Biology 2013 13:221.");
        }

        System.exit(0);
    }


}

