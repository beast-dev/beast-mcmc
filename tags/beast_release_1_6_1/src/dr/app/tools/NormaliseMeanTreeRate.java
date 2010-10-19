package dr.app.tools;

import dr.util.Version;
import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.inference.trace.TraceException;
import dr.inference.model.Parameter;
import dr.evolution.io.TreeImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.Importer;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.FlexibleTree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.NodeTraitLogger;

import java.io.*;
import java.util.*;

/**
 *
 * Date: 27/08/2009
 * Time: 16:43:14
 *
 * Standalone application, that also contains public methods, for normalising the rates (and the times) of trees
 * to a have a mean rate value that is specified. This is required when using discretized branch rates. Setting
 * normaliseMeanRateTo in logTree will trigger the normalisation whenever trees are logged  
 *
 * @author Wai Lok Sibon Li
 *
 */
public class NormaliseMeanTreeRate {

    private final static Version version = new BeastVersion();


    public NormaliseMeanTreeRate(String inputFileName, String outputFileName, double normaliseMeanRateTo) throws java.io.IOException {

        File parentFile = new File(inputFileName);

        if (parentFile.isFile()) {
            System.out.println("Analysing tree file: " + inputFileName);
        } else {
            System.err.println("File " + inputFileName + " does not exist!");
            System.exit(0);
        }

        File outFile = new File(outputFileName);
        if (outputFileName != null) {
            FileOutputStream outputStream = new FileOutputStream(outFile);
            System.setOut(new PrintStream(outputStream));
        }

        if(!outFile.canWrite()) {
            System.err.println("Cannot write to file" + outFile.getAbsolutePath());
            System.exit(0);
        }

        FileReader fileReader = new FileReader(parentFile);
        TreeImporter importer = new NexusImporter(fileReader);
        ArrayList<Tree> treeList = new ArrayList<Tree>();
        ArrayList<String> treeNames = new ArrayList<String>();
        try {
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();
                analyze(tree, normaliseMeanRateTo);
                treeList.add(tree);
                treeNames.add(tree.getId());
            }

            NexusExporter exporter = new NexusExporter(System.out);
            exporter.setSortedTranslationTable(true);
            exporter.exportTrees(treeList.toArray(new Tree[treeList.size()]),
                true, treeNames.toArray(new String[treeNames.size()]));

        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }
    }

    /**
     * Normalises individual trees to the mean rate
     *
     * @param tree                tree to normalise
     * @param normaliseMeanRateTo rate to normalise to
     *          if the trace file is in the wrong format or corrupted
     */
    public static void analyze(Tree tree, double normaliseMeanRateTo) {
        double treeRate = 0;
        double treeTime = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);

            if(tree.getClass().getName().equals("dr.evomodel.tree.TreeModel")) {

                throw new RuntimeException("Does not currently handle TreeModel");
            }

            if(!tree.isRoot(node)) {

                if(tree.getNodeAttribute(node, "rate") == null) {
                    System.out.println("Tree file does not contain rate information. ");
                    System.setOut(System.out);
                    System.err.println("Tree file does not contain rate information. Program terminated");
                    System.exit(0);
                }
                treeRate += tree.getNodeRate(node) * tree.getBranchLength(node);
                treeTime += tree.getBranchLength(node);
            }
        }

        treeRate /= treeTime;

        /* Normalise the rates here */
        FlexibleTree modifiedTree = (FlexibleTree) tree;
        for (int i = 0; i < modifiedTree.getNodeCount(); i++) {
            NodeRef node = modifiedTree.getNode(i);
            if(!modifiedTree.isRoot(node)) {
                double nodeRate = normaliseMeanRateTo * modifiedTree.getNodeRate(node) / treeRate;
                modifiedTree.setNodeAttribute(node, "rate", Double.valueOf(nodeRate));
                double nodeTime = modifiedTree.getBranchLength(node);
                nodeTime = nodeTime * treeRate / normaliseMeanRateTo;
                modifiedTree.setBranchLength(node, nodeTime);

            }
        }
    }

    public static void printTitle() {
        System.out.println();
        centreLine("NormaliseMeanTreeRate " + version.getVersionString() + ", " + version.getDateString(), 60);
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

        arguments.printUsage("normaliseMeanTreeRate", "[-input-file-name <input-file-name>] [-output-file-name <output-file-name>] [-normaliseMeanRateTo <normaliseMeanRateTo>]");
        System.out.println();
        System.out.println("  Example: normaliseMeanTreeRate test.trees out.trees 1.0");
        System.out.println();

    }

    //Main method
    public static void main(String[] args) throws java.io.IOException, TraceException {

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.StringOption("input-file-name", "infile", "Input file name"),
                        new Arguments.StringOption("output-file-name", "outfile", "Output file name"),
                        new Arguments.RealOption("normaliseMeanRateTo", "Mean rate we should normalise to"),
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

        String inputFileName = null;
        if(arguments.hasOption("input-file-name")) {
            inputFileName = arguments.getStringOption("input-file-name");
        }

        String outputFileName = null;
        if(arguments.hasOption("output-file-name")) {
            outputFileName = arguments.getStringOption("output-file-name");
        }

        double normaliseMeanRateTo = Double.NaN;
        if(arguments.hasOption("normaliseMeanRateTo")) {
            normaliseMeanRateTo = arguments.getRealOption("normaliseMeanRateTo");
        }

        if (inputFileName == null) {
            // No input file name was given so throw up a dialog box...
            inputFileName = Utils.getLoadFileName("NormaliseMeanTreeRate " + version.getVersionString() + " - Select log file to analyse");
        }

        if (outputFileName == null) {
            // No input file name was given so throw up a dialog box...
            outputFileName = Utils.getSaveFileName("NormaliseMeanTreeRate " + version.getVersionString() + " - Select file to save to");
        }

        if(Double.isNaN(normaliseMeanRateTo)) {
            System.out.println("Enter rate value to normalise to: ");
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            normaliseMeanRateTo = Double.parseDouble(br.readLine());
        }


        new NormaliseMeanTreeRate(inputFileName, outputFileName, normaliseMeanRateTo);

        System.out.println("Please bear in mind that the trees files are unchanged and results may vary slightly from if you ran them internally");

        System.exit(0);
    }

}
