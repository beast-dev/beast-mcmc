package dr.app.tools;

import dr.util.Version;
import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.inference.trace.TraceException;
import dr.inference.trace.TraceAnalysis;
import dr.evolution.io.TreeImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.Importer;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.FlexibleTree;

import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Sibon
 * Date: 27/08/2009
 * Time: 16:43:14
 * To change this template use File | Settings | File Templates.
 */
public class NormaliseMeanTreeRate {

    private final static Version version = new BeastVersion();


    public NormaliseMeanTreeRate(String inputFileName, String outputFileName, double normaliseMeanRateTo) throws java.io.IOException, TraceException {

        File parentFile = new File(inputFileName);

        if (parentFile.isFile()) {
            System.out.println("Analysing log file: " + inputFileName);
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
            int totalTrees = 0;
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();
                analyze(tree, normaliseMeanRateTo);
                treeList.add(tree);
                treeNames.add(tree.getId());
                totalTrees++;

            }
            new NexusExporter(System.out).exportTrees(treeList.toArray(new Tree[treeList.size()]),
                    true, treeNames.toArray(new String[treeNames.size()]));

        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }

        //analyze(parentFile, outFile, normaliseMeanRateTo);
    }

    /**
     * Recursively analyzes log files.
     *
     * @param tree                tree to normalise
     * @param normaliseMeanRateTo rate to normalise to
     * @throws dr.inference.trace.TraceException
     *          if the trace file is in the wrong format or corrupted
     */
    public static void analyze(Tree tree, double normaliseMeanRateTo) throws TraceException {
        double treeRate = 0;
        double treeTime = 0;
        //int branchCount = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);

            if(!tree.isRoot(node)) {

                if(tree.getNodeAttribute(node, "rate") == null) {
                    System.out.println("Tree file does not contain rate information. ");
                    System.setOut(System.out);
                    System.err.println("Tree file does not contain rate information. Program terminated");
                    System.exit(0);
                }
                //double nodeRate = (Double) tree.getNodeAttribute(node, "rate");
                treeRate += (Double) tree.getNodeAttribute(node, "rate") * tree.getBranchLength(node);
                treeTime += tree.getBranchLength(node);
//                System.out.println(tree.getNodeAttribute(node, "rate") + "\t" + tree.getBranchLength(node));
                //branchCount++;
            }

            /*Iterator iter = tree.getNodeAttributeNames(node);
            if (iter != null) {
                System.out.print(".");
                while (iter.hasNext()) {
                    //System.out.print("|");
                    String name = (String) iter.next();
                    //System.out.println(" fell off " + name);
                    double nodeRate = (Double) tree.getNodeAttribute(node, "rate");
                    //System.out.println("you the fucking best" + nodeRate + "dasf");
                    //if() {
                        treeRate += (Double) tree.getNodeAttribute(node, "rate");
                    //}
                    //attributeNames.add(name);
                }
            }*/
        }

        //treeRate /= branchCount;
        treeRate /= treeTime;
//        System.out.println(" 1 " + tree);

        /* Normalise the rates here */
        FlexibleTree modifiedTree = (FlexibleTree) tree;//.getCopy();
        for (int i = 0; i < modifiedTree.getNodeCount(); i++) {
            NodeRef node = modifiedTree.getNode(i);
            if(!modifiedTree.isRoot(node)) {
                double nodeRate = (Double) modifiedTree.getNodeAttribute(node, "rate");
                nodeRate = normaliseMeanRateTo * nodeRate / treeRate;

                modifiedTree.setNodeAttribute(node, "rate", Double.valueOf(nodeRate));

                //double newValue =(Double) modifiedTree.getNodeAttribute(node, "rate");
                //System.out.println(nodeRate + "\t" + newValue);


                //double nodeTime = (Double) modifiedTree.getNodeAttribute(node, "t");
                double nodeTime = modifiedTree.getBranchLength(node);
                nodeTime = nodeTime * treeRate / normaliseMeanRateTo;
                modifiedTree.setBranchLength(node, nodeTime);

//                System.out.println(nodeRate + "\t" + modifiedTree.getBranchLength(node));
                //modifiedTree.setNodeAttribute(node, "t", Double.valueOf(nodeTime));

            }
        }
        //tree = modifiedTree;


//        System.out.println(" 2 " + tree);

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

        arguments.printUsage("normaliseMeanTreeRate", "[-input-file-name <input-file-name>] [<output-file-name> <output-file-name>] [<normaliseMeanRateTo> <normaliseMeanRateTo>]");
        System.out.println();
        System.out.println("  Example: normaliseMeanTreeRate test.trees out.trees 1.0");
        System.out.println("  Example: loganalyser -input-file-name test.trees out.trees 1.0");
        System.out.println();

    }

    //Main method
    public static void main(String[] args) throws java.io.IOException, TraceException {

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("input-file-name", "Input file name"),
                        new Arguments.IntegerOption("output-file-name", "Output file name"),
                        new Arguments.Option("normaliseMeanRateTo", "Mean rate we should normalise to"),
//				new Arguments.Option("html", "format output as html"),
//				new Arguments.Option("svg", "generate svg graphics"),
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

        /*int burnin = -1;
        if (arguments.hasOption("burnin")) {
            burnin = arguments.getIntegerOption("burnin");
        }

        boolean hpds = arguments.hasOption("hpd");
        boolean ess = arguments.hasOption("ess");
        boolean stdErr = arguments.hasOption("stdErr");
        boolean shortReport = arguments.hasOption("short");

        String marginalLikelihood = null;
        if (arguments.hasOption("marginal")) {
            marginalLikelihood = arguments.getStringOption("marginal");
        }*/

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

        //String outputFileName = arguments.getStringOption("output-file-name");
        //double normaliseMeanRateTo = arguments.getRealOption("normaliseMeanRateTo");

        /*String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 2) {
            System.err.println("Unknown option: " + args2[2]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (args2.length > 0) {
            inputFileName = args2[0];
        }
        if (args2.length > 1) {
            outputFileName = args2[1];
        }

        if (inputFileName == null) {
            // No input file name was given so throw up a dialog box...
            inputFileName = Utils.getLoadFileName("NormaliseMeanTreeRate " + version.getVersionString() + " - Select tree file to normalise");
        }*/

        if (inputFileName == null) {
            // No input file name was given so throw up a dialog box...
            inputFileName = Utils.getLoadFileName("NormaliseMeanTreeRate " + version.getVersionString() + " - Select log file to analyse");
        }

        if (outputFileName == null) {
            // No input file name was given so throw up a dialog box...
            outputFileName = Utils.getSaveFileName("NormaliseMeanTreeRate " + version.getVersionString() + " - Select log file to analyse");
        }

        if(Double.isNaN(normaliseMeanRateTo)) {
            System.out.println("Enter rate value to normalise to: ");
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            normaliseMeanRateTo = Double.parseDouble(br.readLine());
        }


        new NormaliseMeanTreeRate(inputFileName, outputFileName, normaliseMeanRateTo);

        System.exit(0);
    }

}
