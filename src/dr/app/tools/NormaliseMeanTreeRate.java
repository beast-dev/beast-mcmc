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

import java.io.*;
import java.util.Iterator;
import java.util.Set;
import java.util.HashSet;

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

        //BufferedReader br = new BufferedReader(new FileReader(parentFile));
        //String line;
        //while ((line = br.readLine()) != null) {

        //}
        //br.close();

        FileReader fileReader = new FileReader(parentFile);
        TreeImporter importer = new NexusImporter(fileReader);
        try {
            int totalTrees = 0;
            Set<String> attributeNames = new HashSet<String>();
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();
                for (int i = 0; i < tree.getNodeCount(); i++) {
                    NodeRef node = tree.getNode(i);
                    System.out.print(tree.getNodeCount() + " Punk bitch " + tree.getNodeRate(node) + "\t");

                    Iterator iter = tree.getNodeAttributeNames(node);
                    if (iter != null) {
                        while (iter.hasNext()) {
                            String name = (String) iter.next();
                            System.out.println(name);
                            attributeNames.add(name);
                        }
                    }
                }
                System.out.println();
                totalTrees++;
            }
            //System.out.println("Total number of trees: " + totalTrees);

        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }






        analyze(parentFile, outFile, normaliseMeanRateTo);
    }

    /**
     * Recursively analyzes log files.
     *
     * @param inFile              input file
     * @param outFile             output file
     * @param normaliseMeanRateTo rate to normalise to
     * @throws dr.inference.trace.TraceException
     *          if the trace file is in the wrong format or corrupted
     */
    private static void analyze(File inFile, File outFile, double normaliseMeanRateTo) throws TraceException {

        /*try {




        }catch (IOException e) {
            e.printStackTrace();
        }*/
















        /*if (file.isFile()) {
            try {

                String name = file.getCanonicalPath();
                if (verbose) {
                    TraceAnalysis.report(name, burnin, marginalLikelihood);
                } else {
                    TraceAnalysis.shortReport(name, burnin, drawHeader[0], hpds, ess, stdErr, marginalLikelihood);
                    drawHeader[0] = false;
                }
            } catch (IOException e) {
                //e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        } else {
            File[] files = file.listFiles();
            for (File f : files) {
                if (f.isDirectory()) {
                    analyze(f, burnin, verbose, drawHeader, hpds, ess, stdErr, marginalLikelihood);
                } else if (f.getName().endsWith(".log") || f.getName().endsWith(".p")) {
                    analyze(f, burnin, verbose, drawHeader, hpds, ess, stdErr, marginalLikelihood);
                } else {
                    if (verbose) System.out.println("Ignoring file: " + f);
                }
            }
        }*/
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
