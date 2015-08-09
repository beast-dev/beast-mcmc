/*
 * NodeStateAnalyser.java
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

package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.util.Version;
import dr.evolution.io.TreeImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.Importer;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;

import java.io.*;
import java.util.*;


public class NodeStateAnalyser {

    private final static Version version = new BeastVersion();


    public NodeStateAnalyser(int burnin, double mrsd, double scale, String inputFileName, String outputFileName) throws IOException {

        File parentFile = new File(inputFileName);

        if (parentFile.isDirectory()) {
            System.out.println("Analysing all trees files below directory: " + inputFileName);
        } else if (parentFile.isFile()) {
            System.out.println("Analysing tree file: " + inputFileName);
        } else {
            System.err.println(inputFileName + " does not exist!");
            System.exit(0);
        }

        if (outputFileName != null) {
            FileOutputStream outputStream = new FileOutputStream(outputFileName);
            System.setOut(new PrintStream(outputStream));
        }

        analyze(parentFile, burnin, mrsd, scale);
    }

    /**
     * Recursively analyzes log files.
     *
     * @param file       the file to analyze (if this is a directory then the files within it are analyzed)
     * @param burnin     the burnin to use
     * @throws dr.inference.trace.TraceException
     *          if the trace file is in the wrong format or corrupted
     */
    private void analyze(File file, int burnin, double mrsd, double scale) {

        if (file.isFile()) {
            try {

                String name = file.getCanonicalPath();

                report(name, burnin, mrsd, scale);
                //TraceAnalysis.report(name, burnin, marginalLikelihood);

            } catch (IOException e) {
                //e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        } else {
            File[] files = file.listFiles();
            for (File f : files) {
                if (f.isDirectory()) {
                    analyze(f, burnin, mrsd, scale);
                } else if (f.getName().endsWith(".trees")) {
                    analyze(f, burnin, mrsd, scale);
                }
            }
        }
    }


    /**
     * Recursively analyzes trees files.
     *
     * @param name       the file to analyze (if this is a directory then the files within it are analyzed)
     * @param burnin     the burnin to use
     */
    private void report(String name, int burnin, double mrsd, double scale) {
        int count = 0;

        Map<String, List<Double>> heightMap = new HashMap<String, List<Double>>();
        Set<String> states = new TreeSet<String>();

        try {
            FileReader fileReader = new FileReader(new File(name));
            TreeImporter importer = new NexusImporter(fileReader);
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();
                if(count>=burnin) {
                    for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                        NodeRef node = tree.getInternalNode(i);
                        Object value = tree.getNodeAttribute(node, "state");
                        if (value != null) {
                            String state = value.toString();
                            List<Double> heights = heightMap.get(state);
                            if (heights == null) {
                                heights = new ArrayList<Double>();
                                heightMap.put(state, heights);
                                states.add(state);
                            }
                            double h = tree.getNodeHeight(node) * scale;
                            if (Double.isNaN(mrsd)) {
                                heights.add(h);
                            } else {
                                heights.add(mrsd - h);
                            }
                        } else {
                            System.out.println("Node missing state");
                        }
                    }

                }
                count++;
            }

            boolean first = true;
            int maxCount = 0;
            for (String state : states) {
                if (first) {
                    first = false;
                } else {
                    System.out.print("\t");
                }
                System.out.print(state);

                List<Double> heights = heightMap.get(state);
                if (heights.size() > maxCount) {
                    maxCount = heights.size();
                }
            }
            System.out.println();

            for (int i = 0; i < maxCount; i++) {
                first = true;
                for (String state : states) {
                    if (first) {
                        first = false;
                    } else {
                        System.out.print("\t");
                    }

                    List<Double> heights = heightMap.get(state);
                    if (i < heights.size()) {
                        System.out.print(heights.get(i));
                    }
                }
                System.out.println();
            }

        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
        } catch (IOException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
        }
    }

    public static void printTitle() {
        System.out.println();
        centreLine("NodeStateAnalyser " + version.getVersionString() + ", " + version.getDateString(), 60);
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

        arguments.printUsage("nodestateanalyser", "[-burnin <burnin>][<input-file-name> [<output-file-name>]]");
        System.out.println();
        System.out.println("  Example: nodestateanalyser test.tree");
        System.out.println("  Example: nodestateanalyser -burnin 10000 trees.log out.txt");
        System.out.println();

    }

    //Main method
    public static void main(String[] args) throws IOException {

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in'"),
                        new Arguments.RealOption("mrsd","specifies the most recent sampling data in fractional years to rescale time [default=0]"),
                        new Arguments.RealOption("scale","Provide a scaling factor for the node heights [default=1]"),
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

        int burnin = -1;
        if (arguments.hasOption("burnin")) {
            burnin = arguments.getIntegerOption("burnin");
        }

        double scale = 1.0;
        if (arguments.hasOption("scale")) {
            scale = arguments.getRealOption("scale");
        }


        double mrsd = Double.NaN;
        if (arguments.hasOption("mrsd")) {
            mrsd = arguments.getRealOption("mrsd");
        }

        String inputFileName = null;
        String outputFileName = null;

        String[] args2 = arguments.getLeftoverArguments();

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
            inputFileName = Utils.getLoadFileName("NodeStateAnalyser " + version.getVersionString() + " - Select tree file to analyse");
        }

        if(burnin==-1) {
            System.out.println("Enter number of trees to burn-in (integer): ");
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            burnin = Integer.parseInt(br.readLine());
        }

        new NodeStateAnalyser(burnin, mrsd, scale, inputFileName, outputFileName);

        System.exit(0);
    }
}