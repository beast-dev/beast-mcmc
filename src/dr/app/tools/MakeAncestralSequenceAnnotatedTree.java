/*
 * MakeAncestralSequenceAnnotatedTree.java
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
import dr.inference.trace.TraceException;
import dr.util.Version;

import java.io.*;
import java.util.Hashtable;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.inference.trace.TraceException;
import dr.util.Version;

import java.io.*;
import java.util.*;

/**
 * MakeAncestralSequenceAnnotatedTree - Converts BAli-Phy *.out, *.fastas, *.trees files into
 * NEXUS format readable by BEAST and AncestralSequenceAnnotator
 * Based on Marc A. Suchard's make-nexus-file script
 * Reimplemented by Wai Lok Sibon Li in BEAST for compatibility with BAli-Phy v2.1.0 and above
 *
 *
 * @author Wai Lok Sibon Li
 */

public class MakeAncestralSequenceAnnotatedTree {

    private final static Version version = new BeastVersion();


    public MakeAncestralSequenceAnnotatedTree(int thin, int skip, String outInputFileName, String fastasInputFileName,
                                              String treesInputfileName, String outputFileName) throws IOException, TraceException {

        File outInputFile = new File(outInputFileName);

        if (outInputFile.isFile()) {
            System.out.println("Analysing *.out file: " + outInputFileName);
        } else {
            System.err.println(outInputFileName + " does not exist!");
            System.exit(0);
        }

        File fastasInputFile = new File(fastasInputFileName);

        if (fastasInputFile.isFile()) {
            System.out.println("Analysing *.fastas file: " + fastasInputFileName);
        } else {
            System.err.println(fastasInputFileName + " does not exist!");
            System.exit(0);
        }

        File treesInputfile = new File(treesInputfileName);

        if (treesInputfile.isFile()) {
            System.out.println("Analysing *.trees file: " + treesInputfileName);
        } else {
            System.err.println(treesInputfileName + " does not exist!");
            System.exit(0);
        }


        if (outputFileName != null) {
            FileOutputStream outputStream = new FileOutputStream(outputFileName);
            System.setOut(new PrintStream(outputStream));
        }

        analyze(outInputFile, fastasInputFile, treesInputfile, skip, thin);
    }

    /**
     *
     * Analyses the output files
     *
     */
    private void analyze(File outInputFile, File fastasInputFile, File treesInputfileName, int skip, int thin) throws TraceException {

        try {
            BufferedReader outFileReader = new BufferedReader(new FileReader(outInputFile));
            BufferedReader fastasFileReader = new BufferedReader(new FileReader(fastasInputFile));
            BufferedReader treesFileReader = new BufferedReader(new FileReader(treesInputfileName));
            //BufferedReader br = new BufferedReader(fileReader);

            String line;

            int count = 0;
            int printedCount = 0;
            String substModel = null;

            ArrayList<String> trees = new ArrayList<String>();
            Hashtable<String, Boolean> translateHash = new Hashtable<String, Boolean>();
            Hashtable<String, Boolean> taxaHash = new Hashtable<String, Boolean>();




            //int discreteVariableIndex = -1;
            //Hashtable<String, Integer> hash = new Hashtable<String, Integer>();
            line = outFileReader.readLine();
            while(line != null) {



                if(line.matches("subst model\\d* = .*?")) {
                    substModel = line.replaceFirst("subst model\\d* = ", "").replaceFirst("\\(.*\\/(.*)\\.dat\\)/\\(","").replaceFirst("\\.dat\\)", "").replaceAll(" \\+ ", "+");
                    line = outFileReader.readLine();
                }
                else if(line.matches("iterations = \\d+")) {
                    int iterations = Integer.parseInt(line.replaceFirst("iterations = ",""));
                    if(iterations % 10 ==0) { //todo Alignment is only printed every 10 iterations. Figure out a proper way to manage this
                        // parse the model here
                        double likelihood = 0;
                        String model = null;
                        search: while((line=outFileReader.readLine()) != null) {

                            if(line.matches("iterations = \\d+")) {
                                break search;
                            }
                            else if (line.matches(".*?likelihood = -\\d+.*?")) {
                                likelihood = Double.parseDouble(line.replaceFirst(".+likelihood = ", "")
                                        .replaceFirst("\\s.+", ""));
                            }
                            else if (line.matches("\\s+Heat:beta = .+")) {
                                model = line.replaceAll(" = ", "=").replaceAll("::", "_").replaceAll("\\*", "")
                                        .replaceAll("^ ", "").replaceAll("\\s+", ",");
                            }

                        }



                        String tree = treesFileReader.readLine();

                        // parse the sequence here


                        String line2;
                        int previousLineLength=1;
                        search: while((line2=fastasFileReader.readLine()) != null) {
                            if(line2.matches("iterations = \\d+")) {
                                int fastaIterations = Integer.parseInt(line2.replaceFirst("iterations = ",""));
                                if(fastaIterations!=iterations) {
                                    throw new RuntimeException("Iterations is not the same: " +
                                            fastaIterations + ", " + iterations);
                                }
                            }
                            else if(line2.length()==0 && previousLineLength==0) {
                                break search;
                            }
                            else if(line2.matches(">.+")) {

                                String sequence = fastasFileReader.readLine();
                                String sequenceName = line2.replaceFirst(">","").replaceFirst("\\s+.+","");
                                if(tree.matches(".*\\(\\s*"+sequenceName+":.+")) {
                                    taxaHash.put(sequenceName, false);
                                    tree = tree.replaceFirst("\\(\\s*"+sequenceName+":", "("+sequenceName+"[&states="+sequence+"]:");
                                }
                                else if(tree.matches(".*,\\s*"+sequenceName+":.+")) {
                                    taxaHash.put(sequenceName, false);
                                    tree = tree.replaceFirst(",\\s*"+sequenceName+":", ","+sequenceName+"[&states="+sequence+"]:");
                                }
                                else if(tree.matches(".*\\)\\s*"+sequenceName+":.+")) { // This one is for internal nodes
                                    tree = tree.replaceFirst("\\)\\s*"+sequenceName+":", ")"+sequenceName+"[&states="+sequence+"]:");
                                }
                                else if(tree.matches(".*\\)\\s*"+sequenceName+";.*")) { // This one is for internal nodes (root node)
                                    tree = tree.replaceFirst("\\)\\s*"+sequenceName+";", ")"+sequenceName+"[&states="+sequence+"];");
                                }
                                else {
                                    throw new RuntimeException("Sequence name " + sequenceName + " does not exist in tree " + tree);
                                }

                                translateHash.put(sequenceName, false);
                            }
                            previousLineLength = line2.length();
                        }



                        if((count % thin)==0 && count >= skip) {
                            trees.add("tree STATE_"+iterations+" [&lnL="+likelihood+
                                    ",subst="+substModel+"," + model +"] = [&U] "+ tree);
                            printedCount++;
                        }
                        count++;

                    }
                    else {
                        line = outFileReader.readLine();
                    }

                }
                else {
                    line = outFileReader.readLine();
                }

            }
            outFileReader.close();


            report(trees, translateHash.keySet().toArray(new String[translateHash.size()]), taxaHash.keySet().toArray(new String[taxaHash.size()]));
        } catch (IOException e) {
            System.err.println("Error Parsing Input log: " + e.getMessage());
        }

    }

    /*
     * Prints the output
     */
    private void report(ArrayList<String> trees, String[] translateList, String[] taxaList) {
        Arrays.sort(taxaList);
        Arrays.sort(translateList);

        System.out.println("#NEXUS\n\nBegin taxa;");
        System.out.println("\tDimensions ntaxa=" + taxaList.length + ";\n\tTaxlabels");
        for(String taxaName : taxaList) {
            System.out.println("\t\t"+taxaName);
        }
        System.out.println("\t\t;\nEnd;\b\bBegin trees;\n\tTranslate");
        for(int i=0; i<translateList.length; i++) {
            System.out.println("\t\t"+(i+1)+" "+translateList[i]);
        }
        for(String tree : trees) {
            for(int i=0; i<translateList.length; i++) {
                int index = -1;
                if((index = tree.indexOf("("+translateList[i]+"[")) > -1 ) {
                    tree = tree.substring(0, index) + "(" + (i+1) + "[" + tree.substring(index+translateList[i].length()+2, tree.length());
                }
                else if((index = tree.indexOf(","+translateList[i]+"[")) > -1 ) {
                    tree = tree.substring(0, index) + "," + (i+1) + "[" + tree.substring(index+translateList[i].length()+2, tree.length());
                }
                else if((index = tree.indexOf(")"+translateList[i]+"[")) > -1 ) {    // This one is for internal nodes
                    tree = tree.substring(0, index) + ")" + (i+1) + "[" + tree.substring(index+translateList[i].length()+2, tree.length());
                }
                //tree.replaceFirst(translateList[i]+"\\[", String.valueOf(i)+"[");
            }
            System.out.println(tree);
        }
        System.out.println("End;");

    }


    public static void printTitle() {
        System.out.println();
        centreLine("MakeAncestralSequenceAnnotatedTree " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("BAli-Phy MCMC Output analysis", 60);
        centreLine("by", 60);
        centreLine("Wai Lok Sibon Li and Marc A. Suchard", 60);
        System.out.println();
        centreLine("David Geffen School of Medicine", 60);
        centreLine("University of California, Los Angeles", 60);
        centreLine("sibonli@ucla.edu", 60);
        centreLine("and", 60);
        centreLine("msuchard@ucla.edu",60);
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

        arguments.printUsage("makeancestralsequenceannotatedtree", "[<*.out file-name> <*.fastas file-name> <*.trees file-name> [<output-file-name>]]");
        System.out.println();
        System.out.println("  Example: ...");
        System.out.println("  Example: ...");
        System.out.println();

    }

    /**
     * Main method
     */
    public static void main(String[] args) throws IOException, TraceException {

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("thin", "thin"),
                        new Arguments.IntegerOption("skip", "skip"),
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

        int skip = -1;
        if (arguments.hasOption("skip")) {
            skip = arguments.getIntegerOption("skip");
        }
        int thin = -1;
        if (arguments.hasOption("thin")) {
            thin = arguments.getIntegerOption("thin");
        }

//        String discreteVariableName = null;
//        if (arguments.hasOption("discreteVariable")) {
//            discreteVariableName = arguments.getStringOption("discreteVariable");
//        }

        String outInputFileName = null;
        String fastasInputFileName = null;
        String treesInputfileName = null;
        String outputFileName = null;

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 4) {
            System.err.println("Unknown option: " + args2[2]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (args2.length > 0) {
            outInputFileName = args2[0];
        }
        if (args2.length > 0) {
            fastasInputFileName = args2[1];
        }
        if (args2.length > 0) {
            treesInputfileName = args2[2];
        }
        if (args2.length > 3) {
            outputFileName = args2[3];
        }

        if (outInputFileName == null) {
            // No input file name was given so throw up a dialog box...
            outInputFileName = Utils.getLoadFileName("MakeAncestralSequenceAnnotatedTree " + version.getVersionString() + " - Select *.out file to analyse");
        }

        if (fastasInputFileName == null) {
            // No input file name was given so throw up a dialog box...
            fastasInputFileName = Utils.getLoadFileName("MakeAncestralSequenceAnnotatedTree " + version.getVersionString() + " - Select *.fastas file to analyse");
        }

        if (treesInputfileName == null) {
            // No input file name was given so throw up a dialog box...
            treesInputfileName = Utils.getLoadFileName("MakeAncestralSequenceAnnotatedTree " + version.getVersionString() + " - Select *.trees file to analyse");
        }

        if(skip==-1) {
            System.out.println("Enter skip: ");
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            skip = Integer.parseInt(br.readLine());
        }
        if(thin==-1) {
            System.out.println("Enter thin: ");
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            thin = Integer.parseInt(br.readLine());
        }
//        if(discreteVariableName==null) {
//            System.out.println("Enter name of discrete variable: ");
//            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
//            discreteVariableName = br.readLine();
//        }

        new MakeAncestralSequenceAnnotatedTree(thin, skip, outInputFileName, fastasInputFileName, treesInputfileName,  outputFileName/*, discreteVariableName, !shortReport, hpds, ess, stdErr, marginalLikelihood*/);

        System.exit(0);
    }
}

