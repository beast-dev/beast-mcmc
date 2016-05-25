/*
 * GetAncestralSequenceFromSplitTrait.java
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
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeExporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.FlexibleTree;
import dr.inference.trace.TraceException;
import dr.util.Version;
import java.util.*;
import java.io.*;

public class GetAncestralSequenceFromSplitTrait {

    private final static Version version = new BeastVersion();


    public GetAncestralSequenceFromSplitTrait(String treeAnnotatorFileName, String outputFileName) throws IOException, TraceException {

        File outInputFile = new File(treeAnnotatorFileName);

        if (outInputFile.isFile()) {
            System.out.println("Analysing *.tree file: " + treeAnnotatorFileName);
        } else {
            System.err.println(treeAnnotatorFileName + " does not exist!");
            System.exit(0);
        }

        if (outputFileName != null) {
            FileOutputStream outputStream = new FileOutputStream(outputFileName);
            System.setOut(new PrintStream(outputStream));
        }

        analyze(outInputFile);
    }

    /**
     * Recursively analyzes log files.
     *
     * @param treeAnnotatorInputFile       the file to analyze (if this is a directory then the files within it are analyzed)
     * @throws dr.inference.trace.TraceException
     *          if the trace file is in the wrong format or corrupted
     */
    private void analyze(File treeAnnotatorInputFile) throws TraceException {

        try {

            FileReader fileReader = new FileReader(treeAnnotatorInputFile);
            TreeImporter importer = new NexusImporter(fileReader);
            FlexibleTree tree = (FlexibleTree) importer.importNextTree();
            for(int i = 0; i<tree.getNodeCount(); i++) {
                Hashtable<Integer, State> states = new Hashtable<Integer, State>();
                for (Iterator<String> j = tree.getNodeAttributeNames(tree.getNode(i)); j.hasNext();) {
                    String name = j.next();
                    if(name.indexOf("states_")>=0) {
                        Integer d = Integer.parseInt(name.replaceFirst("states_", "").replaceFirst("\\..+", ""));
                        State s; //= new State(name.);
                        if(states.containsKey(d)) {
                            s = states.get(d);
                        }
                        else {
                            s = new State(d);
                        }

//                        if (tree.getNodeAttribute(tree.getNode(0), name) instanceof Object[]) {
//                            Object[] o = (Object[]) tree.getNodeAttribute(tree.getNode(0), name);
//
//                            if(name.matches("states_"+d+".prob")) {
//                                double[] probabilities = new double[o.length];
//                                for(int k=0; k<o.length; k++) {
//                                    probabilities[k] = (Double) o[k];
//                                }
//                                s.setProbabilities(probabilities);
//                            }
//                            else if(name.matches("states_"+d)) {
//                                String[] set = new String[o.length];
//                                for(int k=0; k<o.length; k++) {
//                                    set[k] = (String) o[k];
//                                }
//                                s.setSet(set);
//                            }
//                        }
//                        else {

                        if(name.matches("states_"+d+".prob")) {
                            Object o = tree.getNodeAttribute(tree.getNode(i), name);
                            double probability = (Double) o;
                            s.setProbability(probability);
                        }
                        else if(name.matches("states_"+d)) {
                            Object o = tree.getNodeAttribute(tree.getNode(i), name);
                            String value = (String) o;
                            s.setState(value.replaceAll("\"", ""));
                        }
                        else if(name.matches("states_"+d+".set.prob")) { /* Not necessary but lets parse it anyways */
                            Object[] o = (Object[]) tree.getNodeAttribute(tree.getNode(i), name);
                            double[] probabilities = new double[o.length];
                            for(int k=0; k<o.length; k++) {
                                probabilities[k] = (Double) o[k];
                            }
                            s.setProbabilities(probabilities);
                        }
                        else if(name.matches("states_"+d+".set")) { /* Not necessary but lets parse it anyways */
                            Object[] o = (Object[]) tree.getNodeAttribute(tree.getNode(i), name);
                            String[] set = new String[o.length];
                            for(int k=0; k<o.length; k++) {
                                set[k] = ((String) o[k]).replaceAll("\"", "");
                            }
                            s.setSet(set);
                        }
//                        }
                        states.put(d, s);
                    }
                }

                State[] statesArray = states.values().toArray(new State[states.size()]);
                Arrays.sort(statesArray);
                /* Set the default length to the number of characters that it would need */
                StringBuffer sb = new StringBuffer(statesArray.length * statesArray[0].getState().length());
                for(State s : statesArray) {
                    sb.append(s.getState());
                }

                tree.setNodeAttribute(tree.getNode(i), "seq", sb.toString());
            }
            /* Export the new tree with the new sequences */
            TreeExporter exporter = new NexusExporter(System.out);
            exporter.exportTree(tree);


            System.out.println("Begin trees;");
            System.out.println("\ttree max_tree = " + tree.toString());
            System.out.println("End;");

        } catch (IOException e) {
            System.err.println("Error Parsing Input log: " + e.getMessage());
        }
        catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
        }

    }



    public static void printTitle() {
        System.out.println();
        centreLine("GetAncestralSequenceFromSplitTrait " + version.getVersionString() + ", " + version.getDateString(), 60);
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

        arguments.printUsage("GetAncestralSequenceFromSplitTrait", "[<*.tree file-name> [<output-file-name>]]");
        System.out.println();
        System.out.println("  Example: ...");
        System.out.println("  Example: ...");
        System.out.println();

    }

    //Main method
    public static void main(String[] args) throws IOException, TraceException {

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
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

        String outInputFileName = null;
        String outputFileName = null;
        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 2) {
            System.err.println("Unknown option: " + args2[2]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (args2.length > 0) {
            outInputFileName = args2[0];
        }
        if (args2.length > 1) {
            outputFileName = args2[1];
        }

        if (outInputFileName == null) {
            // No input file name was given so throw up a dialog box...
            outInputFileName = Utils.getLoadFileName("GetAncestralSequenceFromSplitTrait " + version.getVersionString() + " - Select *.tree file to analyse");
        }
        new GetAncestralSequenceFromSplitTrait(outInputFileName, outputFileName);

        System.exit(0);
    }


    private class State implements Comparable<State> {
        private State(int stateIndex) {
            this.stateIndex = stateIndex;
        }

        private int stateIndex;
        private String[] set;
        private double[] probabilities;
        private String state;
        private double probability;

        public String getState() {
            return state;
        }

        public void setState(String state) {
            this.state = state;
        }

        public double getProbability() {
            return probability;
        }

        public void setProbability(double probability) {
            this.probability = probability;
        }

        public int getStateIndex() {
            return stateIndex;
        }

        public void setStateIndex(int stateIndex) {
            this.stateIndex = stateIndex;
        }

        public String[] getSet() {
            return set;
        }

        public void setSet(String[] set) {
            this.set = set;
        }

        public double[] getProbabilities() {
            return probabilities;
        }

        public void setProbabilities(double[] probabilities) {
            this.probabilities = probabilities;
        }


        public int compareTo(State o) {
            return this.getStateIndex()-o.getStateIndex();
        }

//        public String getMaxProbState() {
//            if(set.length != probabilities.length) {
//                throw new RuntimeException("Set and probabilities arrays are not the same length: " + set.length +
//                        ", " + probabilities.length + ". State: " + getStateIndex());
//            }
//            if(set.length ==0 || probabilities.length ==0) {
//                throw new RuntimeException("Array lengths are zero");
//            }
//            double max = 0;
//            int maxIndex =0;
//            for(int i=0; i<probabilities.length; i++) {
//                if(probabilities[i]>max) {
//                    maxIndex=i;
//                }
//            }
//            return set[maxIndex];
//        }
    }
}

