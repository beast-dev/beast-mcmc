/*
 * SkylineReconstructor.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.tools;

import dr.app.util.Arguments;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceException;
import dr.inference.trace.TraceList;
import dr.stats.Variate;
import dr.util.DataTable;
import jam.console.ConsoleApplication;
import jebl.evolution.coalescent.IntervalList;
import jebl.evolution.coalescent.Intervals;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NewickImporter;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.io.TreeImporter;
import jebl.evolution.trees.RootedTree;

import javax.swing.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class JumpHistoryAnalyser {

    private int binCount;
    private double minTime;
    private double maxTime;
    private double ageOfYoungest;

    public JumpHistoryAnalyser(String inputFileName, Set<String> fromStates, Set<String> toStates, boolean iterateFrom, boolean iterateTo, int burnin, int binCount, double minTime, double maxTime, double ageOfYoungest)
            throws IOException, ImportException, TraceException {

        this.binCount = binCount;
        this.minTime = minTime;
        this.maxTime = maxTime;
        this.ageOfYoungest = ageOfYoungest;

        double delta = (maxTime - minTime) / (binCount - 1);

        Set<Set<String>> fromSets = new HashSet<Set<String>>();

        // if a particular set of from states is specified then use them
        if (fromStates.size() > 0) {
            if (iterateFrom) {
                for (String state : fromStates) {
                    Set<String> set = new HashSet<String>();
                    set.add(state);
                    fromSets.add(set);
                }
            } else {
                fromSets.add(fromStates);
            }
        }

        // if a particular set of to states is specified then use them
        Set<Set<String>> toSets = new HashSet<Set<String>>();
        if (toStates.size() > 0) {
            if (iterateTo) {
                for (String state : toStates) {
                    Set<String> set = new HashSet<String>();
                    set.add(state);
                    toSets.add(set);
                }
            } else {
                toSets.add(toStates);
            }
        }

        Pattern pattern = Pattern.compile("\\[([^:]+):([^:]+):([^\\]]+)\\]");
        Set<String> allStates = new HashSet<String>();

        if (fromSets.size() == 0 || toSets.size() == 0) {
            // from or to states have not been specified so parse file to get full set.
            BufferedReader reader = new BufferedReader(new FileReader(inputFileName));

            String line = reader.readLine();
            while (line != null && line.startsWith("#")) {
                line = reader.readLine();
            }

            // read heading line
            line = reader.readLine();

            // read next line
            line = reader.readLine();

            while (line != null) {
                String[] columns = line.split("\t");
                if (columns.length == 3) {
                    Matcher matcher = pattern.matcher(columns[2]);
                    while (matcher.find()) {
                        String fromState = matcher.group(1);
                        allStates.add(fromState);
                        String toState = matcher.group(2);
                        allStates.add(toState);
                    }
                }

                line = reader.readLine();
            }

            reader.close();
        }

        if (fromStates.size() == 0) {
            // if particular from states have not been specified, use them all
            if (iterateFrom) {
                for (String state : allStates) {
                    Set<String> set = new HashSet<String>();
                    set.add(state);
                    fromSets.add(set);
                }
            } else {
                fromSets.add(allStates);
            }
        }

        if (toStates.size() == 0) {
            // if particular to states have not been specified, use them all
            if (iterateTo) {
                for (String state : allStates) {
                    Set<String> set = new HashSet<String>();
                    set.add(state);
                    toSets.add(set);
                }
            } else {
                toSets.add(allStates);
            }
        }

        int fromSetCount = fromSets.size();
        int toSetCount = fromSets.size();

        int[][][] bins = new int[fromSetCount][toSetCount][binCount];

        int index = 0;
        Map<Set<String>, Integer> fromSetIndices = new HashMap<Set<String>, Integer>();
        for (Set<String> fromSet : fromSets) {
            fromSetIndices.put(fromSet, index);
            index++;
        }

        index = 0;
        Map<Set<String>, Integer> toSetIndices = new HashMap<Set<String>, Integer>();
        for (Set<String> toSet : toSets) {
            toSetIndices.put(toSet, index);
            index++;
        }

        BufferedReader reader = new BufferedReader(new FileReader(inputFileName));

        String line = reader.readLine();
        while (line != null && line.startsWith("#")) {
            line = reader.readLine();
        }

        // read heading line
        line = reader.readLine();

        // read next line
        line = reader.readLine();

        while (line != null) {
            String[] columns = line.split("\t");
            if (columns.length == 3) {
                int state = Integer.parseInt(columns[0]);
                int count = (int)Double.parseDouble(columns[1]);

                if (state >= burnin) {
                    Matcher matcher = pattern.matcher(columns[2]);
                    while (matcher.find()) {
                        String fromState = matcher.group(1);
                        String toState = matcher.group(2);
                        String timeString = matcher.group(3);

                        Set<String> fromSet = null;
                        for (Set<String> set : fromSets) {
                            if (set.contains(fromState)) {
                                fromSet = set;
                                break;
                            }
                        }
                        if (fromSet != null) {
                            int fromIndex = fromSetIndices.get(fromSet);

                            Set<String> toSet = null;
                            for (Set<String> set : toSets) {
                                if (set.contains(toState)) {
                                    toSet = set;
                                    break;
                                }
                            }

                            if (toSet != null) {

                                int toIndex = toSetIndices.get(toSet);

                                double time = Double.parseDouble(timeString);
                                time = ageOfYoungest - time;
                                double binTime = minTime;
                                int bin = 0;
                                while (time > binTime) {
                                    binTime += delta;
                                    bin ++;
                                }
                                bins[fromIndex][toIndex][bin] ++;
                            }
                        }
//                        System.out.println(fromState + " ->" + toState + ": " + timeString);
                    }
                }
            }

            line = reader.readLine();
        }

        reader.close();

        System.out.print("time");
        for (Set<String> fromSet : fromSets) {
            String fromLabel = "all";
            if (fromSet.size() == 1) {
                fromLabel = fromSet.toArray(new String[1])[0];
            } else if (!fromSet.equals(allStates)) {
                StringBuilder sb = new StringBuilder("{");
                boolean first = true;
                for (String state : fromSet) {
                    if (!first) {
                        sb.append(",");
                    } else {
                        first = false;
                    }
                    sb.append(state);
                }
                sb.append("}");
                fromLabel = sb.toString();
            }

            for (Set<String> toSet : toSets) {
                String toLabel = "all";
                if (toSet.size() == 1) {
                    toLabel = toSet.toArray(new String[1])[0];
                } else if (!toSet.equals(allStates)) {
                    StringBuilder sb = new StringBuilder("{");
                    boolean first = true;
                    for (String state : toSet) {
                        if (!first) {
                            sb.append(",");
                        } else {
                            first = false;
                        }
                        sb.append(state);
                    }
                    sb.append("}");
                    toLabel = sb.toString();
                }

                if (!fromLabel.equals(toLabel)) {
                    System.out.print("\t" + fromLabel + "->" + toLabel);
                }
            }
        }
        System.out.println();

        double time = minTime;
        for (int bin = 0; bin < binCount; bin++) {
            System.out.print(time);

            for (Set<String> fromSet : fromSets) {
                int fromIndex = fromSetIndices.get(fromSet);

                for (Set<String> toSet : toSets) {
                    int toIndex = toSetIndices.get(toSet);
                    if (fromIndex != toIndex) {
                        System.out.print("\t" + bins[fromIndex][toIndex][bin]);
                    }
                }
            }
            System.out.println();

            time += delta;
        }
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("jumphistoryanalyser", "<input-file-name> [<output-file-name>]");
        System.out.println();
        System.out.println("  Example: jumphistoryanalyser -burnin 100000 -min 1950 -max 2010 -mrsd 2009.9 jumps.txt");
        System.out.println("  Example: jumphistoryanalyser -burnin 100000 -from \"usa\" -to \"uk,fr\" -min 1950 -max 2010 -mrsd 2009.9 jumps.txt");
        System.out.println();
    }


    public static void main(String[] args) {
        String inputFileName = null;


//        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in'"),
                        new Arguments.StringOption("from", "from_states", "set of 'from' states to limit the history [default all states]"),
                        new Arguments.StringOption("to", "to_states", "set of 'to' states to limit the history [default all states]"),
                        new Arguments.Option("iterateFrom", "iterate over 'from' states [default combine states]"),
                        new Arguments.Option("iterateTo", "iterate over 'to' states [default combine states]"),
                        new Arguments.IntegerOption("bins", "the number of discrete bins [default 100]"),
                        new Arguments.RealOption("min", "the minimum bound of the time range"),
                        new Arguments.RealOption("max", "the maximum bound of the time range"),
                        new Arguments.RealOption("mrsd", "the date of the most recently sampled tip"),
                        new Arguments.StringOption("to", "to_states", "set of 'to' states to limit the history [default all states]"),
                        new Arguments.Option("help", "option to print this message"),
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.err.println(ae);
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

        int binCount = 100;
        if (arguments.hasOption("bins")) {
            burnin = arguments.getIntegerOption("bins");
        }

        double minTime = arguments.getRealOption("min");
        double maxTime = arguments.getRealOption("max");

        if (minTime >= maxTime) {
            System.err.println("The minimum time must be less than the maximum time");
            printUsage(arguments);
            System.exit(1);
        }
        double mrsd = arguments.getRealOption("mrsd");

        Set<String> fromStates = new HashSet<String>();
        Set<String> toStates = new HashSet<String>();

        if (arguments.hasOption("from")) {
            String stateString = arguments.getStringOption("from");
            String[] states = stateString.split("[\"\\s,]");
            for (String state : states) {
                if (state.length() > 0) {
                    fromStates.add(state);
                }
            }
        }

        if (arguments.hasOption("to")) {
            String stateString = arguments.getStringOption("to");
            String[] states = stateString.split("[\"\\s,]");
            for (String state : states) {
                if (state.length() > 0) {
                    toStates.add(state);
                }
            }
        }

        boolean iterateFrom = arguments.hasOption("iterateFrom");
        boolean iterateTo = arguments.hasOption("iterateTo");


        final String[] args2 = arguments.getLeftoverArguments();

        switch (args2.length) {
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

        // command line options to follow shortly...
        try {
            JumpHistoryAnalyser jumpHistory = new JumpHistoryAnalyser(
                    inputFileName,
                    fromStates,
                    toStates,
                    iterateFrom,
                    iterateTo,
                    burnin,
                    binCount,
                    minTime,
                    maxTime,
                    mrsd
            );
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        } catch (TraceException e) {
            e.printStackTrace();
        }

        System.exit(0);
    }
}

