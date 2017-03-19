/*
 * SkylineReconstructor.java
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

import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceException;
import dr.inference.trace.TraceList;
import dr.stats.Variate;
import jebl.evolution.coalescent.IntervalList;
import jebl.evolution.coalescent.Intervals;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NewickImporter;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.io.TreeImporter;
import jebl.evolution.trees.RootedTree;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SkylineReconstructor {

    private int binCount;
    private double minTime;
    private double maxTime;
    private double ageOfYoungest;

    private Variate xData = new Variate.D();
    private Variate yDataMean = new Variate.D();
    private Variate yDataMedian = new Variate.D();
    private Variate yDataUpper = new Variate.D();
    private Variate yDataLower = new Variate.D();

    public SkylineReconstructor(File logFile, File treeFile, int burnin,
                                int binCount, double minTime, double maxTime, double ageOfYoungest)
            throws IOException, ImportException, TraceException {

        this.binCount = binCount;
        this.minTime = minTime;
        this.maxTime = maxTime;
        this.ageOfYoungest = ageOfYoungest;

        LogFileTraces traces = new LogFileTraces(logFile.getName(), logFile);
        traces.loadTraces();
        traces.setBurnIn(burnin);

        int stateCount = traces.getStateCount();

        int firstPopSize = findArgument(traces, "popSize");
        int popSizeCount = getTraceRange(traces, firstPopSize);
        int firstGroupSize = findArgument(traces, "groupSize");
        int groupSizeCount = getTraceRange(traces, firstGroupSize);

        boolean isLinear = (groupSizeCount == popSizeCount - 1);
        if (!isLinear && groupSizeCount != popSizeCount) {
            if (isLinear) {
                if (groupSizeCount != popSizeCount - 1) {
                    throw new TraceException("For the stepwise (constant) Bayesian skyline model there should " +
                            "be the same number of group size as population size parameters.");
                }
            }
        }

        ArrayList<ArrayList> popSizes = new ArrayList<ArrayList>();
        ArrayList<ArrayList> groupSizes = new ArrayList<ArrayList>();

        for (int i = 0; i < popSizeCount; i++) {
            popSizes.add(new ArrayList(traces.getValues(firstPopSize + i)));
        }
        for (int i = 0; i < groupSizeCount; i++) {
            groupSizes.add(new ArrayList(traces.getValues(firstGroupSize + i)));
        }

        List heights = traces.getValues(traces.getTraceIndex("treeModel.rootHeight"));
        TraceDistribution distribution = new TraceDistribution(heights,
                traces.getTrace(traces.getTraceIndex("treeModel.rootHeight")).getTraceType());

        double timeMean = distribution.getMean();
        double timeMedian = distribution.getMedian();
        double timeUpper = distribution.getUpperHPD();
        double timeLower = distribution.getLowerHPD();

//        double maxHeight = timeLower;
        double maxHeight = maxTime;
//        switch () {
//            // setting a timeXXXX to -1 means that it won't be displayed...
//            case 0:
//                maxHeight = timeLower;
//                break;
//            case 1:
//                maxHeight = timeMedian;
//                break;
//            case 2:
//                maxHeight = timeMean;
//                break;
//            case 3:
//                maxHeight = timeUpper;
//                break;
//        }

        BufferedReader reader = new BufferedReader(new FileReader(treeFile));

        String line = reader.readLine();

        TreeImporter importer;
        if (line.toUpperCase().startsWith("#NEXUS")) {
            importer = new NexusImporter(reader);
        } else {
            importer = new NewickImporter(reader, false);
        }

        double delta = maxHeight / (binCount - 1);
        int skip = (int) (burnin / traces.getStepSize());
        int state = 0;

        while (importer.hasTree() && state < skip) {
            importer.importNextTree();
            state += 1;
        }

        // the age of the end of this group
        double[][] groupTimes = new double[stateCount][];
        //int tips = 0;
        state = 0;

        while (importer.hasTree()) {
            RootedTree tree = (RootedTree) importer.importNextTree();
            IntervalList intervals = new Intervals(tree);
            int intervalCount = intervals.getIntervalCount();
            //tips = tree.getExternalNodes().size();

            // get the coalescent intervales only
            groupTimes[state] = new double[groupSizeCount];
            double totalTime = 0.0;
            int groupSize = 1;
            int groupIndex = 0;
            int subIndex = 0;
            if (firstGroupSize > 0) {
                double g = (Double) groupSizes.get(groupIndex).get(state);
                if (g != Math.round(g)) {
                    throw new RuntimeException("Group size " + groupIndex + " should be integer but found:" + g);
                } else groupSize = (int) Math.round(g);
            }

            for (int j = 0; j < intervalCount; j++) {

                totalTime += intervals.getInterval(j);

                if (intervals.getIntervalType(j) == IntervalList.IntervalType.COALESCENT) {
                    subIndex += 1;
                    if (subIndex == groupSize) {
                        groupTimes[state][groupIndex] = totalTime;
                        subIndex = 0;
                        groupIndex += 1;
                        if (groupIndex < groupSizeCount) {
                            double g = (Double) groupSizes.get(groupIndex).get(state);
                            if (g != Math.round(g)) {
                                throw new RuntimeException("Group size " + groupIndex + " should be integer but found:" + g);
                            } else groupSize = (int) Math.round(g);
                        }
                    }
                }

                // insert zero-length coalescent intervals
                int diff = intervals.getCoalescentEvents(j) - 1;
                if (diff > 0)
                    throw new RuntimeException("Don't handle multifurcations!");
            }

            state += 1;
        }


        Variate[] bins = new Variate[binCount];
        double height = 0.0;

        for (int k = 0; k < binCount; k++) {
            bins[k] = new Variate.D();

            if (height >= 0.0 && height <= maxHeight) {
                for (state = 0; state < stateCount; state++) {

                    if (isLinear) {
                        double lastGroupTime = 0.0;

                        int index = 0;
                        while (index < groupTimes[state].length && groupTimes[state][index] < height) {
                            lastGroupTime = groupTimes[state][index];
                            index += 1;
                        }

                        if (index < groupTimes[state].length - 1) {
                            double t = (height - lastGroupTime) / (groupTimes[state][index] - lastGroupTime);
                            double p1 = (Double) groupSizes.get(index).get(state);
                            double p2 = (Double) groupSizes.get(index + 1).get(state);
                            double popsize = p1 + ((p2 - p1) * t);
                            bins[k].add(popsize);
                        }
                    } else {
                        int index = 0;
                        while (index < groupTimes[state].length && groupTimes[state][index] < height) {
                            index += 1;
                        }

                        if (index < groupTimes[state].length) {
                            double popSize = (Double) groupSizes.get(index).get(state);
                            if (popSize == 0.0) {
                                throw new RuntimeException("Zero pop size");
                            }

                            bins[k].add(popSize);
                        } else {
                            // Do we really want to do this?
//                                bins[k].add(getPopSize(popSizeCount - 1,state));
                        }
                    }
                }
            }
            height += delta;
        }

        double t = 0.0;

        for (Variate bin : bins) {
            xData.add(t);
            if (bin.getCount() > 0) {
                yDataMean.add(bin.getMean());
                yDataMedian.add(bin.getQuantile(0.5));
                yDataLower.add(bin.getQuantile(0.025));
                yDataUpper.add(bin.getQuantile(0.975));
            } else {
                yDataMean.add(Double.NaN);
                yDataMedian.add(Double.NaN);
                yDataLower.add(Double.NaN);
                yDataUpper.add(Double.NaN);
            }
            t += delta;
        }
    }

    public Variate getXData() {
        return xData;
    }

    public Variate getYDataMean() {
        return yDataMean;
    }

    public Variate getYDataUpper() {
        return yDataUpper;
    }

    public Variate getYDataMedian() {
        return yDataMedian;
    }

    public Variate getYDataLower() {
        return yDataLower;
    }

    private int findArgument(TraceList traceList, String argument) {
        for (int j = 0; j < traceList.getTraceCount(); j++) {
            String statistic = traceList.getTraceName(j);
            String suffix = getNumericalSuffix(statistic);
            if ((suffix.length() == 0 || suffix.equals("1")) && statistic.substring(0, statistic.length() - 1).contains(argument)) {
                return j;
            }
        }
        return -1;
    }

    private String getNumericalSuffix(String argument) {
        int i = argument.length() - 1;

        if (i < 0) return "";

        char ch = argument.charAt(i);

        if (!Character.isDigit(ch)) return "";

        while (i > 0 && Character.isDigit(ch)) {
            i -= 1;
            ch = argument.charAt(i);
        }

        return argument.substring(i + 1, argument.length());
    }

    private int getTraceRange(TraceList traceList, int first) {
        int i = 1;
        int k = first;

        String name = traceList.getTraceName(first);
        String root = name.substring(0, name.length() - 1);
        while (k < traceList.getTraceCount() && traceList.getTraceName(k).equals(root + i)) {
            i++;
            k++;
        }

        return i - 1;
    }

    public static void main(String[] argv) {

        Variate x = null;
        List<Variate> plots = new ArrayList<Variate>();

        for (int i = 1; i <= 200; i++) {
            String stem = "sim" + (i < 10 ? "00" : (i < 100 ? "0" : "")) + i;

            try {
                SkylineReconstructor skyline = new SkylineReconstructor(
                        new File(stem + ".log"),
                        new File(stem + ".trees"),
                        1000000,
                        200,
                        0.0,
                        150000,
                        0.0
                );
                if (x == null) {
                    x = skyline.getXData();
                }
                plots.add(skyline.getYDataMean());

            } catch (IOException e) {
                e.printStackTrace();
            } catch (ImportException e) {
                e.printStackTrace();
            } catch (TraceException e) {
                e.printStackTrace();
            }
            if (i % 10 == 0) {
                System.err.println("Read " + i);
            }
        }

        for (int i = 0; i < x.getCount(); i++) {
            System.out.print(x.get(i));

            for (Variate y : plots) {
                System.out.print("\t" + y.get(i));
            }
            System.out.println();
        }
    }
}
