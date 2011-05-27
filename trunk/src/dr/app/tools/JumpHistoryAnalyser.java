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
import java.util.StringTokenizer;
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

    private Variate xData = new Variate.D();
    private Variate yDataMean = new Variate.D();
    private Variate yDataMedian = new Variate.D();
    private Variate yDataUpper = new Variate.D();
    private Variate yDataLower = new Variate.D();

    public JumpHistoryAnalyser(File logFile, int burnin, int binCount, double minTime, double maxTime, double ageOfYoungest)
            throws IOException, ImportException, TraceException {

        this.binCount = binCount;
        this.minTime = minTime;
        this.maxTime = maxTime;
        this.ageOfYoungest = ageOfYoungest;

        double delta = (maxTime - minTime) / (binCount - 1);

        BufferedReader reader = new BufferedReader(new FileReader(logFile));

        String line = reader.readLine();
        while (line != null && line.startsWith("#")) {
            line = reader.readLine();
        }

        // read heading line
        line = reader.readLine();

        while (line != null) {
           String[] columns = line.split("\t");
            if (columns.length == 3) {
                int state = Integer.parseInt(columns[0]);
                int count = (int)Double.parseDouble(columns[1]);

                if (state >= burnin) {
                    Pattern pattern = Pattern.compile("\\[([^:]+):([^:]+):([^\\]]+)\\]");
                    Matcher matcher = pattern.matcher(columns[2]);
                    while (matcher.find()) {
                        String fromState = matcher.group(1);
                        String toState = matcher.group(2);
                        String time = matcher.group(3);
                        System.out.println(fromState + " ->" + toState + ": " + time);
                    }
                }
            }

            line = reader.readLine();
        }

        reader.close();

//        Variate[] bins = new Variate[binCount];
//        double height = 0.0;
//
//        for (int k = 0; k < binCount; k++) {
//            bins[k] = new Variate.D();
//            bins[k].add(popsize);
//        }
//
//        double t = 0.0;
//
//        for (Variate bin : bins) {
//            xData.add(t);
//            if (bin.getCount() > 0) {
//                yDataMean.add(bin.getMean());
//                yDataMedian.add(bin.getQuantile(0.5));
//                yDataLower.add(bin.getQuantile(0.025));
//                yDataUpper.add(bin.getQuantile(0.975));
//            } else {
//                yDataMean.add(Double.NaN);
//                yDataMedian.add(Double.NaN);
//                yDataLower.add(Double.NaN);
//                yDataUpper.add(Double.NaN);
//            }
//            t += delta;
//        }
    }


    public static void main(String[] argv) {
        // command line options to follow shortly...
        try {
            JumpHistoryAnalyser jumpHistory = new JumpHistoryAnalyser(
                    new File("seg1_ha_type2_SRD06_uncorExp_skyride_asym_traits_expPriors_1_mj.jumpHistory.log"),
                    1000000,
                    200,
                    1950,
                    2010,
                    2010
            );
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        } catch (TraceException e) {
            e.printStackTrace();
        }

//        for (int i = 0; i < x.getCount(); i++) {
//            System.out.print(x.get(i));
//
//            for (Variate y : plots) {
//                System.out.print("\t" + y.get(i));
//            }
//            System.out.println();
//        }
    }
}
