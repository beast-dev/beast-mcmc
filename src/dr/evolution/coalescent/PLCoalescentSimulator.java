/*
 * PLCoalescentSimulator.java
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

package dr.evolution.coalescent;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class PLCoalescentSimulator {
    private static PrintStream debug = System.err;

    public static void main(String[] arg) throws IOException {

        // READ DEMOGRAPHIC FUNCTION

        String filename = arg[0];
        BufferedReader reader = new BufferedReader(new FileReader(filename));

        double popSizeScale = 1.0;
        double generationTime = 1.0;
        if (arg.length > 2) {
            popSizeScale = Double.parseDouble(arg[2]);
        }
        if (arg.length > 3) {
            generationTime = Double.parseDouble(arg[3]);
        }

        PrintWriter populationFuncLogger = null;
        if (arg.length > 5) {
            String logFileName = arg[5];
            if (logFileName.equals("-")) {
                populationFuncLogger = new PrintWriter(System.out);
            } else {
                populationFuncLogger = new PrintWriter(new FileWriter(logFileName));
            }
        }

        List<Double> times = new ArrayList<Double>();

        String line = reader.readLine();
        String[] tokens = line.trim().split("[\t ]+");
        if (tokens.length < 2) throw new RuntimeException();

        ArrayList<ArrayList> popSizes = new ArrayList<ArrayList>();
        while (line != null) {
            double time = Double.parseDouble(tokens[0]) / generationTime;
            times.add(time);
            
            for (int i = 1; i < tokens.length; i++) {
                popSizes.add(new ArrayList<Double>());
                popSizes.get(i - 1).add(Double.parseDouble(tokens[i]));
            }
            line = reader.readLine();
            if (line != null) {
                tokens = line.trim().split("[\t ]+");
                if (tokens.length != popSizes.size() + 1) throw new RuntimeException();
            }
        }

        reader.close();

        // READ SAMPLE TIMES

        String samplesFilename = arg[1];

        reader = new BufferedReader(new FileReader(samplesFilename));

        line = reader.readLine();
        Taxa taxa = new Taxa();
        int id = 0;
        while (line != null) {

            if (!line.startsWith("#")) {

                tokens = line.split("[\t ]+");

                if (tokens.length == 4) {

                    double t0 = Double.parseDouble(tokens[0]);
                    double t1 = Double.parseDouble(tokens[1]);
                    double dt = Double.parseDouble(tokens[2]);
                    int k = Integer.parseInt(tokens[3]);
                    for (double time = t0; time <= t1; time += dt) {

                        double sampleTime = time / generationTime;
                        for (int i = 0; i < k; i++) {
                            Taxon taxon = new Taxon("t" + id);
                            taxon.setAttribute(dr.evolution.util.Date.DATE, new Date(sampleTime, Units.Type.GENERATIONS, true));
                            taxa.addTaxon(taxon);
                            id += 1;
                        }
                    }

                } else {

                    // sample times are in the same units as simulation
                    double sampleTime = Double.parseDouble(tokens[0]) / generationTime;
                    int count = Integer.parseInt(tokens[1]);

                    for (int i = 0; i < count; i++) {
                        Taxon taxon = new Taxon(id + "");
                        taxon.setAttribute(dr.evolution.util.Date.DATE, new Date(sampleTime, Units.Type.GENERATIONS, true));
                        taxa.addTaxon(taxon);
                        id += 1;
                    }
                }
            }
            line = reader.readLine();
        }

        double minTheta = Double.MAX_VALUE;
        double maxTheta = 0.0;

        PrintWriter out;
        if (arg.length > 4) {
            out = new PrintWriter(new FileWriter(arg[4]));
        } else {
            out = new PrintWriter(System.out);
        }

        int pp = 0;
        for (List<Double> popSize : popSizes) {
            double[] thetas = new double[popSize.size()];
            double[] intervals = new double[times.size() - 1];

            if (populationFuncLogger != null) {
                populationFuncLogger.println("# " + pp);
                ++pp;
            }

            // must reverse the direction of the model
            for (int j = intervals.length; j > 0; j--) {
                intervals[intervals.length - j] = times.get(j) - times.get(j - 1);

                final double theta = popSize.get(j) * popSizeScale;
                thetas[intervals.length - j] = theta;
                if (theta < minTheta) {
                    minTheta = theta;
                }
                if (theta > maxTheta) {
                    maxTheta = theta;
                }

                final double t = times.get(intervals.length) - times.get(j);
                if (populationFuncLogger != null) {
                    populationFuncLogger.println(t + "\t" + theta);
                }
            }

            if (debug != null) {
                debug.println("min theta = " + minTheta);
                debug.println("max theta = " + maxTheta);
            }

            PiecewiseLinearPopulation demo = new PiecewiseLinearPopulation(intervals, thetas, Units.Type.GENERATIONS);

            CoalescentSimulator simulator = new CoalescentSimulator();
            Tree tree = simulator.simulateTree(taxa, demo);

            out.println(TreeUtils.newick(tree));
            if (debug != null) {
                debug.println(TreeUtils.newick(tree));
            }
        }

        if (populationFuncLogger != null) {
            populationFuncLogger.flush();
            populationFuncLogger.close();
        }

        out.flush();
        out.close();
    }
}
