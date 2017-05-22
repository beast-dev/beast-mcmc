/*
 * VariableCoalescentSimulator.java
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

package dr.app.vcs;

import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.PiecewiseLinearPopulation;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import ml.options.Options;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * This program simulates a set of coalescent trees from an arbitrary variable population size history expressed as
 * a piecewise-linear function.
 * <p/>
 * Usage:
 * java -jar vcs.jar [options] &lt;infile&gt; &lt;outfile&gt;
 * <p/>
 * Options:
 * <ul>
 * <li> <b> -g &lt;genTime&gt; </b>  sets the generation time to the given value. This is used to scale
 * times in the input file into units of generations.
 * The default value is 1.0
 * <li> <b> -n &lt;numSamples&gt; </b>    sets the number of samples to generate the trees from.
 * The default value is 50
 * <li> <b> -p &lt;popScale&gt; </b>  sets the population size scale factor to the given value. This scale
 * factor transforms the population sizes in the input  file to
 * effective population sizes. This is useful if the population size
 * profiles are expressed, for example, as prevalences, so then scale
 * factor would represent the effective numbers of hosts in the
 * population of interest.
 * The default value is 1.0
 * <li> <b> -se &lt;sampleEnd&gt; </b> the time at which last taxa is sampled, in the time scale provided by the input
 * file. If this differs from the time specified be -ss option then the
 * samples will be evenly spaced between the two times.
 * The default value is the last time if -f is set and the first time
 * otherwise.
 * <li> <b> -ss &lt;sampleStart&gt; </b> the time at which first taxa is sampled, in the time scale provided by the input
 * file. If this differs from the time specified be -se option then the
 * samples will be evenly spaced between the two times.
 * The default value is the last time if -f is set and the first time
 * otherwise.
 * <li> <b> -f </b>          specifies that the population size history proceeds forward in time.
 * Otherwise the population size history is assumed to proceed into the past. <br>
 * <li> <b> -reps &lt;reps&gt; </b> the number of replicate simulations that will be performed.
 * Each replicate will be separated in the output file by a comment line that labels the replicate, e.g. #rep 0.
 * </ul>
 * <p/>
 * <infile>     a whitespace-delimited plain text file. The first column contains the time
 * and should be ascending from zero. Subsequent columns contain the population size
 * histories, for which a tree will be simulated for each.
 * <p/>
 * <outfile>    the file to which the trees will be written in newick format.
 *
 * @author Alexei Drummond
 */
public class VariableCoalescentSimulator {

    public static void main(String[] arg) throws IOException {

        long startTime = System.currentTimeMillis();

        Options options = new Options(arg, 0, 7);

        options.getSet().addOption("g", Options.Separator.EQUALS, Options.Multiplicity.ZERO_OR_ONE);
        options.getSet().addOption("n", Options.Separator.EQUALS, Options.Multiplicity.ZERO_OR_ONE);
        options.getSet().addOption("p", Options.Separator.EQUALS, Options.Multiplicity.ZERO_OR_ONE);
        options.getSet().addOption("se", Options.Separator.EQUALS, Options.Multiplicity.ZERO_OR_ONE);
        options.getSet().addOption("ss", Options.Separator.EQUALS, Options.Multiplicity.ZERO_OR_ONE);
        options.getSet().addOption("reps", Options.Separator.EQUALS, Options.Multiplicity.ZERO_OR_ONE);
        options.getSet().addOption("f", Options.Multiplicity.ZERO_OR_ONE);

        if (!options.check()) {
            System.out.println(options.getCheckErrors());
            System.out.println();
            printUsage();
            System.exit(1);
        }

        double generationTime = 1.0;
        double popSizeScale = 1.0;
        int n = 50;
        double ss = -1;
        double se = -1;
        int reps = 1;

        boolean timeForward = options.getSet().isSet("f");
        if (options.getSet().isSet("g")) {
            String g = options.getSet().getOption("g").getResultValue(0);
            generationTime = Double.parseDouble(g);
            System.out.println("generation time = " + g);
        }
        if (options.getSet().isSet("n")) {
            String sampleSize = options.getSet().getOption("n").getResultValue(0);
            n = Integer.parseInt(sampleSize);
            System.out.println("sample size = " + n);
        }
        if (options.getSet().isSet("p")) {
            String p = options.getSet().getOption("p").getResultValue(0);
            popSizeScale = Double.parseDouble(p);
            System.out.println("population size scale factor = " + p);
        }
        if (options.getSet().isSet("ss")) {
            String sampleStart = options.getSet().getOption("ss").getResultValue(0);
            ss = Double.parseDouble(sampleStart);
            System.out.println("sample start time = " + ss);
        }
        if (options.getSet().isSet("se")) {
            String sampleEnd = options.getSet().getOption("se").getResultValue(0);
            se = Double.parseDouble(sampleEnd);
            System.out.println("sample end time = " + se);
        }
        if (options.getSet().isSet("reps")) {
            String replicates = options.getSet().getOption("reps").getResultValue(0);
            reps = Integer.parseInt(replicates);
            System.out.println("replicates = " + reps);
        }

        String filename = options.getSet().getData().get(0);
        String outfile = options.getSet().getData().get(1);

        // READ DEMOGRAPHIC FUNCTION

        BufferedReader reader = new BufferedReader(new FileReader(filename));

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

        // GENERATE SAMPLE
        double lastTime = times.get(times.size() - 1);

        if (ss == -1) {
            ss = timeForward ? lastTime : times.get(0);
        }
        if (se == -1) {
            se = timeForward ? lastTime : times.get(0);
        }

        double dt = (se - ss) / ((double) n - 1.0);
        double time = ss;

        Taxa taxa = new Taxa();

        for (int i = 0; i < n; i++) {

            double sampleTime;
            if (timeForward) {
                sampleTime = (lastTime - time) / generationTime;
            } else sampleTime = time / generationTime;

            Taxon taxon = new Taxon(i + "");
            taxon.setAttribute(dr.evolution.util.Date.DATE, new Date(sampleTime, Units.Type.GENERATIONS, true));
            taxa.addTaxon(taxon);
            time += dt;
        }

        double minTheta = Double.MAX_VALUE;
        double maxTheta = 0.0;

        PrintWriter out = new PrintWriter(new FileWriter(outfile));

        int popHistory = 0;

        PiecewiseLinearPopulation[] demography = new PiecewiseLinearPopulation[popSizes.size()];

        for (List<Double> popSize : popSizes) {
            double[] thetas = new double[popSize.size()];
            double[] intervals = new double[times.size() - 1];

            for (int i = intervals.length; i >= 0; i--) {

                int j = timeForward ? intervals.length - i : i - 1;
                int k = timeForward ? i : intervals.length - i + 1;

                if (i != 0) intervals[j] = times.get(k) - times.get(k - 1);

                double theta = popSize.get(k) * popSizeScale;
                thetas[j] = theta;
                if (theta < minTheta) {
                    minTheta = theta;
                }
                if (theta > maxTheta) {
                    maxTheta = theta;
                }

                //System.out.println(t + "\t" + theta);
            }

            System.out.println("N" + popHistory + "(t) range = [" + minTheta + ", " + maxTheta + "]");

            demography[popHistory] = new PiecewiseLinearPopulation(intervals, thetas, Units.Type.GENERATIONS);

            popHistory += 1;
        }

        CoalescentSimulator simulator = new CoalescentSimulator();
        for (int i = 0; i < reps; i++) {

            out.println("#rep " + i);
            for (int j = 0; j < demography.length; j++) {
                Tree tree = simulator.simulateTree(taxa, demography[j]);

                out.println(TreeUtils.newick(tree));
                //System.err.println(Tree.Utils.newick(tree));
            }


        }

        out.flush();
        out.close();

        long stopTime = System.currentTimeMillis();
        System.out.println("Took " + (stopTime - startTime) / 1000.0 + " seconds");
    }

    private static Taxa readSampleFile(String fileName, double generationTime) throws IOException {

        BufferedReader reader = new BufferedReader(new FileReader(fileName));

        String line = reader.readLine();
        Taxa taxa = new Taxa();
        int id = 0;
        while (line != null) {

            if (!line.startsWith("#")) {

                String[] tokens = line.split("[\t ]+");

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
            line = reader.readLine();
        }

        return taxa;
    }

    private static void printUsage() {

        System.out.println(
                "Usage: \n" +
                        "  java -jar vcs.jar [options] <infile> <outfile>\n" +
                        " \n" +
                        "Options:\n" +
                        "  -g <value>   sets the generation time to the given value. This is used to scale\n" +
                        "               times in the input file into units of generations. \n" +
                        "               The default value is 1.0\n" +
                        "  -n <int>     sets the number of samples to generate the trees from. \n" +
                        "               The default value is 50\n" +
                        "  -p <value>   sets the population size scale factor to the given value. This scale \n" +
                        "               factor transforms the population sizes in the input  file to\n" +
                        "               effective population sizes. This is useful if the population size\n" +
                        "               profiles are expressed, for example, as prevalences, so then scale\n" +
                        "               factor would represent the effective numbers of hosts in the\n" +
                        "               population of interest. \n" +
                        "               The default value is 1.0\n" +
                        "  -ss <value>  the time at which first taxa is sampled, in the time scale provided by the input\n" +
                        "               file. If this differs from the time specified be -ss option then the \n" +
                        "               samples will be evenly spaced between the two times.\n" +
                        "               The default value is the last time if -f is set and the first time\n" +
                        "               otherwise.\n" +
                        "  -se <value>  the time at which last taxa is sampled, in the time scale provided by the input\n" +
                        "               file. If this differs from the time specified be -se option then the \n" +
                        "               samples will be evenly spaced between the two times.\n" +
                        "               The default value is the last time if -f is set and the first time\n" +
                        "               otherwise.\n" +
                        "  -reps <reps> the number of replicate simulations that will be performed. \n" +
                        "               Each replicate will be separated in the output file by a comment line that\n" +
                        "               labels the replicate, e.g. #rep 0.\n" +
                        "  -f           specifies that the population size history proceeds forward in time.\n" +
                        "               Otherwise the population size history is assumed to proceed into the past.\n" +
                        "\n" +
                        "<infile>  A whitespace-delimited plain text file. The first column contains the time\n" +
                        "          and should be ascending from zero. Subsequent columns contain the population size\n" +
                        "          histories, for which a tree will be simulated for each.\n" +
                        "<outfile> The file to which the trees will be written in newick format.");

    }
}
