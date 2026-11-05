/*
 * TierFidelityDiagnostics.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package test.dr.app.treelime;

import dr.app.beast.BeastParser;
import dr.app.beast.BeastVersion;
import dr.app.treelime.options.PartitionClockModel;
import dr.app.treelime.options.TopologySearchOptions;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.types.ClockType;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.search.NniNeighbourhood;
import dr.evomodel.tree.search.RootedNniMove;
import dr.evomodel.tree.search.TimeTreeSearch;
import dr.evomodel.tree.search.TrialEvaluator;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.math.MathUtils;
import dr.xml.XMLObject;

import java.io.File;
import java.io.PrintWriter;
import java.io.StringReader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * The tier-fidelity diagnostics of plan 7 §6.2: how well the cheap scores of a candidate NNI (T0: one evaluation after
 * the move; T1: its local block re-fitted) predict the score after all heights and rates are re-fitted (T2). From the
 * fitted starting tree and from trees perturbed by random NNIs and re-fitted, every NNI is scored at all three tiers
 * against its incumbent. Reports the Spearman correlations of dF0 and dF1 with dF2, how far below the incumbent T0 and
 * T1 put the moves that T2 finds to gain (what deltaT1 and deltaKeep must cover), and the recovery dF2 - dF1 of those
 * moves; writes one row per candidate.
 * <p>
 * Usage: TierFidelityDiagnostics data (strict|hmc) states perturbation seed out.tsv [t2Iterations [rateWarmStart]];
 * data is a NEXUS or
 * BEAST XML file; dates come from the file, or from two-digit years ending the names (Dengue4). Run from the
 * repository root with BEAGLE on java.library.path. Not a unit test.
 */
public class TierFidelityDiagnostics {

    private static final double GAIN = 0.01; // moveEpsilon: a move that gains more is accepted

    public static void main(String[] args) throws Exception {
        final String path = args[0];
        final boolean relaxed = args[1].equals("hmc");
        final int states = Integer.parseInt(args[2]);
        final int perturbation = Integer.parseInt(args[3]);
        final long seed = Long.parseLong(args[4]);
        final File out = new File(args[5]);
        final int t2Iterations = args.length > 6 ? Integer.parseInt(args[6]) : 20;
        final String warmStart = args.length > 7 ? args[7] : "preserveSubstitutions";

        TreeLimeOptions options = TreeLimeHarness.importData(path);
        if (!useDatesOfTaxa(options)) {
            TreeLimeHarness.setTipDatesFromNames(options);
        }
        if (relaxed) {
            for (PartitionClockModel clock : options.getPartitionClockModels()) {
                clock.setClockType(ClockType.HMC_CLOCK);
            }
        }
        options.maximizationOptions.printToScreen = false;
        options.maximizationOptions.getTopologySearch().mode = TopologySearchOptions.Mode.NNI;
        File directory = Files.createTempDirectory("tierFidelity").toFile();
        options.fileNameStem = new File(directory, "run").getAbsolutePath();
        String xml = TreeLimeSearchGeneratorTest.generateTo(options, new File(directory, "run.xml"));
        xml = xml.replace("<timeTreeSearch id=\"treeSearch\"",
                "<timeTreeSearch id=\"treeSearch\" searchOnParse=\"false\" t2Iterations=\"" + t2Iterations +
                        "\" rateWarmStart=\"" + warmStart + "\"");
        xml = xml.substring(0, xml.indexOf("<!-- The optimum")) + "</beast>\n";
        java.util.logging.Logger.getLogger("dr").setLevel(java.util.logging.Level.WARNING);
        Map<String, XMLObject> objects = new BeastParser(new String[]{}, null, false, false, false,
                new BeastVersion()).parse(new StringReader(xml), false);
        TreeModel tree = (TreeModel) objects.get("treeModel").getNativeObject();
        TimeTreeSearch search = (TimeTreeSearch) objects.get("treeSearch").getNativeObject();
        MaximizerWrtParameter fin = (MaximizerWrtParameter) objects.get("maximizer").getNativeObject();
        MaximizerWrtParameter round = (MaximizerWrtParameter) objects.get("roundMaximizer").getNativeObject();
        TrialEvaluator evaluator = search.getEvaluator();
        MathUtils.setSeed(seed);

        System.out.println("data " + path + ", " + (relaxed ? "HMC relaxed" : "strict") + " clock, " +
                tree.getExternalNodeCount() + " taxa, " + states + " states, perturbation " + perturbation +
                ", T2 " + t2Iterations + " iterations, rate warm start " + warmStart);
        long started = System.nanoTime();
        fin.maximize(new MaximizerWrtParameter.Budget(2000).setQuiet(true));
        System.out.printf("fitted start: F = %.4f (%.1fs)%n", evaluator.evaluate(), seconds(started));

        List<double[]> rows = new ArrayList<double[]>();
        PrintWriter writer = new PrintWriter(out);
        writer.println("state\tmove\tF_inc\tdF0\tdF1\tdF2\tseconds");
        for (int state = 0; state < states; ++state) {
            if (state > 0) {
                // a nearby worse tree, re-fitted: its incumbent is converged like the search's
                for (int k = 0; k < perturbation; ++k) {
                    List<RootedNniMove> moves = NniNeighbourhood.enumerate(tree, search.getConstraints());
                    moves.get(MathUtils.nextInt(moves.size())).apply(tree, null);
                }
                round.maximize(new MaximizerWrtParameter.Budget(500).setLimits(0.0, 0, 1E-7, 10).setQuiet(true));
            }
            final double fInc = evaluator.evaluate();
            final List<RootedNniMove> moves = NniNeighbourhood.enumerate(tree, search.getConstraints());
            final long stateStarted = System.nanoTime();
            for (RootedNniMove move : moves) {
                long t = System.nanoTime();
                TrialEvaluator.Trial trial = evaluator.trial(move, Double.NEGATIVE_INFINITY, true);
                double[] row = {trial.f0 - fInc, trial.f1 - fInc, trial.f2 - fInc};
                rows.add(row);
                writer.printf("%d\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f%n", state, move.describe(), fInc, row[0], row[1],
                        row[2], seconds(t));
            }
            writer.flush();
            System.out.printf("state %d: F_inc = %.4f, %d moves (%.1fs)%n", state, fInc, moves.size(),
                    seconds(stateStarted));
        }
        writer.close();
        summarize(rows);
        System.out.printf("total %.1fs; rows in %s%n", seconds(started), out);
    }

    static void summarize(List<double[]> rows) {
        final int n = rows.size();
        double[] f0 = new double[n], f1 = new double[n], f2 = new double[n];
        for (int i = 0; i < n; ++i) {
            f0[i] = rows.get(i)[0];
            f1[i] = rows.get(i)[1];
            f2[i] = rows.get(i)[2];
        }
        System.out.printf("candidates: %d; Spearman rho(dF0, dF2) = %.3f, rho(dF1, dF2) = %.3f (target > 0.8)%n", n,
                spearman(f0, f2), spearman(f1, f2));

        List<Double> gain0 = new ArrayList<Double>(), gain1 = new ArrayList<Double>(),
                recovery = new ArrayList<Double>();
        int gainAt0 = 0, gainAt1 = 0;
        for (int i = 0; i < n; ++i) {
            if (f2[i] > GAIN) {
                gain0.add(f0[i]);
                gain1.add(f1[i]);
                recovery.add(f2[i] - f1[i]);
                gainAt0 += f0[i] > GAIN ? 1 : 0;
                gainAt1 += f1[i] > GAIN ? 1 : 0;
            }
        }
        System.out.printf("moves that gain after T2: %d (already at T0: %d, at T1: %d)%n", gain0.size(), gainAt0,
                gainAt1);
        if (!gain0.isEmpty()) {
            System.out.printf("  their dF0: min %.3f, 5%% %.3f, median %.3f  -> deltaT1 must exceed %.3f%n",
                    Collections.min(gain0), quantile(gain0, 0.05), quantile(gain0, 0.5), -Collections.min(gain0));
            System.out.printf("  their dF1: min %.3f, 5%% %.3f, median %.3f  -> deltaKeep must exceed %.3f%n",
                    Collections.min(gain1), quantile(gain1, 0.05), quantile(gain1, 0.5), -Collections.min(gain1));
            double p95 = quantile(recovery, 0.95);
            System.out.printf("  recovery dF2 - dF1: median %.3f, 95%% %.3f -> adaptive deltaKeep %.3f " +
                    "(clamped to [1, 50])%n", quantile(recovery, 0.5), p95, Math.max(1.0, Math.min(50.0, p95)));
        }
        // how far below the incumbent T1 leaves the moves overall, to size the T1 tier
        List<Double> all0 = new ArrayList<Double>();
        for (double x : f0) {
            all0.add(x);
        }
        System.out.printf("all dF0: min %.3f, 5%% %.3f, median %.3f, 95%% %.3f%n", Collections.min(all0),
                quantile(all0, 0.05), quantile(all0, 0.5), quantile(all0, 0.95));
    }

    static double quantile(List<Double> values, double q) {
        List<Double> sorted = new ArrayList<Double>(values);
        Collections.sort(sorted);
        double position = q * (sorted.size() - 1);
        int lower = (int) Math.floor(position);
        int upper = Math.min(sorted.size() - 1, lower + 1);
        return sorted.get(lower) + (position - lower) * (sorted.get(upper) - sorted.get(lower));
    }

    static double spearman(double[] x, double[] y) {
        double[] rx = ranks(x), ry = ranks(y);
        double mx = 0, my = 0;
        for (int i = 0; i < x.length; ++i) {
            mx += rx[i];
            my += ry[i];
        }
        mx /= x.length;
        my /= y.length;
        double sxy = 0, sxx = 0, syy = 0;
        for (int i = 0; i < x.length; ++i) {
            sxy += (rx[i] - mx) * (ry[i] - my);
            sxx += (rx[i] - mx) * (rx[i] - mx);
            syy += (ry[i] - my) * (ry[i] - my);
        }
        return sxy / Math.sqrt(sxx * syy);
    }

    // average ranks for ties
    static double[] ranks(final double[] values) {
        Integer[] order = new Integer[values.length];
        for (int i = 0; i < order.length; ++i) {
            order[i] = i;
        }
        Arrays.sort(order, (a, b) -> Double.compare(values[a], values[b]));
        double[] ranks = new double[values.length];
        int i = 0;
        while (i < order.length) {
            int j = i;
            while (j + 1 < order.length && values[order[j + 1]] == values[order[i]]) {
                j++;
            }
            for (int k = i; k <= j; ++k) {
                ranks[order[k]] = (i + j) / 2.0 + 1;
            }
            i = j + 1;
        }
        return ranks;
    }

    // dates already on the taxa (a BEAST XML file): tip heights from the most recent, as the Tip Dates panel sets them
    static boolean useDatesOfTaxa(TreeLimeOptions options) {
        Date mostRecent = null;
        for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
            Date date = options.taxonList.getTaxon(i).getDate();
            if (date == null) {
                return false;
            }
            if (mostRecent == null || date.after(mostRecent)) {
                mostRecent = date;
            }
        }
        options.useTipDates = true;
        options.maximumTipHeight = 0.0;
        for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
            Taxon taxon = options.taxonList.getTaxon(i);
            double height = mostRecent.getTimeValue() - taxon.getDate().getTimeValue();
            taxon.setAttribute("_height", height);
            options.maximumTipHeight = Math.max(options.maximumTipHeight, height);
        }
        return true;
    }

    private static double seconds(long started) {
        return (System.nanoTime() - started) / 1E9;
    }
}
