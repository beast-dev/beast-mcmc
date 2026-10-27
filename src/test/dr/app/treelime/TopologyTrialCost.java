/*
 * TopologyTrialCost.java
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
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.types.ClockType;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.xml.XMLObject;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Measures the cost of a T0 trial of the topology search (plan 7 M0): store the model state, apply a rooted NNI at
 * fixed heights, evaluate the joint density, restore; split into store, edit, tree likelihood, prior and restore, and
 * compared with a full evaluation. Every feasible NNI of the starting tree is timed, after a warm-up.
 * <p>
 * Usage: TopologyTrialCost data.nex (dates|none) (strict|relaxed) repetitions; run from the repository root with
 * BEAGLE on java.library.path. Not a unit test.
 */
public class TopologyTrialCost {

    public static void main(String[] args) throws Exception {
        String path = args[0];
        boolean dates = args[1].equals("dates");
        boolean relaxed = args[2].equals("relaxed");
        TreeLimeOptions options = TreeLimeHarness.importData(path);
        if (dates) {
            TreeLimeHarness.setTipDatesFromNames(options);
        }
        if (relaxed) {
            for (PartitionClockModel clock : options.getPartitionClockModels()) {
                clock.setClockType(ClockType.HMC_CLOCK);
            }
        }
        options.maximizationOptions.printToScreen = false;
        String xml = TreeLimeHarness.generate(options);
        xml = xml.replace("<maximizeWrtParameter id=\"maximizer\"",
                "<maximizeWrtParameter id=\"maximizer\" maximizeOnParse=\"false\"");
        xml = xml.substring(0, xml.indexOf("<!-- The optimum")) + "</beast>\n";
        java.util.logging.Logger.getLogger("dr").setLevel(java.util.logging.Level.WARNING);
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        Map<String, XMLObject> objects = parser.parse(new StringReader(xml), false);
        TreeModel tree = (TreeModel) objects.get("treeModel").getNativeObject();
        Likelihood joint = (Likelihood) objects.get("joint").getNativeObject();
        Likelihood likelihood = (Likelihood) objects.get("likelihood").getNativeObject();
        Likelihood prior = (Likelihood) objects.get("prior").getNativeObject();
        Model model = ((CompoundLikelihood) joint).getModel();
        int n = tree.getExternalNodeCount();

        // feasible NNIs: {P, G, a, U}
        List<int[]> moves = new ArrayList<int[]>();
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef p = tree.getInternalNode(i);
            if (tree.isRoot(p)) continue;
            NodeRef g = tree.getParent(p);
            NodeRef u = tree.getChild(g, 0) == p ? tree.getChild(g, 1) : tree.getChild(g, 0);
            if (tree.getNodeHeight(u) < tree.getNodeHeight(p)) {
                for (int c = 0; c < 2; c++) {
                    moves.add(new int[]{p.getNumber(), g.getNumber(), tree.getChild(p, c).getNumber(), u.getNumber()});
                }
            }
        }
        joint.getLogLikelihood();

        int reps = Integer.parseInt(args[3]);
        // warm-up
        for (int r = 0; r < 2; r++) for (int[] m : moves) trial(tree, model, joint, likelihood, prior, m, null);
        for (int r = 0; r < 50; r++) { joint.makeDirty(); joint.getLogLikelihood(); }

        double[] full = new double[200];
        double[] fullLik = new double[200];
        double[] fullPrior = new double[200];
        for (int r = 0; r < full.length; r++) {
            joint.makeDirty();
            long t0 = System.nanoTime();
            likelihood.getLogLikelihood();
            long t1 = System.nanoTime();
            prior.getLogLikelihood();
            long t2 = System.nanoTime();
            full[r] = (t2 - t0) / 1E3; fullLik[r] = (t1 - t0) / 1E3; fullPrior[r] = (t2 - t1) / 1E3;
        }

        List<double[]> phases = new ArrayList<double[]>();
        double depthSum = 0;
        for (int r = 0; r < reps; r++) {
            for (int[] m : moves) {
                double[] t = new double[6];
                trial(tree, model, joint, likelihood, prior, m, t);
                phases.add(t);
                if (r == 0) depthSum += depth(tree, tree.getNode(m[1])) + 1;
            }
        }
        double meanDepth = depthSum / moves.size();
        String[] names = {"store", "edit", "likelihood", "prior", "restore", "total"};
        System.out.printf("%s n=%d %s %s: %d NNIs, mean nodes recomputed d=%.1f (d/(n-1)=%.2f)%n", path, n,
                dates ? "dated" : "contemporaneous", relaxed ? "HMC relaxed" : "strict", moves.size(), meanDepth,
                meanDepth / (n - 1));
        System.out.printf("  full evaluation (makeDirty): median %.1f us (likelihood %.1f, prior %.1f)%n",
                median(full), median(fullLik), median(fullPrior));
        StringBuilder sb = new StringBuilder("  T0 trial medians (us):");
        for (int k = 0; k < names.length; k++) {
            double[] v = new double[phases.size()];
            for (int i = 0; i < v.length; i++) v[i] = phases.get(i)[k];
            sb.append(String.format(" %s %.1f", names[k], median(v)));
        }
        System.out.println(sb);
        double[] tot = new double[phases.size()];
        for (int i = 0; i < tot.length; i++) tot[i] = phases.get(i)[5];
        System.out.printf("  T0 / full = %.2f; T0 trials per second = %.0f%n", median(tot) / median(full), 1E6 / median(tot));
    }

    static int depth(TreeModel tree, NodeRef node) {
        int d = 0;
        for (NodeRef x = node; !tree.isRoot(x); x = tree.getParent(x)) d++;
        return d + 1;
    }

    static void trial(TreeModel tree, Model model, Likelihood joint, Likelihood likelihood, Likelihood prior, int[] m,
                      double[] t) {
        long t0 = System.nanoTime();
        model.storeModelState();
        long t1 = System.nanoTime();
        NodeRef p = tree.getNode(m[0]), g = tree.getNode(m[1]), a = tree.getNode(m[2]), u = tree.getNode(m[3]);
        tree.beginTreeEdit();
        tree.removeChild(p, a);
        tree.removeChild(g, u);
        tree.addChild(p, u);
        tree.addChild(g, a);
        tree.endTreeEdit();
        long t2 = System.nanoTime();
        likelihood.getLogLikelihood();
        long t3 = System.nanoTime();
        prior.getLogLikelihood();
        joint.getLogLikelihood();
        long t4 = System.nanoTime();
        model.restoreModelState();
        long t5 = System.nanoTime();
        if (t != null) {
            t[0] = (t1 - t0) / 1E3; t[1] = (t2 - t1) / 1E3; t[2] = (t3 - t2) / 1E3; t[3] = (t4 - t3) / 1E3;
            t[4] = (t5 - t4) / 1E3; t[5] = (t5 - t0) / 1E3;
        }
        joint.getLogLikelihood();
    }

    static double median(double[] v) {
        double[] s = v.clone();
        Arrays.sort(s);
        return s[s.length / 2];
    }
}
