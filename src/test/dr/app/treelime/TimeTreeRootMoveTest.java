/*
 * TimeTreeRootMoveTest.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.search.BranchRateWarmStart;
import dr.evomodel.tree.search.LocalOptimizer;
import dr.evomodel.tree.search.RootChangingSprMove;
import dr.evomodel.tree.search.TimeTreeSearch;
import dr.evomodel.tree.search.TopologyConstraints;
import dr.evomodel.tree.search.TopologyHash;
import dr.evomodel.tree.search.TreeStateSnapshot;
import dr.evomodel.tree.search.TrialEvaluator;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static test.dr.app.treelime.TimeTreeSearchTest.Setup;
import static test.dr.app.treelime.TimeTreeSearchTest.joint;
import static test.dr.app.treelime.TimeTreeSearchTest.parse;
import static test.dr.app.treelime.TimeTreeSearchTest.rates;
import static test.dr.app.treelime.TimeTreeSearchTest.roundMaximizer;
import static test.dr.app.treelime.TimeTreeSearchTest.tree;

/**
 * Root-changing SPR (plan 7 M6, §5.4) on the Dengue XML of TimeTreeSearchTest: V6b (every node that keeps a branch
 * keeps its rate through the root change, whose index the rate model shifts; trials leave no trace, gradients
 * included, and replays give the trial's F), the constraint rule against monophyly over every root move, and the
 * search with root moves.
 */
public class TimeTreeRootMoveTest extends TestCase {

    public TimeTreeRootMoveTest(String name) {
        super(name);
    }

    // V6b: each move gives the root the expected node, keeps a valid time tree and keeps every node's rate (the old
    // root's new branch, or P's, takes the new root's old rate); TreeModel's restore undoes it all
    public void testApplyKeepsRatesByNode() throws Exception {
        for (boolean dated : new boolean[]{true, false}) {
            Map<String, XMLObject> objects = parse(new Setup().relaxed(true).dated(dated).xml());
            if (objects == null) {
                return;
            }
            TreeModel tree = tree(objects);
            Likelihood joint = joint(objects);
            Model model = joint.getModel();
            ArbitraryBranchRates branchRates = rates(objects);
            Parameter raw = branchRates.getRateParameter();
            for (int i = 0; i < raw.getDimension(); i++) {
                raw.setParameterValue(i, 0.5 + 0.05 * i);
            }
            final String newick = TreeUtils.newick(tree);
            final double f = joint.getLogLikelihood();
            final double[] raw0 = raw.getParameterValues();
            int above = 0, atRoot = 0;
            for (RootChangingSprMove move : RootChangingSprMove.enumerate(tree, 3, null)) {
                assertTrue(move.describe(), move.isValid(tree));
                final int oldRoot = tree.getRoot().getNumber();
                final double[] before = ratesByNode(tree, branchRates);
                joint.getLogLikelihood();
                model.storeModelState();
                move.apply(tree, null);
                final int newRoot = move.getKind() == RootChangingSprMove.Kind.ABOVE_ROOT ? move.getP() : move.getC();
                assertEquals(move.describe(), newRoot, tree.getRoot().getNumber());
                assertValidTimeTree(tree);
                final double[] after = ratesByNode(tree, branchRates);
                final int gained = move.getKind() == RootChangingSprMove.Kind.ABOVE_ROOT ? oldRoot : move.getP();
                for (int n = 0; n < tree.getNodeCount(); n++) {
                    if (n == newRoot) {
                        assertTrue(Double.isNaN(after[n]));
                    } else if (n == gained) {
                        assertEquals(move.describe(), before[newRoot], after[n], 0.0);
                    } else {
                        assertEquals(move.describe() + " node " + n, before[n], after[n], 0.0);
                    }
                }
                assertTrue(joint.getLogLikelihood() > Double.NEGATIVE_INFINITY);
                model.restoreModelState();
                assertEquals(oldRoot, tree.getRoot().getNumber());
                assertEquals(newick, TreeUtils.newick(tree));
                assertArrayEquals(raw0, raw.getParameterValues(), 0.0);
                assertEquals(f, joint.getLogLikelihood(), 1E-10 * Math.abs(f));
                if (move.getKind() == RootChangingSprMove.Kind.ABOVE_ROOT) {
                    above++;
                } else {
                    atRoot++;
                }
            }
            assertEquals(tree.getNodeCount() - 3, above);
            assertTrue(atRoot > 0);
        }
    }

    // the rule refuses exactly the root moves after which a constrained set is no longer a clade, for every clade of
    // the starting tree as the constraint
    public void testConstraintRule() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().xml());
        if (objects == null) {
            return;
        }
        TreeModel tree = tree(objects);
        Likelihood joint = joint(objects);
        Model model = joint.getModel();
        BitSet[] clades = TopologyConstraints.clades(tree);
        List<RootChangingSprMove> all = RootChangingSprMove.enumerate(tree, 3, null);
        int refused = 0, tested = 0;
        for (int node = 0; node < tree.getNodeCount(); node++) {
            final int size = clades[node].cardinality();
            if (size < 2 || size >= tree.getExternalNodeCount()) {
                continue;
            }
            Set<String> taxa = new HashSet<String>();
            for (int tip = clades[node].nextSetBit(0); tip >= 0; tip = clades[node].nextSetBit(tip + 1)) {
                taxa.add(tree.getNodeTaxon(tree.getNode(tip)).getId());
            }
            TopologyConstraints constraints = new TopologyConstraints(tree, null,
                    Collections.<Set<String>>singletonList(taxa));
            List<RootChangingSprMove> allowed = RootChangingSprMove.enumerate(tree, 3, constraints);
            for (RootChangingSprMove move : all) {
                final boolean violates = constraints.violates(tree, move, clades);
                joint.getLogLikelihood();
                model.storeModelState();
                move.apply(tree, null);
                assertEquals(move.describe() + " " + taxa, violates, !TreeUtils.isMonophyletic(tree, taxa));
                model.restoreModelState();
                assertEquals(!violates, allowed.contains(move));
                refused += violates ? 1 : 0;
                tested++;
            }
        }
        assertTrue(refused > 0 && refused < tested);
    }

    // V6b: trials of root moves (T0, T1, T2) leave the tree, the parameters, the density and a direct gradient call as
    // they were; accepted ones replay to the trial's F, and a snapshot from before restores by clade
    public void testTrialsAndReplays() throws Exception {
        for (boolean relaxed : new boolean[]{false, true}) {
            Map<String, XMLObject> objects = parse(new Setup().relaxed(relaxed).xml());
            if (objects == null) {
                return;
            }
            TreeModel tree = tree(objects);
            Likelihood joint = joint(objects);
            TimeTreeSearchTest.perturb(tree, 2);
            TrialEvaluator evaluator = new TrialEvaluator(tree, joint, roundMaximizer(objects),
                    new BranchRateWarmStart(tree, relaxed ? rates(objects) : null,
                            BranchRateWarmStart.Mode.PRESERVE_SUBSTITUTIONS),
                    new LocalOptimizer(tree, relaxed ? rates(objects) : null, 2), 5, 1);
            GradientWrtParameterProvider gradient =
                    (GradientWrtParameterProvider) objects.get("gradient.nodeHeights").getNativeObject();
            final String newick = TreeUtils.newick(tree);
            final double f = joint.getLogLikelihood();
            final double[] parameters = values(evaluator.getRoundParameters());
            final double[] heightGradient = gradient.getGradientLogDensity().clone();
            final int root = tree.getRoot().getNumber();

            int k = 0;
            for (RootChangingSprMove move : RootChangingSprMove.enumerate(tree, 3, null)) {
                TrialEvaluator.Trial trial = evaluator.trial(move, Double.NEGATIVE_INFINITY, k++ % 4 == 0);
                assertFalse(Double.isNaN(trial.f1));
                assertEquals(root, tree.getRoot().getNumber());
                assertEquals(newick, TreeUtils.newick(tree));
                assertEquals(f, joint.getLogLikelihood(), 1E-10 * Math.abs(f));
                assertArrayEquals(parameters, values(evaluator.getRoundParameters()), 0.0);
                assertArrayEquals(heightGradient, gradient.getGradientLogDensity(), 1E-8);
            }
            assertTrue(evaluator.getVerifications() >= k);

            final List<Parameter> saved = new ArrayList<Parameter>(evaluator.getRoundParameters());
            TreeStateSnapshot start = new TreeStateSnapshot(tree, saved, joint);
            for (int m = 0; m < 3; m++) {
                List<RootChangingSprMove> moves = RootChangingSprMove.enumerate(tree, 3, null);
                TrialEvaluator.Trial trial = evaluator.trial(moves.get((7 * m + 2) % moves.size()),
                        Double.NEGATIVE_INFINITY, m == 1);
                final int before = tree.getRoot().getNumber();
                final long hash = TopologyHash.of(tree);
                final double accepted = evaluator.accept(trial);
                assertEquals(trial.best, accepted, 1E-8 * Math.abs(accepted));
                assertTrue(before != tree.getRoot().getNumber());
                assertTrue(hash != TopologyHash.of(tree));
                joint.makeDirty();
                assertEquals(accepted, joint.getLogLikelihood(), 1E-8 * Math.abs(accepted));
            }
            final int current = tree.getRoot().getNumber();
            start.restore(); // by clade when the root's node changed
            assertEquals(current, tree.getRoot().getNumber());
            assertEquals(TreeUtils.uniqueNewick(tree, tree.getRoot()), uniqueNewick(newick, objects));
            assertEquals(f, joint.getLogLikelihood(), 1E-8 * Math.abs(f));
        }
    }

    // the search with root moves: F never decreases row by row, root rounds run, and a constrained clade stays one
    public void testSearch() throws Exception {
        for (int[] clade : new int[][]{null, {0, 1, 2, 3, 4, 5, 6, 7}}) {
            Map<String, XMLObject> objects = parse(new Setup().relaxed(true).clade(clade)
                    .search("moves=\"nni,spr,root\"").xml());
            if (objects == null) {
                return;
            }
            TreeModel tree = tree(objects);
            TopologyConstraints constraints = new TopologyConstraints(tree, joint(objects), null);
            TimeTreeSearchTest.perturb(tree, 4, constraints);
            TimeTreeSearch search = TimeTreeSearchTest.search(objects);
            search.run();
            double previous = Double.NEGATIVE_INFINITY;
            for (String row : search.getLog().getRows()) {
                String[] fields = row.split("\t", -1);
                double after = Double.parseDouble(fields[8]);
                assertTrue(row, after >= Double.parseDouble(fields[7]) - 1E-6 * Math.abs(after));
                assertTrue(row, after >= previous - 1E-6 * Math.abs(after));
                previous = after;
            }
            assertTrue(search.getRootSprRounds() > 0);
            assertTrue(search.getReport().contains("RootSprRounds: "));
            assertEquals(TimeTreeSearch.StopReason.CONVERGED, search.getStopReason());
            assertTrue(constraints.areSatisfied(tree));
        }
    }

    // -------------------------------------------------------------------------------------------------------------

    // the rate above each node by node number (NaN at the root)
    private static double[] ratesByNode(TreeModel tree, ArbitraryBranchRates rates) {
        double[] result = new double[tree.getNodeCount()];
        for (int n = 0; n < result.length; n++) {
            NodeRef node = tree.getNode(n);
            result[n] = tree.isRoot(node) ? Double.NaN : rates.getBranchRate(tree, node);
        }
        return result;
    }

    // the canonical rooted newick of a newick string with heights, through a fresh tree of the same taxa
    private static String uniqueNewick(String newick, Map<String, XMLObject> objects) throws Exception {
        dr.evolution.io.NewickImporter importer = new dr.evolution.io.NewickImporter(newick);
        dr.evolution.tree.Tree tree = importer.importTree(null);
        return TreeUtils.uniqueNewick(tree, tree.getRoot());
    }

    private static double[] values(List<Parameter> parameters) {
        List<Double> values = new ArrayList<Double>();
        for (Parameter parameter : parameters) {
            for (double value : parameter.getParameterValues()) {
                values.add(value);
            }
        }
        double[] result = new double[values.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = values.get(i);
        }
        return result;
    }

    private static void assertValidTimeTree(TreeModel tree) {
        for (int n = 0; n < tree.getNodeCount(); n++) {
            NodeRef node = tree.getNode(n);
            if (!tree.isRoot(node)) {
                assertTrue(tree.getNodeHeight(node) < tree.getNodeHeight(tree.getParent(node)));
            }
        }
    }

    private static void assertArrayEquals(double[] expected, double[] actual, double tolerance) {
        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals("element " + i, expected[i], actual[i], tolerance * Math.max(1.0, Math.abs(expected[i])));
        }
    }
}
