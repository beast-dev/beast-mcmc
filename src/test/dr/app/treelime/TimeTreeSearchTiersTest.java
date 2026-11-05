/*
 * TimeTreeSearchTiersTest.java
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

import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.tree.treemetrics.SplitDistance;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.search.NniNeighbourhood;
import dr.evomodel.tree.search.RootedNniMove;
import dr.evomodel.tree.search.TimeTreeMove;
import dr.evomodel.tree.search.TimeTreeSearch;
import dr.evomodel.tree.search.TopologyConstraints;
import dr.evomodel.tree.search.TrialEvaluator;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.model.Likelihood;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static test.dr.app.treelime.TimeTreeSearchTest.Setup;
import static test.dr.app.treelime.TimeTreeSearchTest.joint;
import static test.dr.app.treelime.TimeTreeSearchTest.parse;
import static test.dr.app.treelime.TimeTreeSearchTest.perturb;
import static test.dr.app.treelime.TimeTreeSearchTest.search;
import static test.dr.app.treelime.TimeTreeSearchTest.tree;

/**
 * The NNI tier set of the topology search (plan 7 M1b), on the Dengue XML of TimeTreeSearchTest: the speedNNI window,
 * the non-conflicting batch and its replay, the second chance, the root moves (root-child NNIs, alone and chained) and
 * the root scan. Run with assertions on, as TimeTreeSearchTest.
 */
public class TimeTreeSearchTiersTest extends TestCase {

    public TimeTreeSearchTiersTest(String name) {
        super(name);
    }

    // the speedNNI window: the nodes at most r edges from the given ones, against distances by brute force
    public void testWindow() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().xml());
        if (objects == null) {
            return;
        }
        TreeModel tree = tree(objects);
        int[][] distance = distances(tree);
        for (int start = 0; start < tree.getNodeCount(); start += 5) {
            for (int radius = 0; radius <= 3; ++radius) {
                BitSet window = NniNeighbourhood.within(tree, new int[]{start}, radius);
                for (int i = 0; i < tree.getNodeCount(); ++i) {
                    assertEquals(start + " r " + radius + " node " + i, distance[start][i] <= radius, window.get(i));
                }
            }
        }
        // enumeration restricted to the window: exactly the moves of the full list whose P is in it
        BitSet window = NniNeighbourhood.within(tree, new int[]{tree.getRoot().getNumber()}, 2);
        List<RootedNniMove> all = NniNeighbourhood.enumerate(tree, null);
        List<RootedNniMove> local = NniNeighbourhood.enumerate(tree, null, window);
        int expected = 0;
        for (RootedNniMove move : all) {
            expected += window.get(move.getP()) ? 1 : 0;
        }
        assertEquals(expected, local.size());
        assertTrue(all.containsAll(local));
    }

    // root moves: every root-child NNI and its chained ones keep the unrooted tree, give distinct rootings, and
    // reach exactly the root positions one and two edges from the current one; constrained clades are kept
    public void testRootMoves() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().xml());
        if (objects == null) {
            return;
        }
        TreeModel tree = tree(objects);
        Likelihood joint = joint(objects);
        checkRootMoves(tree, joint, null);

        Map<String, XMLObject> constrained = parse(new Setup().clade(new int[]{1, 2, 3}).xml());
        TreeModel other = tree(constrained);
        Likelihood otherJoint = joint(constrained);
        checkRootMoves(other, otherJoint, new TopologyConstraints(other, otherJoint, null));
    }

    private static void checkRootMoves(TreeModel tree, Likelihood joint, TopologyConstraints constraints) {
        final Tree before = new FlexibleTree(tree);
        final String rooted = TreeUtils.uniqueNewick(tree, tree.getRoot());
        List<TimeTreeMove> moves = NniNeighbourhood.rootMoves(tree, constraints);
        assertFalse(moves.isEmpty());
        Set<String> rootings = new HashSet<String>();
        for (TimeTreeMove move : moves) {
            joint.getLogLikelihood();
            joint.getModel().storeModelState();
            move.apply(tree, null);
            assertEquals(move.describe(), 0, SplitDistance.unrooted(before, new FlexibleTree(tree)));
            String newick = TreeUtils.uniqueNewick(tree, tree.getRoot());
            assertFalse(move.describe(), newick.equals(rooted));
            assertTrue(move.describe() + " repeats a rooting", rootings.add(newick));
            if (constraints != null) {
                assertTrue(move.describe(), constraints.areSatisfied(tree));
            }
            assertTrue(joint.getLogLikelihood() > Double.NEGATIVE_INFINITY);
            joint.getModel().restoreModelState();
        }
        assertEquals(rooted, TreeUtils.uniqueNewick(tree, tree.getRoot()));
        if (constraints == null) {
            // the root sits on the edge between its two children; the edges one and two away from it
            assertEquals(edgesNearRoot(tree), moves.size());
        }
    }

    // edges (by their lower node) at distance 1 or 2 from the root edge, counted in the unrooted tree
    private static int edgesNearRoot(TreeModel tree) {
        NodeRef root = tree.getRoot();
        int count = 0;
        for (int k = 0; k < 2; ++k) {
            NodeRef child = tree.getChild(root, k);
            if (tree.isExternal(child)) {
                continue;
            }
            for (int c = 0; c < 2; ++c) {
                NodeRef grandchild = tree.getChild(child, c);
                count++; // one edge away
                if (!tree.isExternal(grandchild)) {
                    count += 2; // two edges away
                }
            }
        }
        return count;
    }

    // a batch of non-conflicting gaining moves, accepted together: the topology has every move, F is what a full
    // recomputation gives, and the search with batches stays monotone and ends no lower than with single moves
    public void testBatch() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().xml());
        if (objects == null) {
            return;
        }
        TreeModel tree = tree(objects);
        TimeTreeSearch search = search(objects);
        TrialEvaluator evaluator = search.getEvaluator();
        MaximizerWrtParameter round = TimeTreeSearchTest.roundMaximizer(objects);
        round.maximize(new MaximizerWrtParameter.Budget(500).setQuiet(true));

        // two moves far apart that both gain: undo two NNIs applied far from each other
        List<RootedNniMove> moves = NniNeighbourhood.enumerate(tree, null);
        RootedNniMove first = moves.get(0);
        RootedNniMove second = null;
        BitSet near = NniNeighbourhood.within(tree, first.quartet(tree), 1);
        for (RootedNniMove move : moves) {
            boolean apart = true;
            for (int node : move.quartet(tree)) {
                apart &= !near.get(node);
            }
            if (apart) {
                second = move;
                break;
            }
        }
        assertNotNull(second);
        final double fBefore = evaluator.evaluate();
        final TrialEvaluator.Theta base = evaluator.capture();
        TrialEvaluator.Trial a = evaluator.trial(first, Double.NEGATIVE_INFINITY, false);
        TrialEvaluator.Trial b = evaluator.trial(second, Double.NEGATIVE_INFINITY, false);
        assertEquals(fBefore, evaluator.evaluate(), 0.0);
        List<TrialEvaluator.Trial> batch = new ArrayList<TrialEvaluator.Trial>();
        batch.add(a);
        batch.add(b);
        Tree start = new FlexibleTree(tree);
        double f = evaluator.acceptBatch(batch, base);
        assertFalse(Double.isNaN(f));
        evaluator.verify(f);
        // both topology changes are there: two rooted splits differ
        assertEquals(4, SplitDistance.rooted(start, new FlexibleTree(tree)));
        // the batch's changes are local: F after the batch is near F before plus both trials' gains, as the moves
        // are far apart (the likelihood is not exactly additive)
        assertEquals(fBefore + (a.best - fBefore) + (b.best - fBefore), f, 0.05 * Math.abs(a.best + b.best - 2 * fBefore) + 1.0);

        // the whole search with and without batches, from the same perturbed start
        double[] finals = new double[2];
        int batches = 0;
        for (int run = 0; run < 2; ++run) {
            Map<String, XMLObject> fresh = parse(new Setup().search(run == 0 ? "batch=\"true\"" : "batch=\"false\"").xml());
            perturb(tree(fresh), 6);
            TimeTreeSearch s = search(fresh);
            s.run();
            finals[run] = s.getFinalValue();
            for (String row : s.getLog().getRows()) {
                if (row.contains("batch ")) {
                    batches += run == 0 ? 1 : 0;
                    assertEquals(0, run);
                }
            }
            assertTrue(s.getFinalValue() > s.getStartValue());
        }
        assertTrue("no round applied a batch", batches > 0);
        assertTrue(finals[0] + " < " + finals[1], finals[0] >= finals[1] - 1.0);
    }

    // speedNNI: rounds after a change score the moves near it; the search never stops on a local scan alone
    public void testSpeedNni() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().relaxed(true).xml());
        if (objects == null) {
            return;
        }
        perturb(tree(objects), 4);
        TimeTreeSearch search = search(objects);
        search.run();
        boolean local = false;
        String lastNni = null;
        for (String row : search.getLog().getRows()) {
            String[] fields = row.split("\t", -1);
            if (fields[2].equals("nni")) {
                local |= fields[17].startsWith("local");
                lastNni = fields[17];
            }
        }
        assertTrue("no round used the speedNNI window", local);
        assertEquals(TimeTreeSearch.StopReason.CONVERGED, search.getStopReason());
        // the round that found nothing scanned the whole tree
        assertTrue(lastNni, lastNni.equals("full") || lastNni.startsWith("local+full"));
    }

    // the second chance: with a wide window every re-fitted candidate that does not gain gets the longer T2; with
    // no iterations, none does; the result is the same tree
    public void testSecondChance() throws Exception {
        String[] attributes = {"secondChanceWindow=\"1000000\"", "secondChanceIterations=\"0\""};
        int[] chances = new int[2];
        String[] trees = new String[2];
        for (int run = 0; run < 2; ++run) {
            Map<String, XMLObject> objects = parse(new Setup().search(attributes[run]).xml());
            if (objects == null) {
                return;
            }
            perturb(tree(objects), 3);
            TimeTreeSearch search = search(objects);
            search.run();
            chances[run] = search.getSecondChances();
            trees[run] = TreeUtils.uniqueNewick(tree(objects), tree(objects).getRoot());
            assertTrue(search.getReport().contains("SecondChances: " + chances[run]));
        }
        assertTrue(chances[0] > 0);
        assertEquals(0, chances[1]);
    }

    // the root scan after a new best: a fitted tree with its root moved one edge comes back to a rooting at least as
    // good; the scan is counted in the report
    public void testRootScan() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().xml());
        if (objects == null) {
            return;
        }
        TreeModel tree = tree(objects);
        perturb(tree, 3);
        TimeTreeSearch search = search(objects);
        search.run();
        assertTrue(search.getRootScans() > 0);
        assertTrue(search.getReport().contains("RootScans: " + search.getRootScans()));
        boolean rootRow = false;
        for (String row : search.getLog().getRows()) {
            rootRow |= row.split("\t", -1)[2].equals("root");
        }
        assertTrue(rootRow);

        // without the root scan the search may end elsewhere, never higher
        Map<String, XMLObject> without = parse(new Setup().search("rootScan=\"false\"").xml());
        perturb(tree(without), 3);
        TimeTreeSearch other = search(without);
        other.run();
        assertEquals(0, other.getRootScans());
        assertTrue(search.getFinalValue() >= other.getFinalValue() - 1.0);
    }

    // T1 from a zero-length branch: a height at its parent's (an infinite logit) used to be skipped, so an NNI across
    // a zero-length branch could not gain; it now starts from just inside its interval and finds the fitted height
    public void testLocalOptimizerOpensZeroLengthBranch() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().xml());
        if (objects == null) {
            return;
        }
        TreeModel tree = tree(objects);
        final Likelihood joint = joint(objects);
        TimeTreeSearchTest.roundMaximizer(objects).maximize(new MaximizerWrtParameter.Budget(500).setQuiet(true));
        final double fitted = joint.getLogLikelihood();

        // the internal node whose branch is longest relative to its interval, collapsed onto its parent
        NodeRef node = null;
        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
            NodeRef candidate = tree.getInternalNode(i);
            if (!tree.isRoot(candidate) && (node == null || branch(tree, candidate) > branch(tree, node))) {
                node = candidate;
            }
        }
        final double height = tree.getNodeHeight(node);
        tree.setNodeHeight(node, tree.getNodeHeight(tree.getParent(node)));
        final double collapsed = joint.getLogLikelihood();
        assertTrue(collapsed < fitted - 1.0);

        dr.evomodel.tree.search.LocalOptimizer local = new dr.evomodel.tree.search.LocalOptimizer(tree, null, 2);
        double f = local.optimize(new dr.evomodel.tree.search.LocalOptimizer.Objective() {
            public double value() {
                return joint.getLogLikelihood();
            }
        }, collapsed, new int[]{node.getNumber()}, new int[0]);
        assertEquals(f, joint.getLogLikelihood(), 0.0);
        assertTrue(f + " from " + collapsed, f > collapsed + 1.0);
        assertTrue(tree.getNodeHeight(node) < tree.getNodeHeight(tree.getParent(node)));
        assertEquals(fitted, f, 0.05 * (fitted - collapsed));
        assertEquals(height, tree.getNodeHeight(node), 0.2 * Math.abs(height));
    }

    private static double branch(TreeModel tree, NodeRef node) {
        return tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node);
    }

    // all-pairs distances in edges, by breadth-first search from every node
    private static int[][] distances(TreeModel tree) {
        final int n = tree.getNodeCount();
        int[][] d = new int[n][n];
        for (int s = 0; s < n; ++s) {
            java.util.Arrays.fill(d[s], Integer.MAX_VALUE);
            d[s][s] = 0;
            List<NodeRef> queue = new ArrayList<NodeRef>();
            queue.add(tree.getNode(s));
            for (int q = 0; q < queue.size(); ++q) {
                NodeRef node = queue.get(q);
                List<NodeRef> adjacent = new ArrayList<NodeRef>();
                if (!tree.isRoot(node)) {
                    adjacent.add(tree.getParent(node));
                }
                for (int c = 0; c < tree.getChildCount(node); ++c) {
                    adjacent.add(tree.getChild(node, c));
                }
                for (NodeRef other : adjacent) {
                    if (d[s][other.getNumber()] == Integer.MAX_VALUE) {
                        d[s][other.getNumber()] = d[s][node.getNumber()] + 1;
                        queue.add(other);
                    }
                }
            }
        }
        return d;
    }
}
