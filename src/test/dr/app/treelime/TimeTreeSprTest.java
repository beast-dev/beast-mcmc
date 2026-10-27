/*
 * TimeTreeSprTest.java
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
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.search.ClampedSprMove;
import dr.evomodel.tree.search.SprNeighbourhood;
import dr.evomodel.tree.search.TopologyConstraints;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
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
import static test.dr.app.treelime.TimeTreeSearchTest.tree;

/**
 * Clamped-height SPR moves of the topology search (plan 7 M4) on the Dengue XML of TimeTreeSearchTest: V1 for SPR
 * (target sets against a brute-force filter for radii 1 to 5, with and without tip dates; the targets that keep P's
 * height are the lineages crossing it) and V8 for SPR (the constraint rule against the monophyly statistic over every
 * SPR within radius 5).
 */
public class TimeTreeSprTest extends TestCase {

    public TimeTreeSprTest(String name) {
        super(name);
    }

    // V1 (SPR): the breadth-first targets equal a brute-force filter of the pruned tree, for each radius window
    public void testTargets() throws Exception {
        for (boolean dated : new boolean[]{true, false}) {
            Map<String, XMLObject> objects = parse(new Setup().dated(dated).xml());
            if (objects == null) {
                return;
            }
            TreeModel tree = tree(objects);
            List<Integer> prunes = SprNeighbourhood.pruneOrder(tree);
            assertEquals(tree.getNodeCount() - 3, prunes.size()); // every node but the root and its children
            for (int i : prunes) {
                for (int rMin = 1; rMin <= 3; ++rMin) {
                    for (int rMax = rMin; rMax <= 5; ++rMax) {
                        Set<Integer> found = new HashSet<Integer>();
                        for (ClampedSprMove move : SprNeighbourhood.enumerate(tree, i, rMin, rMax, null)) {
                            assertTrue(move.describe(), move.isValid(tree));
                            assertTrue(move.describe() + " twice", found.add(move.getJ()));
                        }
                        assertEquals(dated + " i " + i + " [" + rMin + ", " + rMax + "]",
                                bruteForce(tree, i, rMin, rMax), found);
                    }
                }
                // with no limit, the targets that keep P's height are the lineages crossing it
                NodeRef p = tree.getParent(tree.getNode(i));
                double hP = tree.getNodeHeight(p);
                int kept = 0;
                for (ClampedSprMove move : SprNeighbourhood.enumerate(tree, i, 1, tree.getNodeCount(), null)) {
                    kept += move.needsRepair(tree) ? 0 : 1;
                }
                int crossing = 0;
                for (int n = 0; n < tree.getNodeCount(); ++n) {
                    NodeRef node = tree.getNode(n);
                    if (!tree.isRoot(node) && tree.getNodeHeight(node) < hP &&
                            tree.getNodeHeight(tree.getParent(node)) > hP) {
                        crossing++;
                    }
                }
                assertEquals(crossing, kept);
            }
        }
    }

    // every SPR keeps the tree a valid time tree, and its local block exists after the move
    public void testApply() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().xml());
        if (objects == null) {
            return;
        }
        TreeModel tree = tree(objects);
        Likelihood joint = joint(objects);
        Model model = joint.getModel();
        final String before = dr.evolution.tree.TreeUtils.uniqueNewick(tree, tree.getRoot());
        int repaired = 0;
        int count = 0;
        for (int i : SprNeighbourhood.pruneOrder(tree)) {
            for (ClampedSprMove move : SprNeighbourhood.enumerate(tree, i, 1, 3, null)) {
                joint.getLogLikelihood();
                final double hP = tree.getNodeHeight(tree.getNode(move.getP()));
                final boolean repair = move.needsRepair(tree);
                final double newHeight = move.newHeight(tree);
                model.storeModelState();
                move.apply(tree, null);
                count++;
                repaired += repair ? 1 : 0;
                assertEquals(move.getPj(), tree.getParent(tree.getNode(move.getP())).getNumber());
                assertEquals(move.getP(), tree.getParent(tree.getNode(move.getJ())).getNumber());
                assertEquals(move.getQ(), tree.getParent(tree.getNode(move.getC())).getNumber());
                assertEquals(repair ? newHeight : hP, tree.getNodeHeight(tree.getNode(move.getP())), 0.0);
                for (int n = 0; n < tree.getNodeCount(); ++n) {
                    NodeRef node = tree.getNode(n);
                    if (!tree.isRoot(node)) {
                        assertTrue(move.describe(), tree.getNodeHeight(node) < tree.getNodeHeight(tree.getParent(node)));
                    }
                }
                assertFalse(before.equals(dr.evolution.tree.TreeUtils.uniqueNewick(tree, tree.getRoot())));
                assertTrue(joint.getLogLikelihood() > Double.NEGATIVE_INFINITY);
                model.restoreModelState();
            }
        }
        assertEquals(before, dr.evolution.tree.TreeUtils.uniqueNewick(tree, tree.getRoot()));
        assertTrue(count > 0);
        assertTrue(repaired > 0 && repaired < count);
    }

    // V8 (SPR): the rule refuses exactly the SPRs after which the monophyly statistic fails, for a small and a large
    // clade, with constrained clades inside, around and apart from the pruned subtree
    public void testConstraintRule() throws Exception {
        for (int[] clade : new int[][]{{1, 2, 3}, {0, 1, 2, 3, 4, 5, 6, 7}}) {
            Map<String, XMLObject> objects = parse(new Setup().clade(clade).xml());
            if (objects == null) {
                return;
            }
            TreeModel tree = tree(objects);
            Likelihood joint = joint(objects);
            Model model = joint.getModel();
            TopologyConstraints constraints = new TopologyConstraints(tree, joint, null);
            assertEquals(1, constraints.size());
            BitSet[] clades = TopologyConstraints.clades(tree);
            BitSet constraint = null;
            for (BitSet bits : clades) {
                if (bits.cardinality() == clade.length && constraints.isConstrained(tree, indexOf(clades, bits))) {
                    constraint = bits;
                }
            }
            assertNotNull(constraint);
            int refused = 0, inside = 0, around = 0, apart = 0;
            for (int i : SprNeighbourhood.pruneOrder(tree)) {
                BitSet subtree = clades[i];
                List<ClampedSprMove> all = SprNeighbourhood.enumerate(tree, i, 1, 5, null);
                List<ClampedSprMove> allowed = SprNeighbourhood.enumerate(tree, i, 1, 5, constraints);
                for (ClampedSprMove move : all) {
                    final boolean violates = constraints.violates(tree, move, clades);
                    joint.getLogLikelihood();
                    model.storeModelState();
                    move.apply(tree, null);
                    assertEquals(move.describe(), violates, joint.getLogLikelihood() == Double.NEGATIVE_INFINITY);
                    model.restoreModelState();
                    refused += violates ? 1 : 0;
                    assertEquals(move.describe(), !violates, allowed.contains(move));
                    BitSet and = (BitSet) subtree.clone();
                    and.and(constraint);
                    if (and.equals(constraint)) {
                        inside++;
                    } else if (and.equals(subtree)) {
                        around++;
                    } else if (and.isEmpty()) {
                        apart++;
                    }
                }
            }
            assertTrue(refused > 0);
            assertTrue(inside > 0 && around > 0 && apart > 0);
        }
    }

    // the search with SPR: F never decreases row by row, SPR rounds run and accept moves from a perturbed start, the
    // result is at least what NNIs alone reach, and a constrained clade stays a clade
    public void testSearch() throws Exception {
        double[] finals = new double[2];
        for (int run = 0; run < 2; ++run) {
            Map<String, XMLObject> objects = parse(new Setup().clade(new int[]{1, 2, 3})
                    .search(run == 0 ? "moves=\"nni,spr\"" : "moves=\"nni\"").xml());
            if (objects == null) {
                return;
            }
            TreeModel tree = tree(objects);
            TopologyConstraints constraints = new TopologyConstraints(tree, joint(objects), null);
            TimeTreeSearchTest.perturb(tree, 5, constraints);
            dr.evomodel.tree.search.TimeTreeSearch search = TimeTreeSearchTest.search(objects);
            search.run();
            finals[run] = search.getFinalValue();
            double previous = Double.NEGATIVE_INFINITY;
            boolean sprRows = false;
            for (String row : search.getLog().getRows()) {
                String[] fields = row.split("\t", -1);
                double before = Double.parseDouble(fields[7]);
                double after = Double.parseDouble(fields[8]);
                assertTrue(row, after >= before - 1E-6 * Math.abs(before));
                assertTrue(row, after >= previous - 1E-6 * Math.abs(after));
                previous = after;
                sprRows |= fields[2].startsWith("spr_");
            }
            assertTrue(constraints.areSatisfied(tree));
            assertEquals(dr.evomodel.tree.search.TimeTreeSearch.StopReason.CONVERGED, search.getStopReason());
            if (run == 0) {
                assertTrue(sprRows);
                assertTrue(search.getSprRounds() > 0);
                assertTrue(search.getReport().contains("SprRounds: "));
            } else {
                assertEquals(0, search.getSprRounds());
            }
        }
        assertTrue(finals[0] + " < " + finals[1], finals[0] >= finals[1] - 1E-6 * Math.abs(finals[1]));
    }

    // V13 with SPR: the same XML and start give the same log, row for row (but the seconds), and the same tree
    public void testDeterminism() throws Exception {
        String xml = new Setup().relaxed(true).search("moves=\"nni,spr\" seed=\"5\"").xml();
        if (xml == null) {
            return;
        }
        List<String> logs = new ArrayList<String>();
        List<String> trees = new ArrayList<String>();
        for (int run = 0; run < 2; run++) {
            Map<String, XMLObject> objects = parse(xml);
            TimeTreeSearchTest.perturb(tree(objects), 4);
            dr.evomodel.tree.search.TimeTreeSearch search = TimeTreeSearchTest.search(objects);
            search.run();
            assertTrue(search.getSprRounds() > 0);
            StringBuilder sb = new StringBuilder();
            for (String row : search.getLog().getRows()) {
                sb.append(row, 0, row.lastIndexOf('\t')).append('\n');
            }
            logs.add(sb.toString());
            trees.add(dr.evolution.tree.TreeUtils.newick(tree(objects)));
        }
        assertEquals(logs.get(0), logs.get(1));
        assertEquals(trees.get(0), trees.get(1));
    }

    private static int indexOf(BitSet[] clades, BitSet bits) {
        for (int n = 0; n < clades.length; ++n) {
            if (clades[n] == bits) {
                return n;
            }
        }
        return -1;
    }

    // the targets of the subtree below i by a different route: the pruned tree as a parent map, branches adjacent
    // when they share a node, all distances by breadth-first search
    private static Set<Integer> bruteForce(TreeModel tree, int i, int rMin, int rMax) {
        NodeRef nodeI = tree.getNode(i);
        NodeRef p = tree.getParent(nodeI);
        NodeRef c = tree.getChild(p, 0) == nodeI ? tree.getChild(p, 1) : tree.getChild(p, 0);
        NodeRef q = tree.getParent(p);
        BitSet below = TopologyConstraints.clades(tree)[i];
        int n = tree.getNodeCount();
        int[] parent = new int[n];
        boolean[] present = new boolean[n];
        for (int x = 0; x < n; ++x) {
            NodeRef node = tree.getNode(x);
            boolean inSubtree = isBelow(tree, node, nodeI);
            present[x] = !inSubtree && node != p;
            parent[x] = tree.isRoot(node) ? -1 : (node == c ? q.getNumber() : tree.getParent(node).getNumber());
        }
        // branches: present nodes with a parent; distances from c's branch
        int[] distance = new int[n];
        java.util.Arrays.fill(distance, Integer.MAX_VALUE);
        distance[c.getNumber()] = 0;
        List<Integer> queue = new ArrayList<Integer>();
        queue.add(c.getNumber());
        for (int h = 0; h < queue.size(); ++h) {
            int x = queue.get(h);
            for (int y = 0; y < n; ++y) {
                if (present[y] && parent[y] >= 0 && distance[y] == Integer.MAX_VALUE &&
                        (y == parent[x] || parent[y] == x || parent[y] == parent[x])) {
                    distance[y] = distance[x] + 1;
                    queue.add(y);
                }
            }
        }
        Set<Integer> targets = new HashSet<Integer>();
        double hI = tree.getNodeHeight(nodeI);
        for (int y = 0; y < n; ++y) {
            if (present[y] && parent[y] >= 0 && y != c.getNumber() && distance[y] >= rMin && distance[y] <= rMax &&
                    tree.getNodeHeight(tree.getNode(parent[y])) > hI) {
                targets.add(y);
            }
        }
        assertTrue(below.get(nodeI.getNumber()) || !tree.isExternal(nodeI));
        return targets;
    }

    private static boolean isBelow(TreeModel tree, NodeRef node, NodeRef ancestor) {
        for (NodeRef x = node; x != null; x = tree.getParent(x)) {
            if (x == ancestor) {
                return true;
            }
        }
        return false;
    }
}
