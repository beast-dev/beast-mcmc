/*
 * TimeTreeInitializerTest.java
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
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.search.TimeTreeInitializer;
import dr.evomodel.tree.search.TimeTreeSearch;
import dr.evomodel.tree.search.TopologyConstraints;
import dr.inference.model.Likelihood;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.util.Collections;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static test.dr.app.treelime.TimeTreeSearchTest.Setup;
import static test.dr.app.treelime.TimeTreeSearchTest.joint;
import static test.dr.app.treelime.TimeTreeSearchTest.parse;
import static test.dr.app.treelime.TimeTreeSearchTest.rates;
import static test.dr.app.treelime.TimeTreeSearchTest.tree;

/**
 * The NJ start of the topology search (plan 7 M6, §7.1) on the Dengue XML of TimeTreeSearchTest: a clock-like donor
 * is rooted and dated back exactly by root-to-tip regression; the NJ tree of the alignment is adopted by clade with
 * the root's node kept, tips at their dates and a valid time tree; a donor that breaks a constraint is refused; and
 * the search runs from it.
 */
public class TimeTreeInitializerTest extends TestCase {

    public TimeTreeInitializerTest(String name) {
        super(name);
    }

    // substitutions proportional to time on the dated tree: the regression finds its root, its rate and its heights
    public void testClockLikeDonor() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().xml());
        if (objects == null) {
            return;
        }
        TreeModel tree = tree(objects);
        final double rate = 0.004;
        final String rooted = TreeUtils.uniqueNewick(tree, tree.getRoot());
        final double[] heights = sortedInternalHeights(tree);
        FlexibleTree donor = new FlexibleTree(tree);
        for (int n = 0; n < donor.getNodeCount(); n++) {
            donor.setNodeHeight(donor.getNode(n), rate * donor.getNodeHeight(donor.getNode(n)));
        }
        perturb(tree);
        assertFalse(rooted.equals(TreeUtils.uniqueNewick(tree, tree.getRoot())));
        final int root = tree.getRoot().getNumber();

        TimeTreeInitializer initializer = new TimeTreeInitializer(tree, null);
        assertTrue(initializer.adopt(donor, null));
        assertTrue(initializer.isDated());
        assertEquals(rate, initializer.getRate(), 1E-9 * rate);
        assertEquals(root, tree.getRoot().getNumber());
        assertEquals(rooted, TreeUtils.uniqueNewick(tree, tree.getRoot()));
        double[] adopted = sortedInternalHeights(tree);
        for (int i = 0; i < heights.length; i++) {
            assertEquals(heights[i], adopted[i], 1E-8 * heights[heights.length - 1]);
        }
        assertTrue(joint(objects).getLogLikelihood() > Double.NEGATIVE_INFINITY);
    }

    // contemporaneous tips: the root with the least variance of root-to-tip distances, here the clock tree's own, and the
    // rate that keeps the root height
    public void testContemporaneousDonor() throws Exception {
        TreeModel tree = new dr.evomodel.tree.DefaultTreeModel(new dr.evolution.io.NewickImporter(
                "((((A:1,B:1):1,C:2):1,(D:1.5,E:1.5):1.5):1,(F:2.5,G:2.5):1.5);").importTree(null));
        final String rooted = TreeUtils.uniqueNewick(tree, tree.getRoot());
        final double[] heights = sortedInternalHeights(tree);
        FlexibleTree donor = new FlexibleTree(tree);
        for (int n = 0; n < donor.getNodeCount(); n++) {
            donor.setNodeHeight(donor.getNode(n), 0.01 * donor.getNodeHeight(donor.getNode(n)));
        }
        TimeTreeSearchTest.perturb(tree, 2);
        assertFalse(rooted.equals(TreeUtils.uniqueNewick(tree, tree.getRoot())));
        assertEquals(4.0, tree.getNodeHeight(tree.getRoot()), 0.0);
        TimeTreeInitializer initializer = new TimeTreeInitializer(tree, null);
        assertTrue(initializer.adopt(donor, null));
        assertFalse(initializer.isDated());
        assertEquals(0.01, initializer.getRate(), 1E-12);
        assertEquals(rooted, TreeUtils.uniqueNewick(tree, tree.getRoot()));
        double[] adopted = sortedInternalHeights(tree);
        for (int i = 0; i < heights.length; i++) {
            assertEquals(heights[i], adopted[i], 1E-10);
        }
    }

    // the NJ tree of the alignment: adopted with the root's node kept, the tips at their dates, a valid time tree and
    // per-branch rates within range
    public void testNjStart() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().relaxed(true).xml());
        if (objects == null) {
            return;
        }
        TreeModel tree = tree(objects);
        Likelihood joint = joint(objects);
        final int root = tree.getRoot().getNumber();
        final double[] tips = new double[tree.getExternalNodeCount()];
        for (int t = 0; t < tips.length; t++) {
            tips[t] = tree.getNodeHeight(tree.getExternalNode(t));
        }
        TimeTreeInitializer initializer = new TimeTreeInitializer(tree, rates(objects));
        assertTrue(initializer.adopt(TimeTreeInitializer.patterns(tree, joint), null));
        assertTrue(initializer.isDated());
        assertTrue(initializer.getRate() > 0);
        assertEquals(root, tree.getRoot().getNumber());
        for (int t = 0; t < tips.length; t++) {
            assertEquals(tips[t], tree.getNodeHeight(tree.getExternalNode(t)), 0.0);
        }
        for (int n = 0; n < tree.getNodeCount(); n++) {
            NodeRef node = tree.getNode(n);
            if (!tree.isRoot(node)) {
                assertTrue(tree.getNodeHeight(node) < tree.getNodeHeight(tree.getParent(node)));
                final double r = rates(objects).getBranchRate(tree, node);
                assertTrue(r >= initializer.getRate() / 100.0 * (1 - 1E-9) &&
                        r <= initializer.getRate() * 100.0 * (1 + 1E-9));
            }
        }
        final double f = joint.getLogLikelihood();
        assertTrue(f > Double.NEGATIVE_INFINITY);
        joint.makeDirty();
        assertEquals(f, joint.getLogLikelihood(), 1E-8 * Math.abs(f));
    }

    // a donor without a constrained clade is refused, and the tree is left as it was
    public void testConstraintRefusal() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().xml());
        if (objects == null) {
            return;
        }
        TreeModel tree = tree(objects);
        FlexibleTree donor = new FlexibleTree(tree);
        perturb(tree);
        final String before = TreeUtils.newick(tree);
        // a clade of the perturbed tree that the donor does not have
        Set<String> taxa = null;
        for (int n = 0; n < tree.getNodeCount() && taxa == null; n++) {
            NodeRef node = tree.getNode(n);
            Set<String> clade = new HashSet<String>(TreeUtils.getDescendantLeaves(tree, node));
            if (clade.size() >= 2 && clade.size() < tree.getExternalNodeCount() - 1
                    && !TreeUtils.isMonophyletic(donor, clade)) {
                taxa = clade;
            }
        }
        assertNotNull(taxa);
        TopologyConstraints constraints = new TopologyConstraints(tree, null,
                Collections.<Set<String>>singletonList(taxa));
        assertFalse(new TimeTreeInitializer(tree, null).adopt(donor, constraints));
        assertEquals(before, TreeUtils.newick(tree));
    }

    // the search from the NJ tree
    public void testSearchFromNj() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().relaxed(true).search("startingTree=\"nj\" moves=\"nni,spr,root\"").xml());
        if (objects == null) {
            return;
        }
        TimeTreeSearch search = TimeTreeSearchTest.search(objects);
        search.run();
        assertEquals(TimeTreeSearch.StopReason.CONVERGED, search.getStopReason());
        assertTrue(search.getReport().contains("startingTree = nj"));
    }

    // a climb from the tree as it is, then one from the NJ tree: at least as good as either climb alone (the NJ climb
    // starts from the parameters as they were, not the first climb's fit, and re-fits on its own gains)
    public void testSearchFromBothStarts() throws Exception {
        Map<String, XMLObject> objects = parse(new Setup().relaxed(true).search("moves=\"nni,spr,root\"").xml());
        if (objects == null) {
            return;
        }
        TimeTreeSearch current = TimeTreeSearchTest.search(objects);
        current.run();
        final double fCurrent = current.getBestValue();

        objects = parse(new Setup().relaxed(true).search("startingTree=\"nj\" moves=\"nni,spr,root\"").xml());
        TimeTreeSearch nj = TimeTreeSearchTest.search(objects);
        nj.run();
        final double fNj = nj.getBestValue();

        objects = parse(new Setup().relaxed(true).search("startingTree=\"both\" moves=\"nni,spr,root\"").xml());
        TimeTreeSearch both = TimeTreeSearchTest.search(objects);
        both.run();
        assertTrue(both.getStopReason() == TimeTreeSearch.StopReason.CONVERGED ||
                both.getStopReason() == TimeTreeSearch.StopReason.MAX_ROUNDS);
        assertTrue(both.getReport().contains("startingTree = both"));
        assertTrue(both.getBestValue() >= fCurrent - 1E-6 * Math.abs(fCurrent));
        assertTrue(both.getBestValue() + " vs NJ alone " + fNj, both.getBestValue() >= fNj - 1E-6 * Math.abs(fNj));
        assertTrue(both.getFinalValue() >= both.getBestValue() - 1E-6 * Math.abs(fCurrent));
    }

    // -------------------------------------------------------------------------------------------------------------

    private static void perturb(TreeModel tree) {
        TimeTreeSearchTest.perturb(tree, 3);
    }

    private static double[] sortedInternalHeights(TreeModel tree) {
        double[] heights = new double[tree.getInternalNodeCount()];
        for (int i = 0; i < heights.length; i++) {
            heights[i] = tree.getNodeHeight(tree.getInternalNode(i));
        }
        java.util.Arrays.sort(heights);
        return heights;
    }
}
