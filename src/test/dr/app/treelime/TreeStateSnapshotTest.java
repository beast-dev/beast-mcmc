/*
 * TreeStateSnapshotTest.java
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
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.search.TreeStateSnapshot;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.util.Transform;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * The number-keyed tree snapshot of the topology search (plan 7 M0, V14), on the XML TreeLime writes for dated
 * Dengue sequences with the HMC relaxed clock and a constant-size coalescent: after accepted topology moves (and a
 * rejected trial in TreeModel's own store window), a restore brings back the topology with its node numbers, the
 * heights, the per-branch rates and the other parameters, so the node-height proxy, the ratio transform, the
 * coalescent's intervals and the joint density all equal those at the save and those of a fresh parse. Across a root
 * change (plan 7 M6, V14b) the restore is by clade and keeps the root's node.
 */
public class TreeStateSnapshotTest extends TestCase {

    private static String xml = null;

    private Map<String, XMLObject> objects;
    private TreeModel tree;
    private Likelihood joint;
    private Parameter rates;

    public TreeStateSnapshotTest(String name) {
        super(name);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
        if (xml == null) {
            TreeLimeOptions options = TreeLimeHarness.importData("examples/Data/Dengue4.env.nex");
            if (options == null) {
                return;
            }
            TreeLimeHarness.setTipDatesFromNames(options);
            for (PartitionClockModel clock : options.getPartitionClockModels()) {
                clock.setClockType(ClockType.HMC_CLOCK);
            }
            options.maximizationOptions.printToScreen = false;
            String generated = TreeLimeHarness.generate(options);
            generated = generated.replace("<maximizeWrtParameter id=\"maximizer\"",
                    "<maximizeWrtParameter id=\"maximizer\" maximizeOnParse=\"false\"");
            xml = generated.substring(0, generated.indexOf("<!-- The optimum")) + "</beast>\n";
        }
        objects = parse();
        tree = (TreeModel) objects.get("treeModel").getNativeObject();
        joint = (Likelihood) objects.get("joint").getNativeObject();
        rates = (Parameter) objects.get("branchRatesHMC.rates").getNativeObject();
        setDistinctRates(rates);
    }

    public void testRestoreAfterAcceptedMoves() throws Exception {
        if (xml == null) {
            return;
        }
        final Parameter proxy = (Parameter) objects.get("nodeHeights").getNativeObject();
        final Transform ratioTransform = (Transform) objects.get("nodeHeightTransform").getNativeObject();
        final Parameter optimized = (Parameter) objects.get("optimizedParameters").getNativeObject();
        final Likelihood coalescent = (Likelihood) objects.get("coalescent").getNativeObject();
        final Parameter kappa = (Parameter) objects.get("kappa").getNativeObject();

        final double f0 = joint.getLogLikelihood();
        final String newick0 = TreeUtils.newick(tree);
        final double[] heights0 = heights(tree);
        final double[] rates0 = rates.getParameterValues();
        final double[] proxy0 = proxy.getParameterValues();
        final double[] ratios0 = ratioTransform.transform(proxy0, 0, proxy0.length);
        final double coalescent0 = coalescent.getLogLikelihood();
        final int rootNumber = tree.getRoot().getNumber();

        TreeStateSnapshot snapshot = new TreeStateSnapshot(tree, Collections.singletonList(optimized), joint);
        assertEquals(f0, snapshot.getLogDensity(), 0.0);

        // a trial in TreeModel's own store window, rejected
        dr.inference.model.Model model = ((CompoundLikelihood) joint).getModel();
        model.storeModelState();
        assertTrue(rootedNni(tree, 0));
        joint.getLogLikelihood();
        model.restoreModelState();
        assertEquals(newick0, TreeUtils.newick(tree));
        assertEquals(f0, joint.getLogLikelihood(), 1E-8 * Math.abs(f0));

        // accepted moves: NNIs at three places, a height, some rates and kappa
        for (int k = 0; k < 3; k++) {
            assertTrue(rootedNni(tree, 2 * k));
            joint.getLogLikelihood();
        }
        NodeRef node = firstInternalBelowRoot(tree);
        tree.setNodeHeight(node, 0.5 * (tree.getNodeHeight(node) + tree.getNodeHeight(tree.getParent(node))));
        for (int i = 0; i < rates.getDimension(); i += 3) {
            rates.setParameterValue(i, 2.0 * rates.getParameterValue(i));
        }
        kappa.setParameterValue(0, 7.0);
        final double f1 = joint.getLogLikelihood();
        assertFalse(newick0.equals(TreeUtils.newick(tree)));
        assertTrue(Math.abs(f1 - f0) > 1E-3);

        snapshot.restore();

        // TreeModel's own restore swaps its node arrays, so node objects change; the root's number does not
        assertEquals(rootNumber, tree.getRoot().getNumber());
        assertEquals(newick0, TreeUtils.newick(tree));
        assertArrayEquals(heights0, heights(tree), 0.0);
        assertArrayEquals(rates0, rates.getParameterValues(), 0.0);
        assertArrayEquals(proxy0, proxy.getParameterValues(), 0.0);
        assertArrayEquals(ratios0, ratioTransform.transform(proxy.getParameterValues(), 0, proxy0.length), 1E-12);
        assertEquals(coalescent0, coalescent.getLogLikelihood(), 1E-10 * Math.abs(coalescent0));
        assertEquals(f0, joint.getLogLikelihood(), 1E-8 * Math.abs(f0));

        // a fresh parse of the same state
        Map<String, XMLObject> fresh = parse();
        TreeModel freshTree = (TreeModel) fresh.get("treeModel").getNativeObject();
        setDistinctRates((Parameter) fresh.get("branchRatesHMC.rates").getNativeObject());
        assertEquals(newick0, TreeUtils.newick(freshTree));
        assertEquals(f0, ((Likelihood) fresh.get("joint").getNativeObject()).getLogLikelihood(), 1E-8 * Math.abs(f0));
    }

    // V14b: across a root change, the restore is by clade: the root's node is kept, the topology, heights and rates
    // per clade are those at the save, and the node-height proxy, the ratio transform, the coalescent's intervals,
    // the density and a short maximization equal those of a fresh parse; restoring with the saved root's node puts
    // back the saved numbering itself
    public void testRestoreByCladeAcrossARootChange() throws Exception {
        if (xml == null) {
            return;
        }
        final Parameter proxy = (Parameter) objects.get("nodeHeights").getNativeObject();
        final Transform ratioTransform = (Transform) objects.get("nodeHeightTransform").getNativeObject();
        final dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter maximizer =
                (dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter) objects.get("maximizer").getNativeObject();
        final Parameter optimized = maximizer.getParameter(); // everything the maximization below changes
        final Likelihood coalescent = (Likelihood) objects.get("coalescent").getNativeObject();
        final Parameter kappa = (Parameter) objects.get("kappa").getNativeObject();
        final dr.evomodel.branchratemodel.ArbitraryBranchRates branchRates =
                (dr.evomodel.branchratemodel.ArbitraryBranchRates) objects.get("branchRates").getNativeObject();

        final double f0 = joint.getLogLikelihood();
        final String newick0 = TreeUtils.newick(tree);
        final Map<String, Double> byClade0 = byClade(tree, branchRates);
        final double coalescent0 = coalescent.getLogLikelihood();
        final double[] tips0 = tipHeights(tree);
        final int root0 = tree.getRoot().getNumber();
        final String rootSplit0 = rootSplit(tree);
        TreeStateSnapshot snapshot = new TreeStateSnapshot(tree, Collections.singletonList(optimized), joint);

        // a root-changing SPR (a subtree two levels down above the root), an NNI, rates and kappa
        int deep = -1;
        for (int n = 0; n < tree.getNodeCount() && deep < 0; n++) {
            NodeRef node = tree.getNode(n);
            if (!tree.isRoot(node) && !tree.isRoot(tree.getParent(node)) && !tree.isRoot(tree.getParent(tree.getParent(node)))) {
                deep = n;
            }
        }
        dr.evomodel.tree.search.RootChangingSprMove.aboveRoot(tree, deep).apply(tree, null);
        joint.getLogLikelihood();
        assertTrue(rootedNni(tree, 1));
        for (int i = 0; i < rates.getDimension(); i += 3) {
            rates.setParameterValue(i, 2.0 * rates.getParameterValue(i));
        }
        kappa.setParameterValue(0, 7.0);
        final double f1 = joint.getLogLikelihood();
        final int root1 = tree.getRoot().getNumber();
        assertTrue(root1 != root0);
        assertFalse(rootSplit0.equals(rootSplit(tree)));
        assertTrue(Math.abs(f1 - f0) > 1E-3);
        TreeStateSnapshot moved = new TreeStateSnapshot(tree, Collections.singletonList(optimized), joint);

        snapshot.restore();

        assertEquals(root1, tree.getRoot().getNumber()); // the root's node is kept
        assertEquals(rootSplit0, rootSplit(tree));
        assertEquals(byClade0, byClade(tree, branchRates));
        assertArrayEquals(tips0, tipHeights(tree), 0.0);
        assertEquals(coalescent0, coalescent.getLogLikelihood(), 1E-10 * Math.abs(coalescent0));
        assertEquals(f0, joint.getLogLikelihood(), 1E-8 * Math.abs(f0));

        // a fresh parse of the same state, under the parse's own numbering
        Map<String, XMLObject> fresh = parse();
        TreeModel freshTree = (TreeModel) fresh.get("treeModel").getNativeObject();
        setDistinctRates((Parameter) fresh.get("branchRatesHMC.rates").getNativeObject());
        assertEquals(root0, freshTree.getRoot().getNumber());
        assertEquals(byClade0, byClade(freshTree,
                (dr.evomodel.branchratemodel.ArbitraryBranchRates) fresh.get("branchRates").getNativeObject()));
        final Parameter freshProxy = (Parameter) fresh.get("nodeHeights").getNativeObject();
        assertArrayEquals(sorted(freshProxy.getParameterValues()), sorted(proxy.getParameterValues()), 0.0);
        final Transform freshRatios = (Transform) fresh.get("nodeHeightTransform").getNativeObject();
        assertArrayEquals(sorted(freshRatios.transform(freshProxy.getParameterValues(), 0, freshProxy.getDimension())),
                sorted(ratioTransform.transform(proxy.getParameterValues(), 0, proxy.getDimension())), 1E-12);
        // each proxy coordinate is a height of the tree it maps to
        assertProxyMapsHeights(tree, proxy);
        // a short maximization from both ends at the same point
        dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter freshMaximizer =
                (dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter) fresh.get("maximizer").getNativeObject();
        maximizer.maximize(new dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter.Budget(10).setQuiet(true));
        freshMaximizer.maximize(new dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter.Budget(10).setQuiet(true));
        final double fMax = joint.getLogLikelihood();
        assertTrue(fMax > f0);
        assertEquals(((Likelihood) fresh.get("joint").getNativeObject()).getLogLikelihood(), fMax, 1E-6 * Math.abs(fMax));
        assertEquals(TreeUtils.uniqueNewick(freshTree, freshTree.getRoot()), TreeUtils.uniqueNewick(tree, tree.getRoot()));

        // the moved state comes back under its own root number (number-keyed), then the saved one with its root's node
        moved.restore();
        assertEquals(root1, tree.getRoot().getNumber());
        assertEquals(f1, joint.getLogLikelihood(), 1E-8 * Math.abs(f1));
        snapshot.restoreByClade(root0);
        assertEquals(root0, tree.getRoot().getNumber());
        assertEquals(newick0, TreeUtils.newick(tree));
        assertEquals(byClade0, byClade(tree, branchRates));
        assertEquals(f0, joint.getLogLikelihood(), 1E-8 * Math.abs(f0));
    }

    // the height of each internal node's clade and the rate above it, keyed by the clade's taxa
    private static Map<String, Double> byClade(TreeModel tree, dr.evomodel.branchratemodel.ArbitraryBranchRates rates) {
        Map<String, Double> result = new java.util.TreeMap<String, Double>();
        for (int n = 0; n < tree.getNodeCount(); n++) {
            NodeRef node = tree.getNode(n);
            String key = new java.util.TreeSet<String>(TreeUtils.getDescendantLeaves(tree, node)).toString();
            if (!tree.isExternal(node)) {
                result.put("h" + key, tree.getNodeHeight(node));
            }
            if (!tree.isRoot(node)) {
                result.put("r" + key, rates.getBranchRate(tree, node));
            }
        }
        return result;
    }

    private static String rootSplit(TreeModel tree) {
        Set<String> left = new java.util.TreeSet<String>(TreeUtils.getDescendantLeaves(tree, tree.getChild(tree.getRoot(), 0)));
        Set<String> right = new java.util.TreeSet<String>(TreeUtils.getDescendantLeaves(tree, tree.getChild(tree.getRoot(), 1)));
        return left.size() < right.size() || (left.size() == right.size() && left.toString().compareTo(right.toString()) < 0) ?
                left + "|" + right : right + "|" + left;
    }

    private static double[] tipHeights(TreeModel tree) {
        double[] heights = new double[tree.getExternalNodeCount()];
        for (int i = 0; i < heights.length; i++) {
            heights[i] = tree.getNodeHeight(tree.getExternalNode(i));
        }
        return heights;
    }

    private static double[] sorted(double[] values) {
        double[] copy = values.clone();
        java.util.Arrays.sort(copy);
        return copy;
    }

    // every proxy coordinate reads the height of an internal node
    private static void assertProxyMapsHeights(TreeModel tree, Parameter proxy) {
        for (int i = 0; i < proxy.getDimension(); i++) {
            final double value = proxy.getParameterValue(i);
            int matches = 0;
            for (int n = 0; n < tree.getInternalNodeCount(); n++) {
                matches += tree.getNodeHeight(tree.getInternalNode(n)) == value ? 1 : 0;
            }
            assertTrue("proxy " + i, matches >= 1);
        }
    }

    // -------------------------------------------------------------------------------------------------------------

    private static Map<String, XMLObject> parse() throws Exception {
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        return parser.parse(new StringReader(xml), false);
    }

    // distinct rates, so that a rate at the wrong branch changes the density
    private static void setDistinctRates(Parameter rates) {
        for (int i = 0; i < rates.getDimension(); i++) {
            rates.setParameterValue(i, 0.001 * (1.0 + 0.05 * i));
        }
    }

    private static double[] heights(TreeModel tree) {
        double[] heights = new double[tree.getNodeCount()];
        for (int i = 0; i < heights.length; i++) {
            heights[i] = tree.getNodeHeight(tree.getNode(i));
        }
        return heights;
    }

    private static NodeRef firstInternalBelowRoot(TreeModel tree) {
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getInternalNode(i);
            if (!tree.isRoot(node)) {
                return node;
            }
        }
        throw new IllegalStateException();
    }

    /**
     * A feasible rooted NNI at fixed heights, as the search makes them: at a non-root internal node P with parent G
     * and sibling U, swap a child a of P with U, where U is younger than P. Takes the skip-th feasible one in node
     * order.
     *
     * @return whether there was one
     */
    static boolean rootedNni(TreeModel tree, int skip) {
        List<NodeRef[]> moves = new ArrayList<NodeRef[]>();
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef p = tree.getInternalNode(i);
            if (tree.isRoot(p)) {
                continue;
            }
            NodeRef g = tree.getParent(p);
            NodeRef u = tree.getChild(g, 0) == p ? tree.getChild(g, 1) : tree.getChild(g, 0);
            if (tree.getNodeHeight(u) < tree.getNodeHeight(p)) {
                moves.add(new NodeRef[]{p, g, tree.getChild(p, 0), u});
            }
        }
        if (moves.size() <= skip) {
            return false;
        }
        NodeRef[] move = moves.get(skip);
        tree.beginTreeEdit();
        tree.removeChild(move[0], move[2]);
        tree.removeChild(move[1], move[3]);
        tree.addChild(move[0], move[3]);
        tree.addChild(move[1], move[2]);
        tree.endTreeEdit();
        return true;
    }

    private static void assertArrayEquals(double[] expected, double[] actual, double tolerance) {
        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals("element " + i, expected[i], actual[i], tolerance * Math.max(1.0, Math.abs(expected[i])));
        }
    }
}
