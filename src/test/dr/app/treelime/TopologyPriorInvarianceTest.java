/*
 * TopologyPriorInvarianceTest.java
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
import dr.app.treelime.options.PartitionTreePrior;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.types.TreePriorType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.demographicmodel.DemographicModel;
import dr.evomodel.coalescent.TreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.search.TreeStateSnapshot;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * The tree priors TreeLime writes do not change under topology moves that keep every node height (plan 7 M0, V4),
 * so the search's cheap tiers compare topologies by the tree likelihood alone: coalescent models (through
 * TreeIntervals, which builds intervals from heights), skygrid and skyride, Yule and birth-death, and the whole-tree
 * CTMC reference prior (whose tree length, the sum of internal heights plus the root's minus the tips', depends on
 * heights only). Moves: feasible rooted NNIs and fixed-height subtree prune and regraft (FHSPR), each with every
 * height unchanged. A coalescent that excludes a subtree is unchanged as long as the subtree stays a clade.
 */
public class TopologyPriorInvarianceTest extends TestCase {

    private static final double TOLERANCE = 1E-10;

    public TopologyPriorInvarianceTest(String name) {
        super(name);
    }

    // dated tips with an estimated strict clock rate: coalescent priors plus the CTMC reference prior
    public void testCoalescentPriorsWithTipDates() throws Exception {
        for (TreePriorType type : new TreePriorType[]{TreePriorType.CONSTANT, TreePriorType.EXPONENTIAL,
                TreePriorType.SKYGRID, TreePriorType.GMRF_SKYRIDE}) {
            TreeLimeOptions options = TreeLimeHarness.importData("examples/Data/Dengue4.env.nex");
            if (options == null) {
                return;
            }
            TreeLimeHarness.setTipDatesFromNames(options);
            Map<String, XMLObject> objects = parse(options, type);
            assertTrue(type + ": a CTMC reference prior", hasCtmcScalePrior(objects));
            checkInvariance(type.toString(), objects);
        }
    }

    // contemporaneous tips: speciation priors as well
    public void testAllPriorsContemporaneous() throws Exception {
        for (TreePriorType type : new TreePriorType[]{TreePriorType.CONSTANT, TreePriorType.SKYGRID,
                TreePriorType.YULE, TreePriorType.BIRTH_DEATH}) {
            TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
            if (options == null) {
                return;
            }
            checkInvariance(type.toString(), parse(options, type));
        }
    }

    // a constant-size coalescent over the tree without one of its clades: unchanged by the moves that keep the clade
    public void testCoalescentExcludingASubtree() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData("examples/Data/Dengue4.env.nex");
        if (options == null) {
            return;
        }
        TreeLimeHarness.setTipDatesFromNames(options);
        Map<String, XMLObject> objects = parse(options, TreePriorType.CONSTANT);
        TreeModel tree = (TreeModel) objects.get("treeModel").getNativeObject();
        Taxa excluded = cladeOfSize(tree, 4);
        CoalescentLikelihood coalescent = new CoalescentLikelihood(new TreeIntervals(tree, null,
                Collections.<TaxonList>singletonList(excluded)),
                (DemographicModel) objects.get("constant").getNativeObject());
        final double before = coalescent.getLogLikelihood();
        final Set<String> clade = ids(excluded);

        int moves = 0;
        for (int k = 0; k < 40 && moves < 6; k++) {
            TreeStateSnapshot snapshot = new TreeStateSnapshot(tree, null, null);
            boolean moved = k % 2 == 0 ? TreeStateSnapshotTest.rootedNni(tree, k) : fixedHeightSpr(tree, k);
            if (!moved) {
                continue;
            }
            if (!TreeUtils.isMonophyletic(tree, clade)) {
                snapshot.restore(); // the search's constraint filter refuses such moves
                continue;
            }
            moves++;
            assertEquals("move " + k, before, coalescent.getLogLikelihood(), TOLERANCE * Math.abs(before));
        }
        assertTrue(moves >= 4);
    }

    // -------------------------------------------------------------------------------------------------------------

    private static Map<String, XMLObject> parse(TreeLimeOptions options, TreePriorType type) throws Exception {
        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            prior.setNodeHeightPrior(type);
        }
        options.maximizationOptions.printToScreen = false;
        String xml = TreeLimeHarness.generate(options);
        xml = xml.replace("<maximizeWrtParameter id=\"maximizer\"",
                "<maximizeWrtParameter id=\"maximizer\" maximizeOnParse=\"false\"");
        xml = xml.substring(0, xml.indexOf("<!-- The optimum")) + "</beast>\n";
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        return parser.parse(new StringReader(xml), false);
    }

    private static boolean hasCtmcScalePrior(Map<String, XMLObject> objects) {
        for (Likelihood term : terms(objects)) {
            if (term instanceof dr.evomodel.tree.CTMCScalePrior) {
                return true;
            }
        }
        return false;
    }

    private static List<Likelihood> terms(Map<String, XMLObject> objects) {
        return ((CompoundLikelihood) objects.get("prior").getNativeObject()).getLikelihoods();
    }

    // every term of the prior, after each of several NNIs and FHSPRs, equals its value before them
    private static void checkInvariance(String label, Map<String, XMLObject> objects) {
        TreeModel tree = (TreeModel) objects.get("treeModel").getNativeObject();
        List<Likelihood> terms = terms(objects);
        double[] before = values(terms);
        double[] heights = sortedHeights(tree);
        String start = TreeUtils.uniqueNewick(tree, tree.getRoot());

        int nni = 0;
        int spr = 0;
        for (int k = 0; k < 8; k++) {
            if (TreeStateSnapshotTest.rootedNni(tree, k)) {
                nni++;
                check(label + " after NNI " + k, terms, before);
            }
            if (fixedHeightSpr(tree, 2 * k + 1)) {
                spr++;
                check(label + " after FHSPR " + k, terms, before);
            }
        }
        assertTrue(label + ": " + nni + " NNIs and " + spr + " FHSPRs", nni >= 3 && spr >= 3);
        assertFalse(label + ": the topology changed", start.equals(TreeUtils.uniqueNewick(tree, tree.getRoot())));
        double[] after = sortedHeights(tree);
        for (int i = 0; i < heights.length; i++) {
            assertEquals(heights[i], after[i], 0.0);
        }
    }

    private static void check(String label, List<Likelihood> terms, double[] before) {
        double[] after = values(terms);
        for (int i = 0; i < before.length; i++) {
            Likelihood term = terms.get(i);
            String name = label + ": " + term.getClass().getSimpleName() + " " + term.getId();
            if (term instanceof dr.evomodel.tree.CTMCScalePrior) {
                // the tree length is a sum of the same heights in another order
                assertEquals(name, before[i], after[i], 1E-12 * Math.max(1.0, Math.abs(before[i])));
            } else {
                assertEquals(name, before[i], after[i], TOLERANCE * Math.max(1.0, Math.abs(before[i])));
            }
        }
    }

    private static double[] values(List<Likelihood> terms) {
        double[] values = new double[terms.size()];
        for (int i = 0; i < values.length; i++) {
            values[i] = terms.get(i).getLogLikelihood();
        }
        return values;
    }

    private static double[] sortedHeights(TreeModel tree) {
        double[] heights = new double[tree.getNodeCount()];
        for (int i = 0; i < heights.length; i++) {
            heights[i] = tree.getNodeHeight(tree.getNode(i));
        }
        java.util.Arrays.sort(heights);
        return heights;
    }

    /**
     * A fixed-height subtree prune and regraft: prune the subtree under i together with its parent P (not the
     * root), joining P's other child to P's parent, and regraft P onto an edge (j, parent of j) that spans P's height,
     * outside the subtree and not where P was. Takes the skip-th feasible move in node order.
     *
     * @return whether there was one
     */
    static boolean fixedHeightSpr(TreeModel tree, int skip) {
        List<NodeRef[]> moves = new ArrayList<NodeRef[]>();
        for (int a = 0; a < tree.getNodeCount(); a++) {
            NodeRef i = tree.getNode(a);
            if (tree.isRoot(i) || tree.isRoot(tree.getParent(i))) {
                continue;
            }
            NodeRef p = tree.getParent(i);
            NodeRef s = tree.getChild(p, 0) == i ? tree.getChild(p, 1) : tree.getChild(p, 0);
            final double h = tree.getNodeHeight(p);
            for (int b = 0; b < tree.getNodeCount(); b++) {
                NodeRef j = tree.getNode(b);
                if (tree.isRoot(j) || j == s || j == p || isDescendant(tree, j, i)) {
                    continue;
                }
                NodeRef k = tree.getParent(j);
                if (k == p) { // j is s, handled above, or i
                    continue;
                }
                if (tree.getNodeHeight(j) < h && h < tree.getNodeHeight(k)) {
                    moves.add(new NodeRef[]{i, p, s, j});
                }
            }
        }
        if (moves.size() <= skip) {
            return false;
        }
        NodeRef[] move = moves.get(skip);
        NodeRef p = move[1];
        NodeRef s = move[2];
        NodeRef j = move[3];
        NodeRef g = tree.getParent(p);
        tree.beginTreeEdit();
        tree.removeChild(g, p);
        tree.removeChild(p, s);
        tree.addChild(g, s);
        NodeRef k = tree.getParent(j);
        tree.removeChild(k, j);
        tree.addChild(p, j);
        tree.addChild(k, p);
        tree.endTreeEdit();
        return true;
    }

    private static boolean isDescendant(TreeModel tree, NodeRef node, NodeRef ancestor) {
        for (NodeRef n = node; n != null; n = tree.getParent(n)) {
            if (n == ancestor) {
                return true;
            }
        }
        return false;
    }

    // the tips of the first internal node, in node order, with this many tips below it
    private static Taxa cladeOfSize(TreeModel tree, int size) {
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getInternalNode(i);
            Set<String> tips = TreeUtils.getDescendantLeaves(tree, node);
            if (tips.size() == size) {
                Taxa taxa = new Taxa("excluded");
                for (String id : tips) {
                    taxa.addTaxon(tree.getTaxon(tree.getTaxonIndex(id)));
                }
                return taxa;
            }
        }
        throw new IllegalStateException("no clade of " + size + " tips");
    }

    private static Set<String> ids(Taxa taxa) {
        Set<String> ids = new HashSet<String>();
        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            ids.add(taxa.getTaxonId(i));
        }
        return ids;
    }
}
