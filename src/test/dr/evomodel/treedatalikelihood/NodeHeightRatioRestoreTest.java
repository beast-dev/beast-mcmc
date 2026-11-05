/*
 * NodeHeightRatioRestoreTest.java
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

package test.dr.evomodel.treedatalikelihood;

import dr.app.beast.BeastParser;
import dr.app.beast.BeastVersion;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightToRatiosTransformDelegate;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightTransform;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.StringReader;
import java.util.Arrays;
import java.util.Map;

/**
 * The ratio transform's cached ratios across a store and restore (plan 7 M0, V5). A node-height change marks the
 * ratios stale without recomputing them; a store then saves the stale vector, and the restore used to mark it
 * current. Readers that do not go through transform() first (the masks, getRatios, reports) then saw ratios of an
 * older tree. Every tip is at 0, so each ratio is a node's height over its parent's.
 */
public class NodeHeightRatioRestoreTest extends TestCase {

    private static final String XML = "<beast>\n" +
            "  <taxa id=\"taxa\"><taxon id=\"a\"/><taxon id=\"b\"/><taxon id=\"c\"/><taxon id=\"d\"/>" +
            "<taxon id=\"e\"/></taxa>\n" +
            "  <newick id=\"startingTree\">((((a:1,b:1):1,c:2):1,d:3):2,e:5);</newick>\n" +
            "  <treeModel id=\"treeModel\">\n" +
            "    <tree idref=\"startingTree\"/>\n" +
            "    <rootHeight><parameter id=\"treeModel.rootHeight\"/></rootHeight>\n" +
            "    <nodeHeights internalNodes=\"true\"><parameter id=\"treeModel.internalNodeHeights\"/></nodeHeights>\n" +
            "  </treeModel>\n" +
            "  <strictClockBranchRates id=\"branchRates\"><rate><parameter id=\"clock.rate\" value=\"1.0\"/></rate>" +
            "</strictClockBranchRates>\n" +
            "  <nodeHeightTransform id=\"nodeHeightTransform\" withRoot=\"false\">\n" +
            "    <treeModel idref=\"treeModel\"/><strictClockBranchRates idref=\"branchRates\"/>\n" +
            "    <ratios><parameter id=\"ratios\"/></ratios>\n" +
            "  </nodeHeightTransform>\n" +
            "</beast>\n";

    private TreeModel tree;
    private NodeHeightToRatiosTransformDelegate delegate;

    public NodeHeightRatioRestoreTest(String name) {
        super(name);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        Map<String, XMLObject> objects = parser.parse(new StringReader(XML), false);
        tree = (TreeModel) objects.get("treeModel").getNativeObject();
        delegate = (NodeHeightToRatiosTransformDelegate)
                ((NodeHeightTransform) objects.get("nodeHeightTransform").getNativeObject())
                        .getNodeHeightTransformDelegate();
    }

    // node-height events since the last read, a store, more events, a restore: the readers see the ratios of the
    // restored heights
    public void testStaleRatiosAreNotRestoredAsCurrent() {
        NodeRef node = makeStaleThenRestore();
        double[] ratios = delegate.getRatios();
        double[] expected = expectedRatios();
        assertTrue(Arrays.toString(ratios) + " != " + Arrays.toString(expected), equal(expected, ratios));
        // (a,b) at 1.6 under 2: ratio 0.8, where the stale vector had 0.5
        assertEquals(0.8, ratios[index(node)], 1E-12);
        double[] maskByDifference = delegate.setMaskByHeightDifference(0.5);

        // the same after a fresh computation (a tree-wide event)
        tree.pushTreeChangedEvent();
        assertTrue(equal(delegate.getRatios(), ratios));
        assertTrue(Arrays.equals(delegate.setMaskByHeightDifference(0.5), maskByDifference));
    }

    // the mask by ratio, read first after the restore, keeps the ratios strictly between 0.25 and 0.75 of the
    // restored heights: (a,b) is out at 0.8 (it was in at the stale 0.5)
    public void testMaskByRatioAfterRestore() {
        NodeRef node = makeStaleThenRestore();
        double[] mask = delegate.setMaskByRatio(0.25);
        double[] expected = expectedRatios();
        for (int i = 0; i < expected.length; i++) {
            assertEquals("ratio " + i + " = " + expected[i], expected[i] > 0.25 && expected[i] < 0.75 ? 1.0 : 0.0,
                    mask[i], 0.0);
        }
        assertEquals(0.0, mask[index(node)], 0.0);
    }

    // the ratios also come back right when they were current at the store
    public void testCurrentRatiosAreRestored() {
        double[] before = delegate.getRatios();
        tree.storeModelState();
        tree.setNodeHeight(parentOfTip("a"), 1.6);
        assertFalse(equal(before, delegate.getRatios()));
        tree.restoreModelState();
        assertTrue(equal(before, delegate.getRatios()));
        assertTrue(equal(expectedRatios(), delegate.getRatios()));
    }

    // -------------------------------------------------------------------------------------------------------------

    // read the ratios, move (a,b) from 1 to 1.6 (stale ratios), store, move d's parent, restore; returns (a,b)
    private NodeRef makeStaleThenRestore() {
        delegate.getRatios();
        NodeRef node = parentOfTip("a");
        tree.setNodeHeight(node, 1.6);
        tree.storeModelState();
        tree.setNodeHeight(parentOfTip("d"), 4.5);
        tree.restoreModelState();
        assertEquals(1.6, tree.getNodeHeight(node), 0.0);
        assertEquals(3.0, tree.getNodeHeight(parentOfTip("d")), 0.0);
        return node;
    }

    // node height over parent height, by ratio index (the anchors are all at 0)
    private double[] expectedRatios() {
        double[] expected = new double[delegate.getRatios().length];
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getInternalNode(i);
            if (!tree.isRoot(node)) {
                expected[index(node)] = tree.getNodeHeight(node) / tree.getNodeHeight(tree.getParent(node));
            }
        }
        return expected;
    }

    // the ratio index of a non-root internal node: internal nodes in number order, skipping the root
    private int index(NodeRef node) {
        int index = node.getNumber() - tree.getExternalNodeCount();
        return node.getNumber() > tree.getRoot().getNumber() ? index - 1 : index;
    }

    private NodeRef parentOfTip(String taxon) {
        return tree.getParent(tree.getExternalNode(tree.getTaxonIndex(taxon)));
    }

    private static boolean equal(double[] x, double[] y) {
        if (x.length != y.length) {
            return false;
        }
        for (int i = 0; i < x.length; i++) {
            if (Math.abs(x[i] - y[i]) > 1E-12) {
                return false;
            }
        }
        return true;
    }
}
