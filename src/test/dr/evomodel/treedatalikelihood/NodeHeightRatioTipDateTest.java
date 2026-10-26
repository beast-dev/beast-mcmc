/*
 * NodeHeightRatioTipDateTest.java
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
import dr.inference.model.Parameter;
import dr.util.Transform;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.StringReader;
import java.util.Map;

/**
 * The node-height ratio transform with sampled tip dates: each internal node's ratio is measured from the oldest tip
 * below it (its epoch's anchor), so when tip heights change order the epochs must follow, and a rejected change must
 * bring back the epochs of the restored tree. Small ratios put every node just above its anchor, so a stale anchor
 * shows up as a node below one of its tips.
 */
public class NodeHeightRatioTipDateTest extends TestCase {

    // heights: (a,b) 2, ((a,b),c) 3, root 4; every tip at 0 until setUp moves b to 1
    private static final String XML = "<beast>\n" +
            "  <taxa id=\"taxa\"><taxon id=\"a\"/><taxon id=\"b\"/><taxon id=\"c\"/><taxon id=\"d\"/></taxa>\n" +
            "  <newick id=\"startingTree\">(((a:2,b:2):1,c:3):1,d:4);</newick>\n" +
            "  <treeModel id=\"treeModel\">\n" +
            "    <tree idref=\"startingTree\"/>\n" +
            "    <rootHeight><parameter id=\"treeModel.rootHeight\"/></rootHeight>\n" +
            "    <nodeHeights internalNodes=\"true\"><parameter id=\"treeModel.internalNodeHeights\"/></nodeHeights>\n" +
            "    <nodeHeights internalNodes=\"true\" rootNode=\"true\">" +
            "<parameter id=\"treeModel.allInternalNodeHeights\"/></nodeHeights>\n" +
            "    <leafHeight taxon=\"a\"><parameter id=\"age(a)\"/></leafHeight>\n" +
            "    <leafHeight taxon=\"b\"><parameter id=\"age(b)\"/></leafHeight>\n" +
            "  </treeModel>\n" +
            "  <strictClockBranchRates id=\"branchRates\"><rate><parameter id=\"clock.rate\" value=\"1.0\"/></rate>" +
            "</strictClockBranchRates>\n" +
            "  <nodeHeightProxyParameter id=\"nodeHeights\" rootNode=\"true\"><treeModel idref=\"treeModel\"/>" +
            "</nodeHeightProxyParameter>\n" +
            "  <nodeHeightTransform id=\"nodeHeightTransform\" withRoot=\"true\">\n" +
            "    <treeModel idref=\"treeModel\"/><strictClockBranchRates idref=\"branchRates\"/>\n" +
            "    <ratios><parameter id=\"ratios\"/></ratios>\n" +
            "  </nodeHeightTransform>\n" +
            "</beast>\n";

    private TreeModel tree;
    private Parameter nodeHeights;
    private Transform transform;
    private Parameter ageA;
    private Parameter ageB;

    public NodeHeightRatioTipDateTest(String name) {
        super(name);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        Map<String, XMLObject> objects = parser.parse(new StringReader(XML), false);
        tree = (TreeModel) objects.get("treeModel").getNativeObject();
        nodeHeights = (Parameter) objects.get("nodeHeights").getNativeObject();
        transform = (Transform) objects.get("nodeHeightTransform").getNativeObject();
        ageA = (Parameter) objects.get("age(a)").getNativeObject();
        ageB = (Parameter) objects.get("age(b)").getNativeObject();
        ageB.setParameterValue(0, 1.0);
        assertEquals(2.0, tree.getNodeHeight(tree.getParent(tree.getExternalNode(1))), 1E-12);
        // build the epochs for the starting tips: b anchors (a,b)
        assertSmallRatiosGiveAValidTree();
    }

    // a becomes the older of the two tips below (a,b): the epochs must be rebuilt with a as the anchor
    public void testTipsChangeOrder() {
        ageB.setParameterValue(0, 0.0);
        ageA.setParameterValue(0, 1.5);
        assertSmallRatiosGiveAValidTree();
    }

    // the same change, rejected: after the restore b anchors (a,b) again
    public void testRejectedTipChange() {
        tree.storeModelState();
        ageB.setParameterValue(0, 0.0);
        ageA.setParameterValue(0, 1.5);
        assertSmallRatiosGiveAValidTree();
        tree.restoreModelState();
        assertEquals(1.0, ageB.getParameterValue(0), 0.0);
        assertEquals(0.0, ageA.getParameterValue(0), 0.0);
        assertSmallRatiosGiveAValidTree();
    }

    // the restore must not keep the epochs built for the rejected tips, even when inverse() comes next, with no
    // transform() (which sets the heights and so rebuilds the epochs anyway) in between
    public void testRejectedTipChangeThenInverse() {
        tree.storeModelState();
        ageB.setParameterValue(0, 0.0);
        ageA.setParameterValue(0, 1.5);
        assertSmallRatiosGiveAValidTree();
        tree.restoreModelState();
        final int dim = nodeHeights.getDimension();
        double[] y = new double[dim];
        y[0] = tree.getNodeHeight(tree.getRoot());
        for (int i = 1; i < dim; i++) {
            y[i] = 0.01;
        }
        transform.inverse(y, 0, dim);
        assertValidTree();
    }

    // every ratio near 0 puts each internal node just above its epoch's anchor tip
    private void assertSmallRatiosGiveAValidTree() {
        final int dim = nodeHeights.getDimension();
        double[] y = transform.transform(nodeHeights.getParameterValues(), 0, dim);
        for (int i = 1; i < dim; i++) { // y[0] is the root height
            y[i] = 0.01;
        }
        transform.inverse(y, 0, dim);
        assertValidTree();
    }

    private void assertValidTree() {
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getInternalNode(i);
            for (int c = 0; c < tree.getChildCount(node); c++) {
                NodeRef child = tree.getChild(node, c);
                assertTrue("node " + node.getNumber() + " at " + tree.getNodeHeight(node) + " is below its child " +
                                child.getNumber() + " at " + tree.getNodeHeight(child),
                        tree.getNodeHeight(node) > tree.getNodeHeight(child));
            }
        }
    }
}
