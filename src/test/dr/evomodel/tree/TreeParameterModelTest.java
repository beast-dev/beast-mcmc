/*
 * TreeParameterModelTest.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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

package test.dr.evomodel.tree;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.FastTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import junit.framework.TestCase;

/**
 * A TreeParameterModel without the root keeps one value per branch. When a subtree move changes the
 * root, the values must stay with their nodes, and the old root, now an ordinary node, takes the new
 * root's former value (the relaxed clocks rely on this to keep their rates on their branches).
 */
public class TreeParameterModelTest extends TestCase {

    private static final String NEWICK = "((((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,G:3):1):1,H:5);";

    public void testValuesStayWithTheirNodesAcrossRootMoves() throws Exception {
        checkRootMoves(new DefaultTreeModel("treeModel", new NewickImporter(NEWICK).importTree(null), false, false));
        checkRootMoves(new FastTreeModel("treeModel", new NewickImporter(NEWICK).importTree(null), false, false));
    }

    private void checkRootMoves(TreeModel tree) {
        MathUtils.setSeed(7);
        final int nodeCount = tree.getNodeCount();
        TreeParameterModel values = new TreeParameterModel(tree, new Parameter.Default(nodeCount - 1), false);
        // the value each node should hold, by node number (the root's entry is unused)
        double[] expected = new double[nodeCount];
        for (int n = 0; n < nodeCount; n++) {
            NodeRef node = tree.getNode(n);
            if (!tree.isRoot(node)) {
                values.setNodeValue(tree, node, 100 + n);
                expected[n] = 100 + n;
            }
        }

        int moves = 0;
        int rootMovesDown = 0;
        int rootMovesUp = 0;
        for (int iteration = 0; iteration < 5000 && moves < 300; iteration++) {
            NodeRef x = tree.getNode(MathUtils.nextInt(nodeCount));
            NodeRef destination = tree.getNode(MathUtils.nextInt(nodeCount));
            if (tree.isRoot(x) || isAncestorOrSelf(tree, x, destination)) {
                continue;
            }
            NodeRef parent = tree.getParent(x);
            NodeRef sibling = otherChild(tree, parent, x);
            if (destination == parent || destination == sibling) {
                continue;
            }
            NodeRef destinationParent = tree.getParent(destination);
            double top = destinationParent == null ? Double.POSITIVE_INFINITY : tree.getNodeHeight(destinationParent);
            double bottom = Math.max(tree.getNodeHeight(x), tree.getNodeHeight(destination));
            if (top <= bottom) {
                continue;
            }
            double newHeight = destinationParent == null ? bottom + 1.0 + MathUtils.nextDouble() :
                    bottom + MathUtils.nextDouble() * (top - bottom);

            final int oldRoot = tree.getRoot().getNumber();
            subtreePruneRegraft(tree, x, destination, newHeight);
            final int newRoot = tree.getRoot().getNumber();
            if (newRoot != oldRoot) {
                expected[oldRoot] = expected[newRoot];
                if (oldRoot > newRoot) {
                    rootMovesDown++;
                } else {
                    rootMovesUp++;
                }
            }
            for (int n = 0; n < nodeCount; n++) {
                NodeRef node = tree.getNode(n);
                if (!tree.isRoot(node)) {
                    assertEquals(tree.getClass().getSimpleName() + ": node " + n + " after move " + moves +
                            " (root " + oldRoot + " -> " + newRoot + ")", expected[n], values.getNodeValue(tree, node), 0.0);
                }
            }
            moves++;
        }
        assertEquals(300, moves);
        assertTrue("root moves down: " + rootMovesDown, rootMovesDown > 10);
        assertTrue("root moves up: " + rootMovesUp, rootMovesUp > 10);
    }

    private static void subtreePruneRegraft(TreeModel tree, NodeRef x, NodeRef destination, double height) {
        NodeRef parent = tree.getParent(x);
        NodeRef sibling = otherChild(tree, parent, x);
        NodeRef grandParent = tree.getParent(parent);
        NodeRef destinationParent = tree.getParent(destination);
        tree.beginTreeEdit();
        if (grandParent == null) {
            tree.removeChild(parent, sibling);
            tree.setRoot(sibling);
        } else {
            tree.removeChild(parent, sibling);
            tree.removeChild(grandParent, parent);
            tree.addChild(grandParent, sibling);
        }
        if (destinationParent == null) {
            tree.addChild(parent, destination);
            tree.setRoot(parent);
        } else {
            tree.removeChild(destinationParent, destination);
            tree.addChild(parent, destination);
            tree.addChild(destinationParent, parent);
        }
        tree.setNodeHeight(parent, height);
        tree.endTreeEdit();
    }

    private static boolean isAncestorOrSelf(TreeModel tree, NodeRef ancestor, NodeRef node) {
        for (NodeRef n = node; n != null; n = tree.getParent(n)) {
            if (n == ancestor) {
                return true;
            }
        }
        return false;
    }

    private static NodeRef otherChild(TreeModel tree, NodeRef parent, NodeRef child) {
        return tree.getChild(parent, 0) == child ? tree.getChild(parent, 1) : tree.getChild(parent, 0);
    }
}
