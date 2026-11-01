/*
 * CheckpointTreeAdoptionTest.java
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
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import org.junit.Test;

import java.io.StringReader;
import java.util.ArrayList;

import static org.junit.Assert.*;

/**
 * A run resumed from a checkpoint must end up with the same node numbering, tree and
 * node-indexed parameter values as the run that wrote it, whatever starting tree the
 * resumed run began with.
 */
public class CheckpointTreeAdoptionTest {

    private static DefaultTreeModel tree(String newick) throws Exception {
        return new DefaultTreeModel("treeModel", new NewickImporter(new StringReader(newick)).importTree(null));
    }

    @Test
    public void tipsAreNumberedByTaxonNameWhateverTheTopology() throws Exception {
        DefaultTreeModel tree1 = tree("((a:1,b:1):1,(c:1,d:1):1);");
        DefaultTreeModel tree2 = tree("((d:1,a:1):1,(b:1,c:1):1);");
        DefaultTreeModel tree3 = tree("(d:3,(c:2,(b:1,a:1):1):1);");

        String[] expected = {"a", "b", "c", "d"};
        for (DefaultTreeModel tree : new DefaultTreeModel[]{tree1, tree2, tree3}) {
            for (int i = 0; i < expected.length; i++) {
                assertEquals(expected[i], tree.getNodeTaxon(tree.getExternalNode(i)).getId());
                assertEquals(i, tree.getExternalNode(i).getNumber());
            }
        }
    }

    /**
     * Moves the root of the tree so that it is no longer the last node: ((a,b),(c,d)) becomes (a,(b,(c,d)))
     * with the old root's node now inside the tree.
     */
    private static void moveRoot(DefaultTreeModel tree) {
        NodeRef root = tree.getRoot();
        NodeRef u = tree.getChild(root, 0);        // (a,b)
        NodeRef a = tree.getChild(u, 0);
        NodeRef b = tree.getChild(u, 1);

        tree.beginTreeEdit();
        tree.removeChild(root, u);
        tree.removeChild(u, a);
        tree.removeChild(u, b);
        tree.addChild(root, b);                   // root = (b, (c,d))
        tree.addChild(u, a);
        tree.addChild(u, root);                   // u = (a, root)
        tree.setNodeHeight(u, tree.getNodeHeight(root) + 1.0);
        tree.setRoot(u);
        tree.endTreeEdit();
    }

    @Test
    public void adoptedTreeAndTraitsMatchTheSavedRunWhenTheRootHasMoved() throws Exception {
        // the run that writes the checkpoint: its root is not the last node any more
        DefaultTreeModel saved = tree("((a:1,b:1):1,(c:1,d:1):1);");
        moveRoot(saved);
        assertTrue("the test needs a root that is not the last node",
                saved.getRoot().getNumber() != saved.getNodeCount() - 1);

        Parameter savedCategories = new Parameter.Default("categories", saved.getNodeCount() - 1, 0.0);
        TreeParameterModel savedModel = new TreeParameterModel(saved, savedCategories, false);
        for (int i = 0; i < saved.getNodeCount(); i++) {
            NodeRef node = saved.getNode(i);
            if (!saved.isRoot(node)) {
                savedModel.setNodeValue(saved, node, 10.0 * i + 1);
            }
        }

        // what the checkpointer writes: node heights and taxa by node number, edges with the trait values
        int nodeCount = saved.getNodeCount();
        double[] heights = new double[nodeCount];
        String[] taxaNames = new String[saved.getExternalNodeCount()];
        int[] edges = new int[nodeCount];
        int[] childOrder = new int[nodeCount];
        double[][] traitValues = new double[1][nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            NodeRef node = saved.getNode(i);
            heights[i] = saved.getNodeHeight(node);
            if (saved.isExternal(node)) {
                taxaNames[i] = saved.getNodeTaxon(node).getId();
            }
            NodeRef parent = saved.getParent(node);
            if (parent == null) {
                edges[i] = -1;
                childOrder[i] = -1;
            } else {
                edges[i] = parent.getNumber();
                childOrder[i] = (saved.getChild(parent, 0) == node ? 0 : 1);
                traitValues[0][i] = savedModel.getNodeValue(saved, node);
            }
        }

        // the resumed run: a different (random) starting tree and the same model structure
        DefaultTreeModel resumed = tree("(d:3,(c:2,(b:1,a:1):1):1);");
        Parameter resumedCategories = new Parameter.Default("categories", resumed.getNodeCount() - 1, 0.0);
        TreeParameterModel resumedModel = new TreeParameterModel(resumed, resumedCategories, false);

        // the checkpointer restores every parameter by index, including the model's root node number
        for (int k = 0; k < savedCategories.getDimension(); k++) {
            resumedCategories.setParameterValue(k, savedCategories.getParameterValue(k));
        }
        boolean rootNodeNumberRestored = false;
        for (int v = 0; v < resumedModel.getVariableCount(); v++) {
            Variable variable = resumedModel.getVariable(v);
            String name = ((Parameter) variable).getParameterName();
            if (name != null && name.endsWith(".rootNodeNumber")) {
                ((Parameter) variable).setParameterValue(0, saved.getRoot().getNumber());
                rootNodeNumberRestored = true;
            }
        }
        assertTrue("the model's root node number should be a named variable", rootNodeNumberRestored);

        ArrayList<TreeParameterModel> traitModels = new ArrayList<TreeParameterModel>();
        traitModels.add(resumedModel);

        resumed.beginTreeEdit();
        resumed.adoptTreeStructure(edges, heights, childOrder, taxaNames);
        resumed.adoptTraitData(edges, traitModels, traitValues, taxaNames);
        resumed.endTreeEdit();

        assertEquals(saved.getRoot().getNumber(), resumed.getRoot().getNumber());
        assertEquals(TreeUtils.uniqueNewick(saved, saved.getRoot()), TreeUtils.uniqueNewick(resumed, resumed.getRoot()));
        for (int i = 0; i < nodeCount; i++) {
            NodeRef savedNode = saved.getNode(i);
            NodeRef resumedNode = resumed.getNode(i);
            assertEquals(saved.getNodeHeight(savedNode), resumed.getNodeHeight(resumedNode), 0.0);
            if (saved.isExternal(savedNode)) {
                assertEquals(saved.getNodeTaxon(savedNode).getId(), resumed.getNodeTaxon(resumedNode).getId());
            }
            if (!saved.isRoot(savedNode)) {
                assertEquals("value on node " + i, savedModel.getNodeValue(saved, savedNode),
                        resumedModel.getNodeValue(resumed, resumedNode), 0.0);
            }
        }
        // and the parameter itself is laid out identically, so operators that pick an index behave the same
        for (int k = 0; k < savedCategories.getDimension(); k++) {
            assertEquals("parameter index " + k, savedCategories.getParameterValue(k), resumedCategories.getParameterValue(k), 0.0);
        }
    }
}
