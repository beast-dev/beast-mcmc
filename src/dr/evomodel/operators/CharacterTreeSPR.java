/*
 * CharacterTreeSPR.java
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

package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.substmodel.CharacterTreeModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * A root-preserving subtree prune and regraft on the binary skeleton of a character tree: i uniform among the nodes
 * whose parent p is not the root; p is pruned with i and regrafted on the edge above j, uniform among the nodes outside
 * i's subtree other than p, i's sibling and the root. The number of choices is the same before and after the move, so
 * the log Hastings ratio is 0. Node heights are reset to unit depth in the same edit.
 *
 * @author Xiang Ji
 */
public class CharacterTreeSPR extends SimpleMCMCOperator {

    public static final String CHARACTER_TREE_SPR = "characterTreeSPR";

    private final CharacterTreeModel characterTree;

    public CharacterTreeSPR(CharacterTreeModel characterTree, double weight) {
        if (!characterTree.isVariable()) {
            throw new IllegalArgumentException("the character tree " + characterTree.getId() + " is fixed");
        }
        if (characterTree.getExternalNodeCount() < 4) {
            throw new IllegalArgumentException("a subtree prune and regraft needs at least 4 states");
        }
        this.characterTree = characterTree;
        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return CHARACTER_TREE_SPR + "(" + characterTree.getId() + ")";
    }

    private static void mark(DefaultTreeModel tree, NodeRef node, boolean[] inside) {
        inside[node.getNumber()] = true;
        for (int c = 0; c < tree.getChildCount(node); c++) {
            mark(tree, tree.getChild(node, c), inside);
        }
    }

    @Override
    public double doOperation() {
        final DefaultTreeModel skeleton = characterTree.getSkeleton();
        final NodeRef root = skeleton.getRoot();
        final List<NodeRef> candidates = new ArrayList<NodeRef>();
        for (int n = 0; n < skeleton.getNodeCount(); n++) {
            final NodeRef node = skeleton.getNode(n);
            if (node != root && skeleton.getParent(node) != root) {
                candidates.add(node);
            }
        }
        final NodeRef i = candidates.get(MathUtils.nextInt(candidates.size()));
        final NodeRef p = skeleton.getParent(i);
        final NodeRef grandparent = skeleton.getParent(p);
        final NodeRef sibling = skeleton.getChild(p, 0) == i ? skeleton.getChild(p, 1) : skeleton.getChild(p, 0);
        final boolean[] excluded = new boolean[skeleton.getNodeCount()];
        mark(skeleton, i, excluded);
        excluded[p.getNumber()] = true;
        excluded[sibling.getNumber()] = true;
        excluded[root.getNumber()] = true;
        final List<NodeRef> targets = new ArrayList<NodeRef>();
        for (int n = 0; n < skeleton.getNodeCount(); n++) {
            if (!excluded[n]) {
                targets.add(skeleton.getNode(n));
            }
        }
        final NodeRef j = targets.get(MathUtils.nextInt(targets.size()));
        skeleton.beginTreeEdit();
        skeleton.removeChild(p, sibling);
        skeleton.removeChild(grandparent, p);
        skeleton.addChild(grandparent, sibling);
        final NodeRef jParent = skeleton.getParent(j);
        skeleton.removeChild(jParent, j);
        skeleton.addChild(jParent, p);
        skeleton.addChild(p, j);
        characterTree.resetUnitHeightsQuietly();
        skeleton.endTreeEdit();
        return 0.0;
    }
}
