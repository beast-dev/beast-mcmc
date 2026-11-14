/*
 * CharacterTreeNNI.java
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

/**
 * A rooted nearest-neighbour interchange on the binary skeleton of a character tree: u uniform among the non-root
 * internal nodes, one of its two children exchanged with its sibling. Node heights are reset to unit depth in the
 * same edit, since the character tree reads only its topology. The move is symmetric (log Hastings ratio 0) and
 * connects all binary trees; at a collapsed u the contracted tree does not change.
 *
 * @author Xiang Ji
 */
public class CharacterTreeNNI extends SimpleMCMCOperator {

    public static final String CHARACTER_TREE_NNI = "characterTreeNNI";

    private final CharacterTreeModel characterTree;
    private final int[] nonRootInternal;

    public CharacterTreeNNI(CharacterTreeModel characterTree, double weight) {
        if (!characterTree.isVariable()) {
            throw new IllegalArgumentException("the character tree " + characterTree.getId() + " is fixed");
        }
        this.characterTree = characterTree;
        final DefaultTreeModel skeleton = characterTree.getSkeleton();
        final int tips = skeleton.getExternalNodeCount();
        final int root = skeleton.getRoot().getNumber();
        this.nonRootInternal = new int[tips - 2];
        int k = 0;
        for (int n = tips; n < skeleton.getNodeCount(); n++) {
            if (n != root) {
                nonRootInternal[k++] = n;
            }
        }
        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return CHARACTER_TREE_NNI + "(" + characterTree.getId() + ")";
    }

    @Override
    public double doOperation() {
        final DefaultTreeModel skeleton = characterTree.getSkeleton();
        final NodeRef u = skeleton.getNode(nonRootInternal[MathUtils.nextInt(nonRootInternal.length)]);
        final NodeRef parent = skeleton.getParent(u);
        final NodeRef sibling = skeleton.getChild(parent, 0) == u ? skeleton.getChild(parent, 1) : skeleton.getChild(parent, 0);
        final NodeRef child = skeleton.getChild(u, MathUtils.nextInt(2));
        skeleton.beginTreeEdit();
        skeleton.removeChild(u, child);
        skeleton.removeChild(parent, sibling);
        skeleton.addChild(u, sibling);
        skeleton.addChild(parent, child);
        characterTree.resetUnitHeightsQuietly();
        skeleton.endTreeEdit();
        return 0.0;
    }
}
