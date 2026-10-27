/*
 * NniNeighbourhood.java
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

package dr.evomodel.tree.search;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * The rooted NNIs of a time tree, 2 (n - 2) of them, in the search's priority order (plan 7 §5.8): nodes whose
 * branch is short relative to its parent's height first, (h(G) - h(P)) / h(G) increasing, computed from the current
 * heights; ties go to the lower node number, and at each node the move of P's first child comes first. A branch
 * collapsing towards zero is where a topology change is most likely to help.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class NniNeighbourhood {

    private NniNeighbourhood() {
    }

    /**
     * @return the moves, in priority order, leaving out those the constraints forbid
     */
    public static List<RootedNniMove> enumerate(final TreeModel tree, TopologyConstraints constraints) {
        List<NodeRef> nodes = new ArrayList<NodeRef>();
        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
            NodeRef node = tree.getInternalNode(i);
            if (!tree.isRoot(node)) {
                nodes.add(node);
            }
        }
        final double[] priority = new double[tree.getNodeCount()];
        for (NodeRef node : nodes) {
            final double parent = tree.getNodeHeight(tree.getParent(node));
            priority[node.getNumber()] = parent > 0.0 ? (parent - tree.getNodeHeight(node)) / parent : 0.0;
        }
        Collections.sort(nodes, new Comparator<NodeRef>() {
            public int compare(NodeRef x, NodeRef y) {
                int c = Double.compare(priority[x.getNumber()], priority[y.getNumber()]);
                return c != 0 ? c : Integer.compare(x.getNumber(), y.getNumber());
            }
        });

        List<RootedNniMove> moves = new ArrayList<RootedNniMove>();
        for (NodeRef p : nodes) {
            final NodeRef g = tree.getParent(p);
            final NodeRef u = tree.getChild(g, 0) == p ? tree.getChild(g, 1) : tree.getChild(g, 0);
            for (int c = 0; c < tree.getChildCount(p); ++c) {
                RootedNniMove move = new RootedNniMove(p.getNumber(), g.getNumber(),
                        tree.getChild(p, c).getNumber(), u.getNumber());
                if (constraints == null || !constraints.violates(tree, move)) {
                    moves.add(move);
                }
            }
        }
        return moves;
    }
}
