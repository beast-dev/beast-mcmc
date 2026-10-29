/*
 * TopologyView.java
 *
 * Copyright © 2026, the BEAST Development Team.
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
 */

package dr.evomodel.treedatalikelihood.MAT;

/**
 * A read-only topology by node number, with the length (in expected substitutions) of the branch above
 * each node: either the tree as it is, or the tree as it was before a subtree move (so that the states
 * under the old topology can be scored after the tree has changed).
 */
interface TopologyView {

    /** @return the node's parent, or -1 for the root */
    int getParent(int node);

    int getChildCount(int node);

    int getChild(int node, int i);

    boolean isExternal(int node);

    int getRoot();

    /** @return the length of the branch above the node in expected substitutions (0 for the root) */
    double getBranchLength(int node);
}
