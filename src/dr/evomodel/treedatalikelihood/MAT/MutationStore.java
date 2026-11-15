/*
 * MutationStore.java
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
 * Read access to a mutation annotated tree: the state of every pattern at every node, and the
 * mutations (including ambiguities) this implies on each branch. Node numbers refer to the tree
 * the store was built for, and the store records that tree's topology as it was when the states
 * were last made consistent with it.
 */
public interface MutationStore {

    int getNodeCount();

    int getPatternCount();

    /**
     * @return the parent of the node in the topology the store currently describes, -1 for the root
     */
    int getParent(int nodeNumber);

    int getRoot();

    /**
     * @return the state of the pattern at the node, or -1 if it is ambiguous in every tip below the node
     */
    int getState(int nodeNumber, int pattern);

    /**
     * @return the substitutions and ambiguities on the branch above the node
     */
    BranchMutations getBranchMutations(int nodeNumber);

    /**
     * @return the substitutions and ambiguities on the branch above the node as Mutation objects
     * (a convenience; use {@link #getBranchMutations} in tight loops)
     */
    default Mutation[] getMutationsForNode(int nodeNumber) {
        return getBranchMutations(nodeNumber).toMutations(nodeNumber);
    }

    /**
     * @return the total number of substitutions (excluding ambiguities)
     */
    int getMutationCount();

    /**
     * @return the total number of ambiguities
     */
    int getAmbiguityCount();
}
