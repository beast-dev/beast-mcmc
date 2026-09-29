/*
 * EpochEvolutionaryProcessDelegate.java
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

package dr.evomodel.treedatalikelihood;

/**
 * An evolutionary process that changes at a fixed sequence of epoch transition times. Instead of convolving the
 * transition matrices of the epochs that a branch spans, a traversal adds a degree-2 node to the branch at each
 * transition time, so that every branch segment lies within a single epoch and uses a single substitution model.
 *
 * @author Marc A Suchard
 */
public interface EpochEvolutionaryProcessDelegate extends EvolutionaryProcessDelegate {

    /**
     * @return the heights at which the epoch changes, in strictly ascending order. Epoch e lies between
     * transition time e - 1 and e, so there is one more epoch (and substitution model) than transition times.
     * The array must not be modified.
     */
    double[] getEpochTransitionTimes();

    /**
     * @return the bookkeeping of the degree-2 nodes on the branches, which also identifies the epoch of the
     * transition matrix of every original and augmented node
     */
    AugmentedNodeRegistry getAugmentedNodeRegistry();

    /**
     * @return the index of a transition matrix that is the identity, which lets a degree-2 node be computed
     * as a node with two children where the second child contributes nothing. The matrix is set by
     * the first call of updateTransitionMatrices().
     */
    int getIdentityMatrixIndex();
}
