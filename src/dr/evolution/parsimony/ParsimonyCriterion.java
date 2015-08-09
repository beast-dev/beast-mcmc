/*
 * ParsimonyCriterion.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evolution.parsimony;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;

/**
 * @author rambaut
 *         Date: Jun 20, 2005
 *         Time: 4:56:34 PM
 */
public interface ParsimonyCriterion {

    /**
     * Calculates the minimum number of steps for the parsimony reconstruction for the given tree.
     * It is expected that the implementation's constructor will be set up with the characters so
     * that repeated calls can be made to this function to evaluate different trees.
     * @param tree a tree object to reconstruct the characters on
     * @return an array containing the parsimony score for each site
     */
    double[] getSiteScores(Tree tree);

    /**
     * Calculates the minimum number of steps for the parsimony reconstruction for the given tree.
     * It is expected that the implementation's constructor will be set up with the characters so
     * that repeated calls can be made to this function to evaluate different trees.
     * @param tree a tree object to reconstruct the characters on
     * @return the total score
     */
    double getScore(Tree tree);


    /**
     * Returns the reconstructed character states for a given node in the tree.
     * @param tree a tree object to reconstruct the characters on
     * @param node the node of the tree
     * @return an array containing the reconstructed states for this node
     */
    int[] getStates(Tree tree, NodeRef node);
}
