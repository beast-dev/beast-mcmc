/*
 * LikelihoodCore.java
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

package dr.oldevomodel.treelikelihood;

import dr.evomodel.treelikelihood.LikelihoodPartialsProvider;

/**
 * LikelihoodCore - An interface describing the core likelihood functions.
 *
 * @author Andrew Rambaut
 * @version $Id: LikelihoodCore.java,v 1.15 2006/08/29 18:07:23 rambaut Exp $
 */

@Deprecated // Switching to BEAGLE
public interface LikelihoodCore extends LikelihoodPartialsProvider {

    /**
     * initializes partial likelihood arrays.
     */
    void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories);

    /**
     * cleans up and deallocates arrays.
     */
    void finalize() throws java.lang.Throwable;

    /**
     * Allocates partials for a node
     */
    void createNodePartials(int nodeIndex);

    /**
     * Sets partials for a node
     */
    void setNodePartials(int nodeIndex, double[] partials);

    /**
     * Sets states for a node
     */
    void setNodeStates(int nodeIndex, int[] states);

    /**
     * Sets states for a node and a pattern
     * Allocates states for a node
     */
    void createNodeStates(int nodeIndex);

    /**
     * Specify that the matrices for the given node are about to be updated
     *
     * @param nodeIndex
     */
    void setNodeMatrixForUpdate(int nodeIndex);

    /**
     * Sets probability matrix for a node
     */
    void setNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix);

    /**
     * Specify that the partials for the given node are about to be updated
     *
     * @param nodeIndex
     */
    void setNodePartialsForUpdate(int nodeIndex);

    /**
     * Set the currently updating node partials from the array...
     *
     * @param nodeIndex
     */
    void setCurrentNodePartials(int nodeIndex, double[] partials);

    /**
     * Calculates partial likelihoods at a node using a matrixMap.
     *
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     */
    void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3);

    /**
     * Calculates partial likelihoods at a node using a matrixMap.
     *
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     * @param matrixMap  a map of which matrix to use for each pattern
     */
    void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3, int[] matrixMap);

    /**
     * Gets the partials for a particular node.
     *
     * @param nodeIndex   the node
     * @param outPartials an array into which the partials will go
     */
    public void getPartials(int nodeIndex, double[] outPartials);

    /**
     * Integrates partials across categories.
     *
     * @param nodeIndex   the node at which to calculate the likelihoods
     * @param proportions the proportions of sites in each category
     * @param outPartials an array into which the partials will go
     */
    void integratePartials(int nodeIndex, double[] proportions, double[] outPartials);

    /**
     * Calculates pattern log likelihoods at a node.
     *
     * @param partials          the partials used to calculate the likelihoods
     * @param frequencies       an array of state frequencies
     * @param outLogLikelihoods an array into which the log likelihoods will go
     */
    void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods);

    void setUseScaling(boolean useScaling);

    double getLogScalingFactor(int pattern);

    boolean arePartialsRescaled();

    void getLogScalingFactors(int nodeIndex, double[] buffer);

    /**
     * Store current state
     */
    void storeState();

    /**
     * Restore the stored state
     */
    void restoreState();
}
