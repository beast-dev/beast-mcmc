/*
 * OldAbstractLikelihoodCore.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.treelikelihood;

import java.util.logging.Logger;

/**
 * AbstractLikelihoodCore - An abstract base class for LikelihoodCores
 *
 * @version $Id: OldAbstractLikelihoodCore.java,v 1.1 2006/08/31 14:57:24 rambaut Exp $
 *
 * @author Andrew Rambaut
 */

public abstract class OldAbstractLikelihoodCore implements LikelihoodCore {

    protected int stateCount;
    protected int nodeCount;
    protected int patternCount;
    protected int partialsSize;
    protected int matrixSize;
    protected int matrixCount;

    protected boolean integrateCategories;

    protected double[][] partials;
    protected double[][] storedPartials;

    protected int[][] states;

    protected double[][] matrices;
    protected double[][] storedMatrices;

    protected boolean useScaling;

    protected double[][] scalingFactors;
    protected double[][] storedScalingFactors;

    private final double scalingThreshold = 1.0E-100;
    private int scalingCheckCount;

    /**
     * Constructor
     *
     * @param stateCount number of states
     */
    public OldAbstractLikelihoodCore(int stateCount) {
        this.stateCount = stateCount;
    }

    /**
     * initializes partial likelihood arrays.
     * @param nodeCount the number of nodes in the tree
     * @param patternCount the number of patterns
     * @param matrixCount the number of matrices (i.e., number of categories)
     * @param integrateCategories whether sites are being integrated over all matrices
     */
    public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories, boolean useScaling) {

        this.nodeCount = nodeCount;
        this.patternCount = patternCount;
        this.matrixCount = matrixCount;

        this.integrateCategories = integrateCategories;

        if (integrateCategories) {
            partialsSize = patternCount * stateCount * matrixCount;
        } else {
            partialsSize = patternCount * stateCount;
        }

        partials = new double[nodeCount][];
        storedPartials = new double[nodeCount][];

        states = new int[nodeCount][];

        for (int i = 0; i < nodeCount; i++) {
            partials[i] = null;
            storedPartials[i] = null;

            states[i] = null;
        }

        matrixSize = stateCount * stateCount;

        matrices = new double[nodeCount][matrixCount * matrixSize];
        storedMatrices = new double[nodeCount][matrixCount * matrixSize];

        this.useScaling = useScaling;

        if (useScaling) {
            scalingFactors = new double[nodeCount][patternCount];
            storedScalingFactors = new double[nodeCount][patternCount];

            scalingCheckCount = 0;
        }
    }

    /**
     * cleans up and deallocates arrays.
     */
    public void finalize() {

        nodeCount = 0;
        patternCount = 0;
        matrixCount = 0;

        partials = null;
        storedPartials = null;
        states = null;
        matrices = null;
        storedMatrices = null;

        scalingFactors = null;
        storedScalingFactors = null;
    }

    /**
     * Allocates partials for a node
     */
    public void createNodePartials(int nodeIndex) {

        this.partials[nodeIndex] = new double[partialsSize];
        this.storedPartials[nodeIndex] = new double[partialsSize];
    }

    /**
     * Sets partials for a node
     */
    public void setNodePartials(int nodeIndex, double[] partials) {

        if (this.partials[nodeIndex] == null) {
            createNodePartials(nodeIndex);
        }
        if (partials.length < partialsSize) {
            int k = 0;
            for (int i = 0; i < matrixCount; i++)  {
                System.arraycopy(partials, 0, this.partials[nodeIndex], k, partials.length);
                k += partials.length;
            }
        } else {
            System.arraycopy(partials, 0, this.partials[nodeIndex], 0, partials.length);
        }
    }

    /**
     * Allocates states for a node
     */
    public void createNodeStates(int nodeIndex) {

        this.states[nodeIndex] = new int[patternCount];
    }

    /**
     * Sets states for a node
     */
    public void setNodeStates(int nodeIndex, int[] states) {

        if (this.states[nodeIndex] == null) {
            createNodeStates(nodeIndex);
        }
        System.arraycopy(states, 0, this.states[nodeIndex], 0, patternCount);
    }

    public void setNodeMatrixForUpdate(int nodeIndex) {
        // nothing to do
    }

    /**
     * Sets probability matrix for a node
     */
    public void setNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {

        System.arraycopy(matrix, 0, matrices[nodeIndex],
                matrixIndex * matrixSize, matrixSize);
    }


    public void setNodePartialsForUpdate(int nodeIndex) {
        // nothing to do
    }

    /**
     * Calculates partial likelihoods at a node.
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     */
    public void calculatePartials( int nodeIndex1, int nodeIndex2, int nodeIndex3 )
    {
        if (states[nodeIndex1] != null) {
            if (states[nodeIndex2] != null) {
                calculateStatesStatesPruning(states[nodeIndex1], matrices[nodeIndex1],
                        states[nodeIndex2], matrices[nodeIndex2],
                        partials[nodeIndex3]);
            } else {
                calculateStatesPartialsPruning(states[nodeIndex1], matrices[nodeIndex1],
                        partials[nodeIndex2], matrices[nodeIndex2],
                        partials[nodeIndex3]);
            }
        } else {
            if (states[nodeIndex2] != null) {
                calculateStatesPartialsPruning(states[nodeIndex2], matrices[nodeIndex2],
                        partials[nodeIndex1], matrices[nodeIndex1],
                        partials[nodeIndex3]);
            } else {
                calculatePartialsPartialsPruning(partials[nodeIndex1], matrices[nodeIndex1],
                        partials[nodeIndex2], matrices[nodeIndex2],
                        partials[nodeIndex3]);
            }
        }

        if (useScaling) {
            scalePartials(nodeIndex3);
        }
    }

    /**
     * Calculates partial likelihoods at a node when both children have states.
     */
    protected abstract void calculateStatesStatesPruning(int[] states1, double[] matrices1,
                                                         int[] states2, double[] matrices2,
                                                         double[] partials3);

    /**
     * Calculates partial likelihoods at a node when one child has states and one has partials.
     */
    protected abstract void calculateStatesPartialsPruning(	int[] states1, double[] matrices1,
                                                               double[] partials2, double[] matrices2,
                                                               double[] partials3);

    /**
     * Calculates partial likelihoods at a node when both children have partials.
     */
    protected abstract void calculatePartialsPartialsPruning(double[] partials1, double[] matrices1,
                                                             double[] partials2, double[] matrices2,
                                                             double[] partials3);

    /**
     * Calculates partial likelihoods at a node.
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     * @param matrixMap a map of which matrix to use for each pattern (can be null if integrating over categories)
     */
    public void calculatePartials( int nodeIndex1, int nodeIndex2, int nodeIndex3, int[] matrixMap ) {
        if (states[nodeIndex1] != null) {
            if (states[nodeIndex2] != null) {
                calculateStatesStatesPruning(states[nodeIndex1], matrices[nodeIndex1],
                        states[nodeIndex2], matrices[nodeIndex2],
                        partials[nodeIndex3], matrixMap);
            } else {
                calculateStatesPartialsPruning(states[nodeIndex1], matrices[nodeIndex1],
                        partials[nodeIndex2], matrices[nodeIndex2],
                        partials[nodeIndex3], matrixMap);
            }
        } else {
            if (states[nodeIndex2] != null) {
                calculateStatesPartialsPruning(states[nodeIndex2], matrices[nodeIndex2],
                        partials[nodeIndex1], matrices[nodeIndex1],
                        partials[nodeIndex3], matrixMap);
            } else {
                calculatePartialsPartialsPruning(partials[nodeIndex1], matrices[nodeIndex1],
                        partials[nodeIndex2], matrices[nodeIndex2],
                        partials[nodeIndex3], matrixMap);
            }
        }

        if (useScaling) {
            scalePartials(nodeIndex3);
        }
    }

    /**
     * Calculates partial likelihoods at a node when both children have states.
     */
    protected abstract void calculateStatesStatesPruning(int[] states1, double[] matrices1,
                                                         int[] states2, double[] matrices2,
                                                         double[] partials3, int[] matrixMap);

    /**
     * Calculates partial likelihoods at a node when one child has states and one has partials.
     */
    protected abstract void calculateStatesPartialsPruning(	int[] states1, double[] matrices1,
                                                               double[] partials2, double[] matrices2,
                                                               double[] partials3, int[] matrixMap);

    /**
     * Calculates partial likelihoods at a node when both children have partials.
     */
    protected abstract void calculatePartialsPartialsPruning(double[] partials1, double[] matrices1,
                                                             double[] partials2, double[] matrices2,
                                                             double[] partials3, int[] matrixMap);


    public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials)
    {
        calculateIntegratePartials(partials[nodeIndex], proportions, outPartials);
    }

    /**
     * Integrates partials across categories.
     * @param inPartials the partials at the node to be integrated
     * @param proportions the proportions of sites in each category
     * @param outPartials an array into which the integrated partials will go
     */
    protected abstract void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials);

    /**
     * Scale the partials at a given node. This uses a scaling suggested by Ziheng Yang in
     * Yang (2000) J. Mol. Evol. 51: 423-432
     *
     * This function looks over the partial likelihoods for each state at each pattern
     * and finds the largest. If this is less than the scalingThreshold (currently set
     * to 1E-40) then it rescales the partials for that pattern by dividing by this number
     * (i.e., normalizing to between 0, 1). It then stores the log of this scaling.
     * This is called for every internal node after the partials are calculated so provides
     * most of the performance hit. Ziheng suggests only doing this on a proportion of nodes
     * but this sounded like a headache to organize (and he doesn't use the threshold idea
     * which improves the performance quite a bit).
     * @param nodeIndex
     */
    protected void scalePartials(int nodeIndex) {
        int u = 0;

        for (int i = 0; i < patternCount; i++) {

            double scaleFactor = 0.0;
            int v = u;
            for (int k = 0; k < matrixCount; k++) {
                for (int j = 0; j < stateCount; j++) {
                    if (partials[nodeIndex][v] > scaleFactor) {
                        scaleFactor = partials[nodeIndex][v];
                    }
                    v++;
                }
                v += patternCount;
            }

            if (scaleFactor < scalingThreshold) {

                v = u;
                for (int k = 0; k < matrixCount; k++) {
                    for (int j = 0; j < stateCount; j++) {
                        partials[nodeIndex][v] /= scaleFactor;
                        v++;
                    }
                    v += patternCount;
                }
                scalingFactors[nodeIndex][i] = Math.log(scaleFactor);
            } else {
                scalingFactors[nodeIndex][i] = 0.0;
            }
            u += stateCount;
        }
    }

    /**
     * This function returns the scaling factor for that pattern by summing over
     * the log scalings used at each node. If scaling is off then this just returns
     * a 0.
     * @return the log scaling factor
     */
    protected double getLogScalingFactor(int pattern) {
        double logScalingFactor = 0.0;
        if (useScaling) {
            for (int i = 0; i < scalingFactors.length; i++) {
                logScalingFactor += scalingFactors[i][pattern];
            }
        }
        return logScalingFactor;
    }

    /**
     * Check whether the scaling is still required. If the sum of all the logScalingFactors
     * is zero then we simply turn off the useScaling flag. This will speed up the likelihood
     * calculations when scaling is not required.
     */
    public void checkScaling() {
        if (useScaling) {
            if (scalingCheckCount % 1000 == 0) {
                double totalScalingFactor = 0.0;
                for (int i = 0; i < nodeCount; i++) {
                    for (int j = 0; j < patternCount; j++) {
                        totalScalingFactor += scalingFactors[i][j];
                    }
                }
                useScaling = totalScalingFactor < 0.0;
                Logger.getLogger("dr.evomodel").info("LikelihoodCore total log scaling factor: " + totalScalingFactor);
                if (!useScaling) {
                    Logger.getLogger("dr.evomodel").info("LikelihoodCore scaling turned off.");
                }
            }
            scalingCheckCount++;
        }
    }

    /**
     * Gets the partials for a particular node.
     * @param nodeIndex the node
     * @param outPartials an array into which the partials will go
     */
    public void getPartials(int nodeIndex, double[] outPartials)
    {
        double[] partials1 = partials[nodeIndex];

        for (int k = 0; k < partialsSize; k++) {
            outPartials[k] = partials1[k];
        }
    }

    /**
     * Store current state
     */
    public void storeState() {

        for (int i = 0; i < nodeCount; i++) {
            if (partials[i] != null) {
                System.arraycopy(partials[i], 0, storedPartials[i], 0, partialsSize);
    }
            System.arraycopy(matrices[i], 0, storedMatrices[i], 0, matrixSize * matrixCount);

            if (useScaling) {
                System.arraycopy(scalingFactors[i], 0, storedScalingFactors[i], 0, patternCount);
            }
        }

    }

    /**
     * Restore the stored state
     */
    public void restoreState() {
        // Rather than copying the stored stuff back, just swap the pointers...
        double[][] tmp1 = partials;
        partials = storedPartials;
        storedPartials = tmp1;

        double[][] tmp2 = matrices;
        matrices = storedMatrices;
        storedMatrices = tmp2;

        if (useScaling) {
            double[][] tmp3 = scalingFactors;
            scalingFactors = storedScalingFactors;
            storedScalingFactors = tmp3;
    }
}
}
