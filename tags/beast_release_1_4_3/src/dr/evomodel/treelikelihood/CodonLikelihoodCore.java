/*
 * CodonLikelihoodCore.java
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
 * CodonLikelihoodCore - An implementation of LikelihoodCore for codon data
 *
 * @version $Id: CodonLikelihoodCore.java,v 1.13 2006/08/29 18:07:23 rambaut Exp $
 *
 * @author Andrew Rambaut
 */

public class CodonLikelihoodCore implements LikelihoodCore {

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

    private final double scalingThreshold = 1.0E-30;
    private int scalingCheckCount;

    /**
     * Constructor
     *
     * @param stateCount number of states
     */
    public CodonLikelihoodCore(int stateCount) {
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

        if (!integrateCategories) {
            throw new RuntimeException("Integrate Categories must be used with CodonLikelihoodCore");
        }

        partialsSize = patternCount * stateCount * matrixCount;

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
        throw new RuntimeException("createNodeStates not implemented for CodonLikelihoodCore");
    }

    /**
     * Specify that the matrices for the given node are about to be updated
     *
     * @param nodeIndex
     */
    public void setNodeMatrixForUpdate(int nodeIndex) {
        // nothing to do
    }

    /**
     * Sets states for a node
     */
    public void setNodeStates(int nodeIndex, int[] states) {
        throw new RuntimeException("createNodeStates not implemented for CodonLikelihoodCore");
    }

    /**
     * Sets probability matrix for a node
     */
    public void setNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {

        System.arraycopy(matrix, 0, matrices[nodeIndex],
                                matrixIndex * matrixSize, matrixSize);
    }

    /**
     * Specify that the partials for the given node are about to be updated
     *
     * @param nodeIndex
     */
    public void setNodePartialsForUpdate(int nodeIndex) {
        // nothing to do
    }

    /**
     * Scale the partials at a given node. This uses a scaling suggested by Ziheng Yang in
     * Yang (2000) J. Mol. Evol. 51: 423-432
     * @param nodeIndex
     */
    protected void scalePartials(int nodeIndex) {
        int k = 0;
        int l = 0;

        for (int i = 0; i < patternCount; i++) {
            double scaleFactor = partials[nodeIndex][k];
            k++;
            for (int j = 1; j < stateCount; j++) {
                if (partials[nodeIndex][k] > scaleFactor) {
                    scaleFactor = partials[nodeIndex][k];
                }
                k++;
            }

            if (scaleFactor < scalingThreshold) {

                for (int j = 0; j < stateCount; j++) {
                    partials[nodeIndex][l] /= scaleFactor;
                    l++;
                }
                scalingFactors[nodeIndex][i] = Math.log(scaleFactor);
            } else {
                scalingFactors[nodeIndex][i] = 0.0;
                l += stateCount;
            }
        }
    }

    /**
     * Get the per pattern log scaling factor. This should be added to the log likelihood of
     * each pattern. This will generally be non-zero when a likelihood core is scaling the
     * partial likelihoods in order to avoid underflow numerical errors.
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
     * Calculates partial likelihoods at a node.
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     */
    public void calculatePartials( int nodeIndex1, int nodeIndex2, int nodeIndex3 )
    {
        double[] partials1 = partials[nodeIndex1];
        double[] matrices1 = matrices[nodeIndex1];
        double[] partials2 = partials[nodeIndex2];
        double[] matrices2 = matrices[nodeIndex2];
        double[] partials3 = partials[nodeIndex3];

        double sum1, sum2;

        int u = 0;
        int v = 0;

        for (int l = 0; l < matrixCount; l++) {

            for (int k = 0; k < patternCount; k++) {

                int w = l * matrixCount;

                for (int i = 0; i < stateCount; i++) {

                    sum1 = sum2 = 0.0;

                    for (int j = 0; j < stateCount; j++) {
                        sum1 += matrices1[w] * partials1[v + j];
                        sum2 += matrices2[w] * partials2[v + j];

                        w++;
                    }

                    partials3[u] = sum1 * sum2;
                    u++;
                }
                v += stateCount;
            }
        }


        if (useScaling) {
            scalePartials(nodeIndex3);
        }
    }

    /**
     * Calculates partial likelihoods at a node.
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     * @param matrixMap a map of which matrix to use for each pattern (can be null if integrating over categories)
     */
    public void calculatePartials( int nodeIndex1, int nodeIndex2, int nodeIndex3, int[] matrixMap )
    {
        throw new RuntimeException("calculatePartials not implemented for CodonLikelihoodCore");
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
     * Integrates partials across categories.
     * @param nodeIndex the node
     * @param proportions the proportions of sites in each category
     * @param outPartials an array into which the partials will go
     */
    public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials)
    {
        double[] partials1 = partials[nodeIndex];

        int u = 0;
        int v = 0;
        for (int k = 0; k < patternCount; k++) {

            for (int i = 0; i < stateCount; i++) {

                outPartials[u] = partials1[v] * proportions[0];
                u++;
                v++;
            }
        }


        for (int l = 1; l < matrixCount; l++) {
            u = 0;

            for (int k = 0; k < patternCount; k++) {

                for (int i = 0; i < stateCount; i++) {

                    outPartials[u] += partials1[v] * proportions[l];
                    u++;
                    v++;
                }
            }
        }
    }

    /**
     * Calculates site likelihoods at a node.
     * @param partials the partials used to calculate the likelihoods
     * @param frequencies an array of state frequencies
     * @param outLogLikelihoods an array into which the likelihoods will go
     */
    public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods)
    {
        int v = 0;
        for (int k = 0; k < patternCount; k++) {

            double logScalingFactor = getLogScalingFactor(k);

            double sum = 0.0;
            for (int i = 0; i < stateCount; i++) {

                sum += frequencies[i] * partials[v];
                v++;
            }
            outLogLikelihoods[k] = Math.log(sum) + logScalingFactor;
        }
        checkScaling();
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
