/*
 * ArbitraryPrecisionLikelihoodCore.java
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

import java.math.BigDecimal;
import java.math.MathContext;

/**
 * GeneralLikelihoodCore - An implementation of LikelihoodCore for any data
 *
 * @author Andrew Rambaut
 * @version $Id: GeneralLikelihoodCore.java,v 1.28 2006/08/31 14:57:24 rambaut Exp $
 */

@Deprecated // Switching to BEAGLE
public class ArbitraryPrecisionLikelihoodCore implements LikelihoodCore {

    private MathContext precision;

    private int stateCount;
    private int nodeCount;
    private int patternCount;
    private int partialsSize;
    private int matrixSize;
    private int matrixCount;

    private boolean integrateCategories;

    private BigDecimal[][][] partials;

    private int[][] states;

    private BigDecimal[][][] matrices;

    private int[] currentMatricesIndices;
    private int[] storedMatricesIndices;
    private int[] currentPartialsIndices;
    private int[] storedPartialsIndices;


    /**
     * Constructor
     *
     * @param stateCount number of states
     */
    public ArbitraryPrecisionLikelihoodCore(int stateCount, int precision) {
        this.stateCount = stateCount;
        this.precision = new MathContext(precision);
    }

    /**
     * initializes partial likelihood arrays.
     *
     * @param nodeCount           the number of nodes in the tree
     * @param patternCount        the number of patterns
     * @param matrixCount         the number of matrices (i.e., number of categories)
     * @param integrateCategories whether sites are being integrated over all matrices
     */
    public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories) {

        this.nodeCount = nodeCount;
        this.patternCount = patternCount;
        this.matrixCount = matrixCount;

        this.integrateCategories = integrateCategories;

        if (integrateCategories) {
            partialsSize = patternCount * stateCount * matrixCount;
        } else {
            partialsSize = patternCount * stateCount;
        }

        partials = new BigDecimal[2][nodeCount][];

        currentMatricesIndices = new int[nodeCount];
        storedMatricesIndices = new int[nodeCount];

        currentPartialsIndices = new int[nodeCount];
        storedPartialsIndices = new int[nodeCount];

        states = new int[nodeCount][];

        for (int i = 0; i < nodeCount; i++) {
            partials[0][i] = null;
            partials[1][i] = null;

            states[i] = null;
        }

        matrixSize = stateCount * stateCount;

        matrices = new BigDecimal[2][nodeCount][matrixCount * matrixSize];
    }

    /**
     * cleans up and deallocates arrays.
     */
    public void finalize() throws java.lang.Throwable {
        super.finalize();

        nodeCount = 0;
        patternCount = 0;
        matrixCount = 0;

        partials = null;
        currentPartialsIndices = null;
        storedPartialsIndices = null;
        states = null;
        matrices = null;
        currentMatricesIndices = null;
        storedMatricesIndices = null;
    }

    /**
     * Allocates partials for a node
     */
    public void createNodePartials(int nodeIndex) {

        this.partials[0][nodeIndex] = new BigDecimal[partialsSize];
        this.partials[1][nodeIndex] = new BigDecimal[partialsSize];
    }

    /**
     * Sets partials for a node
     */
    public void setNodePartials(int nodeIndex, double[] partials) {

        if (this.partials[0][nodeIndex] == null) {
            createNodePartials(nodeIndex);
        }
        if (partials.length < partialsSize) {
            int k = 0;
            for (int i = 0; i < matrixCount; i++) {
                for (int j = 0; j < partials.length; j++) {
                    this.partials[0][nodeIndex][k] = new BigDecimal(partials[j], precision);
                    k++;
                }
            }
        } else {
            for (int j = 0; j < partials.length; j++) {
                this.partials[0][nodeIndex][j] = new BigDecimal(partials[j], precision);
            }
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

    /**
     * Gets states for a node
     */
    public void getNodeStates(int nodeIndex, int[] states) {
        System.arraycopy(this.states[nodeIndex], 0, states, 0, patternCount);
    }

    public void setNodeMatrixForUpdate(int nodeIndex) {
        currentMatricesIndices[nodeIndex] = 1 - currentMatricesIndices[nodeIndex];

    }

    /**
     * Sets probability matrix for a node
     */
    public void setNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {

        for (int j = 0; j < matrixSize; j++) {
            matrices[currentMatricesIndices[nodeIndex]][nodeIndex][(matrixIndex * matrixSize) + j] = new BigDecimal(matrix[j], precision);
        }
    }

    /**
     * Gets probability matrix for a node
     */
    public void getNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {
        for (int j = 0; j < matrixSize; j++) {
            matrix[j] = matrices[currentMatricesIndices[nodeIndex]][nodeIndex][(matrixIndex * matrixSize) + j].doubleValue();
        }
    }

    public void setNodePartialsForUpdate(int nodeIndex) {
        currentPartialsIndices[nodeIndex] = 1 - currentPartialsIndices[nodeIndex];
    }

    /**
     * Get the currently updating node partials for direct access
     *
     * @param nodeIndex
     */
    public void setCurrentNodePartials(int nodeIndex, double[] partials) {
        throw new UnsupportedOperationException("setCurrentNodePartials is not supported by ArbitraryPrecisionLikelihoodCore");
    }

    /**
     * Calculates partial likelihoods at a node.
     *
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     */
    public void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3) {
        if (states[nodeIndex1] != null) {
            if (states[nodeIndex2] != null) {
                calculateStatesStatesPruning(
                        states[nodeIndex1], matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
                        states[nodeIndex2], matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndices[nodeIndex3]][nodeIndex3]);
            } else {
                calculateStatesPartialsPruning(states[nodeIndex1], matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
                        partials[currentPartialsIndices[nodeIndex2]][nodeIndex2], matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndices[nodeIndex3]][nodeIndex3]);
            }
        } else {
            if (states[nodeIndex2] != null) {
                calculateStatesPartialsPruning(states[nodeIndex2], matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndices[nodeIndex1]][nodeIndex1], matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
                        partials[currentPartialsIndices[nodeIndex3]][nodeIndex3]);
            } else {
                calculatePartialsPartialsPruning(partials[currentPartialsIndices[nodeIndex1]][nodeIndex1], matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
                        partials[currentPartialsIndices[nodeIndex2]][nodeIndex2], matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndices[nodeIndex3]][nodeIndex3]);
            }
        }
    }

    /**
     * Calculates partial likelihoods at a node when both children have states.
     */
    private void calculateStatesStatesPruning(int[] states1, BigDecimal[] matrices1,
                                              int[] states2, BigDecimal[] matrices2,
                                              BigDecimal[] partials3) {
        int v = 0;

        for (int l = 0; l < matrixCount; l++) {

            for (int k = 0; k < patternCount; k++) {

                int state1 = states1[k];
                int state2 = states2[k];

                int w = l * matrixSize;

                if (state1 < stateCount && state2 < stateCount) {

                    for (int i = 0; i < stateCount; i++) {

                        partials3[v] = matrices1[w + state1].multiply(matrices2[w + state2], precision);

                        v++;
                        w += stateCount;
                    }

                } else if (state1 < stateCount) {
                    // child 2 has a gap or unknown state so treat it as unknown

                    for (int i = 0; i < stateCount; i++) {

                        partials3[v] = matrices1[w + state1];

                        v++;
                        w += stateCount;
                    }
                } else if (state2 < stateCount) {
                    // child 2 has a gap or unknown state so treat it as unknown

                    for (int i = 0; i < stateCount; i++) {

                        partials3[v] = matrices2[w + state2];

                        v++;
                        w += stateCount;
                    }
                } else {
                    // both children have a gap or unknown state so set partials to 1

                    for (int j = 0; j < stateCount; j++) {
                        partials3[v] = BigDecimal.ONE;
                        v++;
                    }
                }
            }
        }
    }

    /**
     * Calculates partial likelihoods at a node when one child has states and one has partials.
     */
    private void calculateStatesPartialsPruning(int[] states1, BigDecimal[] matrices1,
                                                BigDecimal[] partials2, BigDecimal[] matrices2,
                                                BigDecimal[] partials3) {

        BigDecimal sum, tmp;

        int u = 0;
        int v = 0;

        for (int l = 0; l < matrixCount; l++) {
            for (int k = 0; k < patternCount; k++) {

                int state1 = states1[k];

                int w = l * matrixSize;

                if (state1 < stateCount) {


                    for (int i = 0; i < stateCount; i++) {

                        tmp = matrices1[w + state1];

                        sum = BigDecimal.ZERO;
                        for (int j = 0; j < stateCount; j++) {
                            sum = sum.add(matrices2[w].multiply(partials2[v + j], precision), precision);
                            w++;
                        }

                        partials3[u] = tmp.multiply(sum, precision);
                        u++;
                    }

                    v += stateCount;
                } else {
                    // Child 1 has a gap or unknown state so don't use it

                    for (int i = 0; i < stateCount; i++) {

                        sum = BigDecimal.ZERO;
                        for (int j = 0; j < stateCount; j++) {
                            sum = sum.add(matrices2[w].multiply(partials2[v + j], precision), precision);
                            w++;
                        }

                        partials3[u] = sum;
                        u++;
                    }

                    v += stateCount;
                }
            }
        }
    }

    /**
     * Calculates partial likelihoods at a node when both children have partials.
     */
    private void calculatePartialsPartialsPruning(BigDecimal[] partials1, BigDecimal[] matrices1,
                                                  BigDecimal[] partials2, BigDecimal[] matrices2,
                                                  BigDecimal[] partials3) {
        BigDecimal sum1, sum2;

        int u = 0;
        int v = 0;

        for (int l = 0; l < matrixCount; l++) {

            for (int k = 0; k < patternCount; k++) {

                int w = l * matrixSize;

                for (int i = 0; i < stateCount; i++) {

                    sum1 = sum2 = BigDecimal.ZERO;

                    for (int j = 0; j < stateCount; j++) {
                        sum1 = sum1.add(matrices1[w].multiply(partials1[v + j], precision), precision);
                        sum2 = sum2.add(matrices2[w].multiply(partials2[v + j], precision), precision);

                        w++;
                    }

                    partials3[u] = sum1.multiply(sum2, precision);
                    u++;
                }
                v += stateCount;
            }
        }
    }

    /**
     * Calculates partial likelihoods at a node.
     *
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     * @param matrixMap  a map of which matrix to use for each pattern (can be null if integrating over categories)
     */
    public void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3, int[] matrixMap) {
        throw new UnsupportedOperationException("calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3, int[] matrixMap) is not implemented in this likelihood core");
    }


    public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {
        throw new UnsupportedOperationException("integratePartials(int nodeIndex, double[] proportions, double[] outPartials) is not implemented in this likelihood core");
    }

    public void integratePartials(int nodeIndex, double[] proportions, BigDecimal[] outPartials) {
        BigDecimal[] prop = new BigDecimal[proportions.length];
        for (int i = 0; i < proportions.length; i++) {
            prop[i] = new BigDecimal(proportions[i], precision);
        }
        calculateIntegratePartials(partials[currentPartialsIndices[nodeIndex]][nodeIndex], prop, outPartials);
    }

    /**
     * Integrates partials across categories.
     *
     * @param inPartials  the array of partials to be integrated
     * @param proportions the proportions of sites in each category
     * @param outPartials an array into which the partials will go
     */
    private void calculateIntegratePartials(BigDecimal[] inPartials, BigDecimal[] proportions, BigDecimal[] outPartials) {

        int u = 0;
        int v = 0;
        for (int k = 0; k < patternCount; k++) {

            for (int i = 0; i < stateCount; i++) {

                outPartials[u] = inPartials[v].multiply(proportions[0], precision);
                u++;
                v++;
            }
        }


        for (int l = 1; l < matrixCount; l++) {
            u = 0;

            for (int k = 0; k < patternCount; k++) {

                for (int i = 0; i < stateCount; i++) {

                    outPartials[u] = outPartials[u].add(inPartials[v].multiply(proportions[l], precision), precision);
                    u++;
                    v++;
                }
            }
        }
    }

    /**
     * Calculates patten log likelihoods at a node.
     *
     * @param partials          the partials used to calculate the likelihoods
     * @param frequencies       an array of state frequencies
     * @param outLogLikelihoods an array into which the likelihoods will go
     */
    public void calculateLogLikelihoods(BigDecimal[] partials, double[] frequencies, double[] outLogLikelihoods) {
        BigDecimal[] freqs = new BigDecimal[frequencies.length];
        for (int i = 0; i < freqs.length; i++) {
            freqs[i] = new BigDecimal(frequencies[i], precision);
        }

        int v = 0;
        for (int k = 0; k < patternCount; k++) {

            BigDecimal sum = BigDecimal.ZERO;
            for (int i = 0; i < stateCount; i++) {

                sum = sum.add(freqs[i].multiply(partials[v], precision), precision);
                v++;
            }

            double scale = sum.scale();
            double value = sum.unscaledValue().doubleValue();

            outLogLikelihoods[k] = Math.log(value) - (scale * Math.log(10));
        }
    }

    public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) {
        throw new UnsupportedOperationException("calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) is not implemented in this likelihood core");
    }

    public void setUseScaling(boolean useScaling) {
    }

    public double getLogScalingFactor(int pattern) {
        return 0;
    }

    public boolean arePartialsRescaled() {
        return false;
    }

    public void getLogScalingFactors(int nodeIndex, double[] buffer) {
        throw new RuntimeException("Not yet implemented.");
    }

    /**
     * Gets the partials for a particular node.
     *
     * @param nodeIndex   the node
     * @param outPartials an array into which the partials will go
     */
    public void getPartials(int nodeIndex, double[] outPartials) {
        for (int i = 0; i < outPartials.length; i++) {
            outPartials[i] = partials[currentPartialsIndices[nodeIndex]][nodeIndex][i].doubleValue();
        }
    }

    /**
     * Store current state
     */
    public void storeState() {

        System.arraycopy(currentMatricesIndices, 0, storedMatricesIndices, 0, nodeCount);
        System.arraycopy(currentPartialsIndices, 0, storedPartialsIndices, 0, nodeCount);
    }

    /**
     * Restore the stored state
     */
    public void restoreState() {
        // Rather than copying the stored stuff back, just swap the pointers...
        int[] tmp1 = currentMatricesIndices;
        currentMatricesIndices = storedMatricesIndices;
        storedMatricesIndices = tmp1;

        int[] tmp2 = currentPartialsIndices;
        currentPartialsIndices = storedPartialsIndices;
        storedPartialsIndices = tmp2;
    }

    public void checkScaling() {
        // do nothing
    }

}