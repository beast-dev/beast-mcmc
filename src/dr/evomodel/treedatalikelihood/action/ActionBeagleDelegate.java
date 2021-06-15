/*
 * ActionBeagleDelegate.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.action;

import beagle.Beagle;
import beagle.InstanceDetails;
import org.newejml.data.DMatrixSparseCSC;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class ActionBeagleDelegate implements Beagle {

    public static final boolean DEBUG = false;
    public static final boolean SCALING = true;
    public static final int SCALING_FACTOR_COUNT = 254;
    public static final int SCALING_FACTOR_OFFSET = 126;
    private static final int SCALING_EXPONENT_THRESHOLD = 2;
    protected final int tipCount;
    protected final int partialsBufferCount;
    protected final int stateCount;
    protected final int patternCount;
    protected final int matrixBufferCount;
    protected final int categoryCount;
    protected int partialsSize;
    protected int matrixSize;
    protected double[][] cMatrices;
    protected double[][] stateFrequencies;
    protected double[] categoryRates;
    protected double[] categoryWeights;
    protected double[] patternWeights;
    protected double[][] partials;
    protected int[][] scalingFactorCounts;
    protected double[][] matrices;
    double[] tmpPartials;
    protected double[] scalingFactors;
    protected double[] logScalingFactors;
    protected DMatrixSparseCSC[] instantaneousMatrices;

    public ActionBeagleDelegate(int tipCount,
                                int partialsBufferCount,
                                int patternCount,
                                int stateCount,
                                int categoryCount,
                                int matrixBufferCount,
                                int partialsSize,
                                DMatrixSparseCSC[] instantaneousMatrices) {
        this.tipCount = tipCount;
        this.partialsBufferCount = partialsBufferCount;
        this.patternCount = patternCount;
        this.stateCount = stateCount;
        this.categoryCount = categoryCount;
        this.matrixBufferCount = matrixBufferCount;
        this.partialsSize = partialsSize;
        this.categoryWeights = new double[categoryCount];
        this.categoryRates = new double[categoryCount];
        this.patternWeights = new double[patternCount];
        partials = new double[partialsBufferCount][];
        for (int i = 0; i < partialsBufferCount; i++) {
            partials[i] = new double[partialsSize];
        }
        this.instantaneousMatrices = instantaneousMatrices;
    }

    @Override
    public void finalize() throws Throwable {

    }

    @Override
    public void setCPUThreadCount(int i) {

    }

    @Override
    public void setPatternWeights(double[] doubles) {
        assert(doubles.length == patternCount);
        System.arraycopy(doubles, 0, patternWeights, 0, doubles.length);
    }

    @Override
    public void setPatternPartitions(int i, int[] ints) {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    @Override
    public void setTipStates(int i, int[] ints) {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    @Override
    public void getTipStates(int i, int[] ints) {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    @Override
    public void setTipPartials(int i, double[] doubles) {
        assert(i >= 0 && i < tipCount);
        assert(doubles.length == partialsSize);
        if (this.partials[i] == null) {
            this.partials[i] = new double[partialsSize];
        }

        int partialIndex = 0;
        for (int category = 0; category < categoryCount; category++) {
            System.arraycopy(doubles, 0, partials[i], partialIndex, doubles.length);
            partialIndex += doubles.length;
        }
    }
//
//    private int getPartialIndex(int i) {  // TODO: use BufferIndexHelper by XJ
//        return i * categoryCount * patternCount * stateCount * partialsPerNode;
//    }

    @Override
    public void setRootPrePartials(int[] ints, int[] ints1, int i) {
        throw new RuntimeException("Should not be called with action likelihood");
    }

    @Override
    public void setPartials(int i, double[] doubles) { //TODO: check for double buffering by XJ
        assert this.partials[i] != null;

        System.arraycopy(doubles, 0, this.partials[i], 0, this.partialsSize);
    }

    @Override
    public void getPartials(int i, int i1, double[] doubles) {
        System.arraycopy(this.partials[i], 0, doubles, 0, this.partialsSize);
    }

    @Override
    public void getLogScaleFactors(int i, double[] doubles) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public void setEigenDecomposition(int i, double[] doubles, double[] doubles1, double[] doubles2) {
        throw new RuntimeException("Should not be called with action likelihood");
    }

    @Override
    public void setStateFrequencies(int i, double[] doubles) {
        System.arraycopy(doubles, 0, this.stateFrequencies[i], 0, this.stateCount);
    }

    @Override
    public void setCategoryWeights(int i, double[] doubles) {
        System.arraycopy(doubles, 0, this.categoryWeights[i], 0, this.categoryCount);
    }

    @Override
    public void setCategoryRates(double[] doubles) {
        System.arraycopy(doubles, 0, this.categoryRates, 0, this.categoryRates.length);
    }

    @Override
    public void setCategoryRatesWithIndex(int i, double[] doubles) {
        throw new UnsupportedOperationException("setCategoryRatesWithIndex not implemented in ActionBeagleDelegate.");
    }

    @Override
    public void convolveTransitionMatrices(int[] ints, int[] ints1, int[] ints2, int i) {
        throw new UnsupportedOperationException("convolveTransitionMatrices not implemented in ActionBeagleDelegate.");
    }

    @Override
    public void addTransitionMatrices(int[] ints, int[] ints1, int[] ints2, int i) {
        throw new UnsupportedOperationException("addTransitionMatrices not implemented in ActionBeagleDelegate.");
    }

    @Override
    public void transposeTransitionMatrices(int[] ints, int[] ints1, int i) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void updateTransitionMatrices(int i, int[] ints, int[] ints1, int[] ints2, double[] doubles, int i1) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void updateTransitionMatricesWithMultipleModels(int[] ints, int[] ints1, int[] ints2, int[] ints3, int[] ints4, double[] doubles, int i) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void setTransitionMatrix(int i, double[] doubles, double v) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void setDifferentialMatrix(int i, double[] doubles) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void getTransitionMatrix(int i, double[] doubles) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void updatePrePartials(int[] ints, int i, int i1) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void updatePrePartialsByPartition(int[] ints, int i) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void calculateEdgeDerivative(int[] ints, int[] ints1, int i, int[] ints2, int[] ints3, int i1, int i2, int i3, int[] ints4, int i4, double[] doubles, double[] doubles1) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void calculateEdgeDifferentials(int[] ints, int[] ints1, int[] ints2, int[] ints3, int i, double[] doubles, double[] doubles1, double[] doubles2) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void calculateCrossProductDifferentials(int[] ints, int[] ints1, int[] ints2, int[] ints3, double[] doubles, int i, double[] doubles1, double[] doubles2) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void updatePartials(int[] ints, int i, int i1) {

    }

    @Override
    public void updatePartialsByPartition(int[] ints, int i) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void accumulateScaleFactors(int[] ints, int i, int i1) {

    }

    @Override
    public void accumulateScaleFactorsByPartition(int[] ints, int i, int i1, int i2) {

    }

    @Override
    public void removeScaleFactors(int[] ints, int i, int i1) {

    }

    @Override
    public void removeScaleFactorsByPartition(int[] ints, int i, int i1, int i2) {

    }

    @Override
    public void copyScaleFactors(int i, int i1) {

    }

    @Override
    public void resetScaleFactors(int i) {

    }

    @Override
    public void resetScaleFactorsByPartition(int i, int i1) {

    }

    @Override
    public void calculateRootLogLikelihoods(int[] ints, int[] ints1, int[] ints2, int[] ints3, int i, double[] doubles) {

    }

    @Override
    public void calculateRootLogLikelihoodsByPartition(int[] ints, int[] ints1, int[] ints2, int[] ints3, int[] ints4, int i, int i1, double[] doubles, double[] doubles1) {

    }

    @Override
    public void getSiteLogLikelihoods(double[] doubles) {

    }

    @Override
    public InstanceDetails getDetails() {
        return null;
    }
}
