/*
 * ActionDataLikelihoodDelegate.java
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

    private final int partialsPerNode;
    private final int patternCount;
    private final int stateCount;
    private final int categoryCount;
    private double[] partials;
    private double[] weights;
    private double[] categoryRates;
    private double[] patternWeights;
    private DMatrixSparseCSC[] instantaneousMatrices;

    public ActionBeagleDelegate(int nodeCount,
                                int patternCount,
                                int stateCount,
                                int categoryCount,
                                int partialsPerNode,
                                DMatrixSparseCSC[] instantaneousMatrices) {
        this.partialsPerNode = partialsPerNode;
        this.patternCount = patternCount;
        this.stateCount = stateCount;
        this.categoryCount = categoryCount;
        this.weights = new double[categoryCount];
        this.categoryRates = new double[categoryCount];
        this.patternWeights = new double[patternCount];
        partials = new double[categoryCount * patternCount * stateCount * nodeCount * partialsPerNode];
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
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public void setTipStates(int i, int[] ints) {
        throw new RuntimeException("Should not be called with action likelihood");
    }

    @Override
    public void getTipStates(int i, int[] ints) {
        throw new RuntimeException("Should not be called with action likelihood");
    }

    @Override
    public void setTipPartials(int i, double[] doubles) {
        assert(doubles.length == categoryCount * patternCount * stateCount * partialsPerNode);
        System.arraycopy(doubles, 0, partials, getPartialIndex(i), doubles.length);
    }

    private int getPartialIndex(int i) {  // TODO: use BufferIndexHelper by XJ
        return i * categoryCount * patternCount * stateCount * partialsPerNode;
    }

    @Override
    public void setRootPrePartials(int[] ints, int[] ints1, int i) {
        throw new RuntimeException("Should not be called with action likelihood");
    }

    @Override
    public void setPartials(int i, double[] doubles) { //TODO: check for double buffering by XJ
        System.arraycopy(doubles, 0, partials, getPartialIndex(i), doubles.length);
    }

    @Override
    public void getPartials(int i, int i1, double[] doubles) {
        System.arraycopy(partials, getPartialIndex(i), doubles, 0, doubles.length);
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

    }

    @Override
    public void setCategoryWeights(int i, double[] doubles) {

    }

    @Override
    public void setCategoryRates(double[] doubles) {

    }

    @Override
    public void setCategoryRatesWithIndex(int i, double[] doubles) {

    }

    @Override
    public void convolveTransitionMatrices(int[] ints, int[] ints1, int[] ints2, int i) {

    }

    @Override
    public void addTransitionMatrices(int[] ints, int[] ints1, int[] ints2, int i) {

    }

    @Override
    public void transposeTransitionMatrices(int[] ints, int[] ints1, int i) {

    }

    @Override
    public void updateTransitionMatrices(int i, int[] ints, int[] ints1, int[] ints2, double[] doubles, int i1) {

    }

    @Override
    public void updateTransitionMatricesWithMultipleModels(int[] ints, int[] ints1, int[] ints2, int[] ints3, int[] ints4, double[] doubles, int i) {

    }

    @Override
    public void setTransitionMatrix(int i, double[] doubles, double v) {

    }

    @Override
    public void setDifferentialMatrix(int i, double[] doubles) {

    }

    @Override
    public void getTransitionMatrix(int i, double[] doubles) {

    }

    @Override
    public void updatePrePartials(int[] ints, int i, int i1) {

    }

    @Override
    public void updatePrePartialsByPartition(int[] ints, int i) {

    }

    @Override
    public void calculateEdgeDerivative(int[] ints, int[] ints1, int i, int[] ints2, int[] ints3, int i1, int i2, int i3, int[] ints4, int i4, double[] doubles, double[] doubles1) {

    }

    @Override
    public void calculateEdgeDifferentials(int[] ints, int[] ints1, int[] ints2, int[] ints3, int i, double[] doubles, double[] doubles1, double[] doubles2) {

    }

    @Override
    public void calculateCrossProductDifferentials(int[] ints, int[] ints1, int[] ints2, int[] ints3, double[] doubles, int i, double[] doubles1, double[] doubles2) {

    }

    @Override
    public void updatePartials(int[] ints, int i, int i1) {

    }

    @Override
    public void updatePartialsByPartition(int[] ints, int i) {

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
